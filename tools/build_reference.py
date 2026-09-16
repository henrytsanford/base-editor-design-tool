"""Builds a local Ensembl reference bundle from the release FTP files.

The design script needs five things per transcript, all of them pure functions
of (transcript ID, Ensembl release). This script distills them out of the bulk
release files once, so a design run needs no network at all:

    <out>/ensembl-<release>/transcripts.db     transcript records, exon and CDS
                                               blocks, CDS and protein sequences
    <out>/ensembl-<release>/genome.fa.bgz      primary assembly, bgzip-compressed
                          + .fai / .gzi        indexes for random access
    <out>/ensembl-<release>/manifest.json      release, source URLs, build date

Run once per release:

    python tools/build_reference.py --release 116

It is safe to re-run; each stage is skipped if its output is already present
and newer than its input. Expect ~30-60 minutes and ~5 GB of transient disk on
the first run.
"""
import argparse
import gzip
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import date, datetime

from Bio import bgzf
from Bio.SeqIO.FastaIO import SimpleFastaParser

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from transcript_source import GENOME_FA, TRANSCRIPTS_DB

FTP_BASE = 'https://ftp.ensembl.org/pub'
SPECIES = 'homo_sapiens'
ASSEMBLY = 'GRCh38'

# Soft-masked (dna_sm) and hard-masked (dna_rm) assemblies would pass the
# script's sequence check but silently break base-edit detection, which
# compares against uppercase 'C'/'A'. Only the plain assembly is acceptable.
DNA_FLAVOUR = 'dna.primary_assembly'

# ftp.ensembl.org throttles each connection; a handful in parallel is much faster.
DEFAULT_STREAMS = 8


def source_urls(release):
    """The four release files the bundle is built from."""
    base = '%s/release-%s' % (FTP_BASE, release)
    fa = '%s/fasta/%s' % (base, SPECIES)
    return {
        'gtf': '%s/gtf/%s/Homo_sapiens.%s.%s.gtf.gz' % (base, SPECIES, ASSEMBLY, release),
        'dna': '%s/dna/Homo_sapiens.%s.%s.fa.gz' % (fa, ASSEMBLY, DNA_FLAVOUR),
        'cds': '%s/cds/Homo_sapiens.%s.cds.all.fa.gz' % (fa, ASSEMBLY),
        'pep': '%s/pep/Homo_sapiens.%s.pep.all.fa.gz' % (fa, ASSEMBLY),
    }


def remote_size(url):
    """Content-Length, or None if the server will not say."""
    request = urllib.request.Request(url, method='HEAD')
    try:
        with urllib.request.urlopen(request) as r:
            length = r.headers.get('Content-Length')
            return int(length) if length else None
    except Exception:
        return None


def _fetch_range(args):
    """Streams one byte range to its own part file, so peak memory stays small
    however large the file is."""
    url, start, end, path = args
    request = urllib.request.Request(url, headers={'Range': 'bytes=%d-%d' % (start, end)})
    with urllib.request.urlopen(request) as r, open(path, 'wb') as out:
        shutil.copyfileobj(r, out, 1024 * 1024)
    got = os.path.getsize(path)
    if got != end - start + 1:
        raise IOError('range %d-%d returned %d bytes' % (start, end, got))
    return path


def download(url, dest, streams=DEFAULT_STREAMS):
    """Fetches a release file, in parallel byte ranges where possible.

    ftp.ensembl.org throttles hard per connection -- a single stream pulls the
    880 MB assembly at around 250 kB/s, an hour of waiting -- but it serves
    range requests and several connections together go several times faster.
    Falls back to one plain stream if the server will not cooperate.
    """
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print('  have %s' % os.path.basename(dest))
        return dest
    tmp = dest + '.part'
    size = remote_size(url) if streams > 1 else None
    started = time.time()
    if size and size > 32 * 1024 * 1024:
        print('  fetching %s (%.0f MB, %d streams)'
              % (os.path.basename(url), size / 1e6, streams))
        chunk = -(-size // streams)
        ranges = [(url, i * chunk, min(size, (i + 1) * chunk) - 1,
                   '%s.range%02d' % (tmp, i)) for i in range(streams)]
        try:
            with ThreadPoolExecutor(max_workers=streams) as pool:
                parts = list(pool.map(_fetch_range, ranges))
            with open(tmp, 'wb') as out:
                for part in parts:
                    with open(part, 'rb') as fh:
                        shutil.copyfileobj(fh, out, 1024 * 1024)
            if os.path.getsize(tmp) != size:
                raise IOError('assembled %d bytes, expected %d'
                              % (os.path.getsize(tmp), size))
        except Exception as e:
            print('  parallel fetch failed (%s); falling back to one stream' % e)
            size = None
        finally:
            for _, _, _, part in ranges:
                if os.path.exists(part):
                    os.remove(part)
    if not size or not os.path.exists(tmp):
        print('  fetching %s' % url)
        with urllib.request.urlopen(url) as r, open(tmp, 'wb') as out:
            shutil.copyfileobj(r, out, 1024 * 1024)
    os.rename(tmp, dest)
    print('  %.0f MB in %.0fs' % (os.path.getsize(dest) / 1e6, time.time() - started))
    return dest


# ---------------------------------------------------------------- GTF parsing

def parse_attributes(field):
    """GTF attributes: key "value"; key "value"; ..."""
    attrs = {}
    for part in field.rstrip().rstrip(';').split(';'):
        part = part.strip()
        if not part:
            continue
        key, _, value = part.partition(' ')
        attrs[key] = value.strip().strip('"')
    return attrs


WANTED_FEATURES = ('transcript', 'exon', 'CDS', 'stop_codon')


def read_gtf(path):
    """Streams the GTF into per-transcript records.

    Coordinates stay exactly as the GTF gives them: 1-based, inclusive, and
    always start <= end regardless of strand — the same convention the Ensembl
    REST API uses, so downstream code needs no adjustment.
    """
    transcripts = {}
    assembly = ASSEMBLY
    seen = 0
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                if line.startswith('#!genome-version'):
                    # "#!genome-version GRCh38" -- what the REST API reports as
                    # assembly_name, and what lands in the output's assembly
                    # column. Not #!genome-build, which carries a patch suffix,
                    # nor #!genome-build-accession.
                    assembly = line.split()[1]
                continue
            cols = line.split('\t')
            feature = cols[2]
            if feature not in WANTED_FEATURES:
                continue
            attrs = parse_attributes(cols[8])
            tr = attrs.get('transcript_id')
            if not tr:
                continue
            start, end = int(cols[3]), int(cols[4])
            rec = transcripts.get(tr)
            if rec is None:
                rec = transcripts[tr] = {
                    'exons': [], 'cds': [], 'stop': [],
                    'seq_region': cols[0],
                    'strand': 1 if cols[6] == '+' else -1,
                }
            if feature == 'transcript':
                rec['start'] = start
                rec['end'] = end
                rec['gene_id'] = attrs.get('gene_id', '')
                rec['gene_name'] = attrs.get('gene_name', '')
                # Ensembl's REST display_name for a transcript is its
                # transcript_name (e.g. PPP1R21-201). Not every transcript has
                # one; fall back to the ID so downstream string handling works.
                rec['display_name'] = attrs.get('transcript_name') or tr
                rec['biotype'] = attrs.get('transcript_biotype', '')
            elif feature == 'exon':
                rec['exons'].append((int(attrs.get('exon_number', 0)), start, end))
            elif feature == 'CDS':
                rec['cds'].append((int(attrs.get('exon_number', 0)), start, end))
            elif feature == 'stop_codon':
                rec['stop'].append((int(attrs.get('exon_number', 0)), start, end))
            seen += 1
            if seen % 2000000 == 0:
                print('    %d GTF features' % seen)
    return transcripts, assembly


def order_blocks(blocks, strand):
    """Transcript order: 5' to 3'.

    The REST API returns exons and CDS mappings in transcript order, and
    get_utrs relies on index 0 being the first exon of the transcript, so a
    minus-strand transcript comes back in descending genomic order.
    """
    blocks = sorted(blocks, key=lambda b: b[1])
    if strand == -1:
        blocks.reverse()
    return [(s, e) for _, s, e in blocks]


def merge_stop_codon(cds, stop, strand):
    """Folds the stop codon into the CDS blocks.

    The GTF stores stop_codon as a feature separate from CDS, but Ensembl's CDS
    coordinate system — what /map/cds returns and what cds.all.fa contains —
    includes it. Leaving it out would shift every downstream coordinate by 3 bp.

    Both lists are already in transcript order; the stop codon abuts the last
    CDS block (or, when a stop codon is split across an intron, the last two).
    """
    if not cds or not stop:
        return cds
    blocks = list(cds)
    for s, e in stop:
        last_s, last_e = blocks[-1]
        if strand == 1 and s == last_e + 1:
            blocks[-1] = (last_s, e)
        elif strand == -1 and e == last_s - 1:
            blocks[-1] = (s, last_e)
        else:
            blocks.append((s, e))
    return blocks


# -------------------------------------------------------------- FASTA parsing

def read_fasta_by_transcript(path, key_attr=None):
    """Maps transcript ID (version stripped) -> sequence.

    cds.all.fa is keyed by transcript ID; pep.all.fa is keyed by protein ID and
    names its transcript in a `transcript:ENST...` header attribute.
    """
    seqs = {}
    with gzip.open(path, 'rt') as fh:
        for title, seq in SimpleFastaParser(fh):
            fields = title.split()
            if key_attr:
                keys = [f.split(':', 1)[1] for f in fields[1:]
                        if f.startswith(key_attr + ':')]
                if not keys:
                    continue
                key = keys[0]
            else:
                key = fields[0]
            seqs[key.split('.')[0]] = seq
    return seqs


# ----------------------------------------------------------------- the bundle

SCHEMA = """
CREATE TABLE transcript (
    transcript_id    TEXT PRIMARY KEY,
    display_name     TEXT NOT NULL,
    gene_id          TEXT NOT NULL,
    gene_name        TEXT,
    biotype          TEXT,
    seq_region       TEXT NOT NULL,
    strand           INTEGER NOT NULL,
    start            INTEGER NOT NULL,
    end              INTEGER NOT NULL,
    exons            TEXT NOT NULL,   -- JSON [[start,end],...] in transcript order
    cds              TEXT,            -- JSON, transcript order, stop codon included
    cds_sequence     TEXT,
    protein_sequence TEXT
);
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT);
"""


def build_db(path, transcripts, assembly, release, cds_seqs, pep_seqs):
    tmp = path + '.part'
    if os.path.exists(tmp):
        os.remove(tmp)
    db = sqlite3.connect(tmp)
    db.executescript(SCHEMA)
    rows = []
    skipped = 0
    for tr, rec in transcripts.items():
        if 'start' not in rec:
            # exon/CDS lines without a transcript line: not usable
            skipped += 1
            continue
        strand = rec['strand']
        exons = order_blocks(rec['exons'], strand)
        cds = merge_stop_codon(order_blocks(rec['cds'], strand),
                               order_blocks(rec['stop'], strand), strand)
        rows.append((
            tr, rec['display_name'], rec['gene_id'], rec['gene_name'], rec['biotype'],
            rec['seq_region'], strand, rec['start'], rec['end'],
            json.dumps(exons, separators=(',', ':')),
            json.dumps(cds, separators=(',', ':')) if cds else None,
            cds_seqs.get(tr), pep_seqs.get(tr),
        ))
    db.executemany('INSERT INTO transcript VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)', rows)
    db.executemany('INSERT INTO meta VALUES (?,?)', [
        ('release', str(release)),
        ('assembly', assembly),
        ('built', datetime.now().isoformat(timespec='seconds')),
        ('transcripts', str(len(rows))),
    ])
    db.commit()
    db.execute('VACUUM')
    db.close()
    os.rename(tmp, path)
    print('  %d transcripts (%d skipped)' % (len(rows), skipped))


def _bgzip_with_htslib(src_gz, tmp):
    """Recompresses via the bgzip binary, which is far faster than Biopython's
    pure-Python writer. Returns False if bgzip is not installed."""
    bgzip = shutil.which('bgzip')
    if not bgzip:
        return False
    print('    using %s' % bgzip)
    with open(src_gz, 'rb') as fin, open(tmp, 'wb') as fout:
        unzip = subprocess.Popen(['gzip', '-dc'], stdin=fin, stdout=subprocess.PIPE)
        compress = subprocess.Popen([bgzip, '-c', '-@', '4'],
                                    stdin=unzip.stdout, stdout=fout)
        unzip.stdout.close()
        compress.communicate()
        unzip.wait()
    if compress.returncode or unzip.returncode:
        raise IOError('bgzip pipeline failed (%s/%s)' % (unzip.returncode,
                                                         compress.returncode))
    return True


def build_genome(src_gz, dest):
    """Recompresses the assembly from gzip to bgzip so it can be seeked into.

    Bytes are copied through unchanged; case normalisation happens on read in
    LocalSource, which also covers a bundle built from a masked assembly by
    mistake.
    """
    tmp = dest + '.part'
    started = time.time()
    if not _bgzip_with_htslib(src_gz, tmp):
        # No htslib on PATH: Biopython can write bgzf too, just slowly.
        with gzip.open(src_gz, 'rb') as fin, bgzf.BgzfWriter(tmp, 'wb') as fout:
            copied = 0
            while True:
                block = fin.read(4 * 1024 * 1024)
                if not block:
                    break
                fout.write(block)
                copied += len(block)
                if copied % (512 * 1024 * 1024) < 4 * 1024 * 1024:
                    print('    %.1f GB in %.0fs' % (copied / 1e9, time.time() - started))
    print('    compressed in %.0fs' % (time.time() - started))
    os.rename(tmp, dest)
    print('  indexing %s' % os.path.basename(dest))
    from pyfaidx import Fasta
    Fasta(dest)  # writes .fai and .gzi alongside
    return dest


def verify(db_path, genome_path):
    """Checks the annotation and the assembly actually describe each other.

    A GTF paired with the wrong assembly, or with a FASTA that uses UCSC-style
    "chr1" names, would otherwise produce empty sequences transcript by
    transcript rather than failing outright.
    """
    from pyfaidx import Fasta
    genome = Fasta(genome_path, as_raw=True)
    contigs = set(genome.keys())
    db = sqlite3.connect(db_path)
    regions = [r[0] for r in db.execute('SELECT DISTINCT seq_region FROM transcript')]
    missing = sorted(set(regions) - contigs)
    if missing:
        raise SystemExit(
            'Bundle is inconsistent: %d of %d sequence regions in the GTF are not '
            'in the assembly (e.g. %s).\nThe GTF and FASTA are probably from '
            'different releases or assemblies.'
            % (len(missing), len(regions), ', '.join(missing[:5])))
    # Every transcript must fit inside its contig, or the flanked slice is wrong.
    overruns = 0
    for region, longest in db.execute(
            'SELECT seq_region, MAX(end) FROM transcript GROUP BY seq_region'):
        if longest > len(genome[region]):
            overruns += 1
    if overruns:
        raise SystemExit('Bundle is inconsistent: %d contigs are shorter than the '
                         'transcripts annotated on them.' % overruns)
    print('  checked %d sequence regions against the assembly' % len(regions))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--release', default='116',
                        help='Ensembl release to build from (default: 116)')
    parser.add_argument('--out', default='refdata',
                        help='Directory to write the bundle into (default: refdata)')
    parser.add_argument('--download-dir', default=None,
                        help='Where to keep the downloaded release files '
                             '(default: <out>/_download)')
    parser.add_argument('--force', action='store_true',
                        help='Rebuild outputs that already exist')
    parser.add_argument('--streams', type=int, default=DEFAULT_STREAMS,
                        help='Parallel connections per download (default: %d). '
                             'Use 1 if a proxy mishandles range requests.'
                             % DEFAULT_STREAMS)
    parser.add_argument('--stage', choices=['all', 'db', 'genome'], default='all',
                        help='Build only part of the bundle (default: all). The '
                             'database needs the GTF and the CDS/protein FASTA; '
                             'the genome needs only the assembly FASTA.')
    args = parser.parse_args(argv)

    bundle = os.path.join(args.out, 'ensembl-%s' % args.release)
    downloads = args.download_dir or os.path.join(args.out, '_download')
    os.makedirs(bundle, exist_ok=True)
    os.makedirs(downloads, exist_ok=True)

    urls = source_urls(args.release)
    needed = {'all': list(urls), 'db': ['gtf', 'cds', 'pep'], 'genome': ['dna']}[args.stage]
    print('Fetching release %s' % args.release)
    paths = {k: download(urls[k], os.path.join(downloads, os.path.basename(urls[k])),
                         streams=args.streams)
             for k in needed}

    db_path = os.path.join(bundle, TRANSCRIPTS_DB)
    genome_path = os.path.join(bundle, GENOME_FA)

    if args.stage == 'genome':
        print('Skipping %s (--stage genome)' % db_path)
    elif args.force or not os.path.exists(db_path):
        print('Reading GTF')
        transcripts, assembly = read_gtf(paths['gtf'])
        print('  %d transcripts, assembly %s' % (len(transcripts), assembly))
        print('Reading CDS and protein FASTA')
        cds_seqs = read_fasta_by_transcript(paths['cds'])
        pep_seqs = read_fasta_by_transcript(paths['pep'], key_attr='transcript')
        print('  %d CDS, %d protein sequences' % (len(cds_seqs), len(pep_seqs)))
        print('Writing %s' % db_path)
        build_db(db_path, transcripts, assembly, args.release, cds_seqs, pep_seqs)
    else:
        print('Have %s' % db_path)

    if args.stage == 'db':
        print('Skipping %s (--stage db)' % genome_path)
    elif args.force or not os.path.exists(genome_path + '.fai'):
        print('Writing %s' % genome_path)
        build_genome(paths['dna'], genome_path)
    else:
        print('Have %s' % genome_path)

    n_tr = None
    if os.path.exists(db_path):
        n_tr = int(sqlite3.connect(db_path).execute(
            "SELECT value FROM meta WHERE key = 'transcripts'").fetchone()[0])

    if args.stage == 'all' or (os.path.exists(db_path)
                               and os.path.exists(genome_path + '.fai')):
        verify(db_path, genome_path)

    manifest = os.path.join(bundle, 'manifest.json')
    with open(manifest, 'w') as fh:
        json.dump({
            'release': args.release,
            'assembly': ASSEMBLY,
            'species': SPECIES,
            'dna_flavour': DNA_FLAVOUR,
            'built': date.today().isoformat(),
            'sources': urls,
            'transcripts': n_tr,
        }, fh, indent=2, sort_keys=True)
        fh.write('\n')
    print('Wrote %s' % manifest)
    print('\nDone. Use it with:\n  --source local --refdata %s' % args.out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
