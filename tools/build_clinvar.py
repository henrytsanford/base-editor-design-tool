"""Builds the ClinVar variant database from NCBI's variant_summary.txt.

The design script needs one thing from ClinVar: the GRCh38 SNVs for a gene. Loading
the whole 4.1-million-row table to answer that costs ~40 s and ~2 GB on every run, so
this script distills it once into an indexed SQLite database:

    <out>/clinvar-<date>.db    one `variant` table, indexed on GeneSymbol

Run once per ClinVar release (NCBI publishes weekly):

    python tools/build_clinvar.py                          # downloads the latest
    python tools/build_clinvar.py --variant-summary variant_summary.txt

Alleles come from the VCF-normalized columns: `ReferenceAllele` and `AlternateAllele`
are `na` on effectively every row, so the alleles worth having are in
`ReferenceAlleleVCF` / `AlternateAlleleVCF`. They are renamed to the internal names
here so nothing downstream has to know.
"""
import argparse
import gzip
import os
import shutil
import sqlite3
import sys
import urllib.request
from datetime import date, datetime

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from transcript_source import VARIANT_COLUMNS, VARIANT_TABLE

SOURCE_URL = ('https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/'
              'variant_summary.txt.gz')

# The columns read out of variant_summary.txt: what the output keeps, plus the four
# used only to select GRCh38 SNVs on a real chromosome.
USECOLS = ['#AlleleID', 'GeneSymbol', 'Name', 'ClinicalSignificance', 'PhenotypeList',
           'Assembly', 'Chromosome', 'Type', 'PositionVCF', 'ReferenceAlleleVCF',
           'AlternateAlleleVCF', 'ReviewStatus']

SCHEMA = '''
CREATE TABLE variant (
    "#AlleleID"           TEXT,
    RefSeqID              TEXT,
    Name                  TEXT,
    GeneSymbol            TEXT,
    ClinicalSignificance  TEXT,
    PhenotypeList         TEXT,
    ClinVar_SNP_Position  INTEGER,
    ReferenceAllele       TEXT,
    AlternateAllele       TEXT,
    ReviewStatus          TEXT
);
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT);
'''


def parse_variant_df(variant_df):
    """Cleans the ClinVar table down to the GRCh38 SNVs the design script can use.

    Keeps the internal column names `ClinVar_SNP_Position` / `ReferenceAllele` /
    `AlternateAllele`, which is what `get_snps` and both output files are written
    against.
    """
    # Remove non-GRCh38 rows, non-SNPs, and chromosomes other than 1-22, X, and Y.
    # A missing VCF allele (8 rows out of 4.1 M) leaves nothing to compare against.
    temp_variant_df = variant_df[(variant_df.Assembly == 'GRCh38')
                                & (variant_df.Type == 'single nucleotide variant')
                                & (variant_df.Chromosome != 'MT')
                                & (variant_df.Chromosome != 'na')
                                & (variant_df.ReferenceAlleleVCF != 'na')
                                & (variant_df.AlternateAlleleVCF != 'na')]
    parsed_variant_df = temp_variant_df.copy()
    parsed_variant_df = parsed_variant_df.reset_index(drop=True)
    # PositionVCF agrees with Start on all but 8 SNV rows, but it is the coordinate
    # that belongs with the VCF alleles.
    parsed_variant_df = parsed_variant_df.rename(columns={
        'PositionVCF': 'ClinVar_SNP_Position',
        'ReferenceAlleleVCF': 'ReferenceAllele',
        'AlternateAlleleVCF': 'AlternateAllele'})
    parsed_variant_df = parsed_variant_df.assign(
        RefSeqID=parsed_variant_df['Name'].str.split(pat='(', n=1).str[0]
    )
    return parsed_variant_df[list(VARIANT_COLUMNS)]


def read_variant_summary(path):
    """The subset of variant_summary.txt this builder needs, as a DataFrame."""
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        return pd.read_table(fh, usecols=USECOLS, dtype=str)


def build(variant_summary_path, db_path):
    """Writes the indexed variant database. Returns the number of rows kept."""
    parsed = parse_variant_df(read_variant_summary(variant_summary_path))
    # get_snps matches positions with .isin() against a list of ints, so the column
    # has to arrive as an integer and not as text.
    parsed['ClinVar_SNP_Position'] = parsed['ClinVar_SNP_Position'].astype('int64')

    tmp_path = db_path + '.tmp'
    if os.path.exists(tmp_path):
        os.remove(tmp_path)
    db = sqlite3.connect(tmp_path)
    try:
        db.executescript(SCHEMA)
        parsed.to_sql(VARIANT_TABLE, db, if_exists='append', index=False,
                      chunksize=100000)
        db.execute('CREATE INDEX idx_variant_gene ON variant(GeneSymbol)')
        db.executemany('INSERT INTO meta (key, value) VALUES (?, ?)', [
            ('built', datetime.now().isoformat(timespec='seconds')),
            ('source', os.path.abspath(variant_summary_path)),
            ('rows', str(len(parsed))),
        ])
        db.commit()
    finally:
        db.close()
    os.replace(tmp_path, db_path)
    return len(parsed)


def download(url, dest):
    """Fetches variant_summary.txt.gz, decompressing it on the way in."""
    print('Fetching %s' % url)
    with urllib.request.urlopen(url) as r, gzip.GzipFile(fileobj=r) as gz, \
            open(dest, 'wb') as out:
        shutil.copyfileobj(gz, out)
    return dest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--variant-summary', default=None,
                        help='Path to variant_summary.txt (or .gz). Downloaded from '
                             'NCBI if not given.')
    parser.add_argument('--out', default='refdata',
                        help='Directory to write the database into (default: refdata)')
    parser.add_argument('--download-dir', default=None,
                        help='Where to keep the downloaded variant_summary.txt '
                             '(default: <out>/_download)')
    parser.add_argument('--date', default=None,
                        help='Version stamp for the filename, clinvar-<date>.db '
                             "(default: today). variant_summary.txt does not name "
                             'its own release, so this is the date it was fetched.')
    parser.add_argument('--force', action='store_true',
                        help='Rebuild the database even if it already exists')
    args = parser.parse_args(argv)

    os.makedirs(args.out, exist_ok=True)
    source = args.variant_summary
    if source is None:
        downloads = args.download_dir or os.path.join(args.out, '_download')
        os.makedirs(downloads, exist_ok=True)
        source = download(SOURCE_URL, os.path.join(downloads, 'variant_summary.txt'))

    db_path = os.path.join(args.out, 'clinvar-%s.db'
                           % (args.date or date.today().isoformat()))
    if os.path.exists(db_path) and not args.force:
        print('Have %s (use --force to rebuild)' % db_path)
        return 0

    print('Reading %s' % source)
    print('Writing %s' % db_path)
    rows = build(source, db_path)
    print('  %d variants' % rows)
    print('\nDone. Use it with:\n  --clinvar-db %s' % db_path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
