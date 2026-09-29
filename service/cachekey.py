"""The cache key.

sha256 of canonical JSON over the resolved parameters, so a preset and the explicit
parameters it resolves to share one result, and a new reference release or a bumped
engine version gives new keys rather than stale answers.

The digest is the only thing that ever becomes a path segment, which is what keeps
user input out of the filesystem.
"""
import dataclasses
import hashlib
import json

RESULT_FILES = {
    'designs': 'designs.tsv.gz',
    'errors': 'errors.tsv.gz',
    'clinvar': 'clinvar.tsv.gz',
}
MANIFEST = 'manifest.json'
# The prefix every result lives under, and what storage.evict sweeps.
RESULTS = 'results'


@dataclasses.dataclass(frozen=True)
class Target:
    """What a run designs over: a transcript, or a sequence given inline.

    Cache-key vocabulary, which is why it lives here: `key` is the part of the digest
    that says what was designed, and `files` names keys of RESULT_FILES. Both front
    doors and the pool worker need it, and this is the module below all of them --
    a Target travels to a worker process by pickle, so it must import cleanly there.
    """
    kind: str
    transcript_id: str = ''
    name: str = ''
    sequence: str = ''

    @property
    def label(self):
        """What the output files are named after, as the CLI names them."""
        return self.transcript_id if self.kind == 'transcript' else self.name

    @property
    def key(self):
        """The part of the cache key that says what was designed.

        A transcript is its own ID, so every key computed before sequences existed is
        unchanged. A sequence is keyed by its digest rather than by itself: the key is
        canonical JSON that gets hashed anyway, and a 100 kb sequence has no business
        being built into a string first.
        """
        if self.kind == 'transcript':
            return self.transcript_id
        digest = hashlib.sha256(
            ('%s\n%s' % (self.name, self.sequence)).encode('utf-8')).hexdigest()
        return 'seq:%s' % digest

    @property
    def files(self):
        """Which of the three files this kind of run writes.

        Nucleotide input has no gene to look up, so it never annotates -- the CLI does
        not create the file at all, and neither does this.
        """
        return ('designs', 'errors', 'clinvar') if self.kind == 'transcript' \
            else ('designs', 'errors')


def key_fields(transcript_id, params, ensembl_release, clinvar_version, engine_version):
    """Ordering is handled by sort_keys.

    The design parameters come from the dataclass rather than a second list of
    names, so a parameter added to `DesignParams` changes the key instead of being
    silently left out of it -- which would serve a cached result computed under a
    different value.
    """
    return dict(
        dataclasses.asdict(params),
        transcript_id=transcript_id,
        ensembl_release=str(ensembl_release),
        clinvar_version=str(clinvar_version),
        engine_version=str(engine_version),
    )


def cache_key(transcript_id, params, ensembl_release, clinvar_version, engine_version):
    fields = key_fields(transcript_id, params, ensembl_release, clinvar_version,
                        engine_version)
    canonical = json.dumps(fields, sort_keys=True, separators=(',', ':'),
                           ensure_ascii=True)
    return hashlib.sha256(canonical.encode('utf-8')).hexdigest()


def result_prefix(key):
    """Where a key's objects live. Two hex characters of fan-out keeps any one
    directory from collecting every result the service has ever computed."""
    return '%s/%s/%s' % (RESULTS, key[:2], key)


def result_key(key, name):
    """Maps a result name through a fixed table. The name never reaches a path."""
    return '%s/%s' % (result_prefix(key), RESULT_FILES[name])


def manifest_key(key):
    return '%s/%s' % (result_prefix(key), MANIFEST)
