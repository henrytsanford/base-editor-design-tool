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
    return 'results/%s/%s' % (key[:2], key)


def result_key(key, name):
    """Maps a result name through a fixed table. The name never reaches a path."""
    return '%s/%s' % (result_prefix(key), RESULT_FILES[name])


def manifest_key(key):
    return '%s/%s' % (result_prefix(key), MANIFEST)
