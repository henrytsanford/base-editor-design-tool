"""What a design looks like before anyone asks for a row.

A finished run is thousands of guides, and handing those to a model is neither
readable nor affordable. What a researcher actually decides from is the shape of the
result: how many guides knock the protein out, which ClinVar classifications are
reachable, which exons are covered, which residue changes are on offer. All of that
is already computed -- ResultTable for the filter vocabularies, Coverage for the
panels the web page draws -- so this reads it off rather than counting anything twice.

Rows come from query_guides, filtered by exactly the values reported here.
"""
from service.results import EMPTY_FACETS, ROWS_PER_PAGE, page_of

# The matrix is 20x21 squares, most of them empty. Enough of the busiest to show what
# a screen could target, without turning a summary into a table.
TOP_SUBSTITUTIONS = 15


def exon_coverage(table):
    """Guides per exon, with each exon's length and whether it codes.

    Numbered from 1, the way the exon filter takes it. Empty when the bundle could not
    answer for the transcript, which is also when the web page draws no map.
    """
    if table.coverage is None:
        return []
    geometry = table.coverage.geometry
    counts = table.coverage.exon_guides
    return [{'exon': i + 1,
             'bp': int(high - low + 1),
             'coding': any(part[2] for part in geometry.parts[i]),
             'guides': int(counts[i])}
            for i, (low, high) in enumerate(geometry.exons)]


def top_substitutions(table, limit=TOP_SUBSTITUTIONS):
    """The residue changes this run can make, busiest first.

    `sites` is how many distinct residues the change is available at and `pathogenic`
    how many of its guides recreate a ClinVar variant classified pathogenic -- which
    is usually the number that decides whether a substitution is worth a screen.
    """
    if table.coverage is None:
        return []
    counts = table.coverage.matrix_counts
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:limit]
    return [{'substitution': key,
             'guides': int(count),
             'sites': int(table.coverage.matrix_sites.get(key, 0)),
             'pathogenic': int(table.coverage.matrix_pathogenic.get(key, 0))}
            for key, count in ranked]


def build(manifest, table):
    """The design_guides payload: everything but the guides.

    `counts.errors` is reported even when it is zero because a non-zero one is load
    bearing -- the engine abandons a transcript mid-way on a CDS error, so a run with
    errors is a partial answer that otherwise looks complete.
    """
    facets = table.facets() if table is not None else EMPTY_FACETS
    designs = manifest['designs']
    _, pages, _ = page_of(designs, 1)
    return {
        'run_id': manifest['run_id'],
        'target': {'kind': manifest['kind'],
                   'transcript_id': manifest['transcript_id'],
                   'gene': manifest['gene'],
                   'name': manifest['label']},
        'params': manifest['params'],
        'engine_version': manifest['engine_version'],
        'reference': manifest['reference'],
        'clinvar': manifest['clinvar'],
        'counts': {'designs': designs,
                   'errors': manifest['errors'],
                   'annotations': manifest['annotations']},
        'consequence': facets['consequence'],
        'significance': [{'value': value, 'count': count}
                         for value, count in facets['significance']],
        'any_match': facets['any_match'],
        'edits': facets['edits'],
        'exon_coverage': exon_coverage(table) if table is not None else [],
        'substitutions_top': top_substitutions(table) if table is not None else [],
        'runtime_seconds': manifest['runtime_seconds'],
        'rows_per_page': ROWS_PER_PAGE,
        'pages': pages,
        'files': manifest['files'],
    }
