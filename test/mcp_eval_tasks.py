"""Design tasks, and the route through the tools each one should take.

Each task is a question a researcher would actually ask, paired with the sequence of
calls a competent agent should make to answer it and the things that must be true of
what comes back. No model runs: the plan is written out, so the eval measures the
tool surface rather than a model's mood, and it can run as an ordinary test.

The prompts are not decoration. They are what a live evaluation would send, and they
are what makes it checkable that the route below actually answers the question asked.

Every task names a golden case, and every task that returns rows has them checked
against that case's frozen output -- filtering is a way of reading a result, never a
way of inventing rows.
"""
from collections import namedtuple

from test_mcp_golden import GOLDEN_PARAMS, gfp_sequence

Task = namedtuple('Task', 'name prompt golden plan expect')

# A step is (tool, build), where build turns the results so far into arguments.
# `r` is keyed by tool name, most recent call winning.


def _mane(r):
    return r['list_transcripts']['transcripts'][0]['transcript_id']


def _run(r):
    return r['design_guides']['run_id']


TASKS = [
    Task(
        name='map2k1_pathogenic_missense',
        prompt=('I want to model the MAP2K1 missense variants that ClinVar already '
                'calls pathogenic, using a cytosine base editor. Which guides '
                'should I order?'),
        golden='map2k1',
        plan=[
            ('list_editors', lambda r: {}),
            ('resolve_gene', lambda r: {'symbol': 'MAP2K1'}),
            ('list_transcripts', lambda r: {'gene': r['resolve_gene']['gene']}),
            ('design_guides', lambda r: dict(GOLDEN_PARAMS, transcript_id=_mane(r))),
            ('query_guides', lambda r: {'run_id': _run(r), 'consequence': 'mis',
                                        'significance': 'any-match',
                                        'deaminase': 'C-T'}),
        ],
        expect=[
            ('the editor list offers a cytosine editor',
             lambda r: any(e['edit'] == 'C-T' for e in r['list_editors']['editors'])),
            ('the gene resolved exactly',
             lambda r: r['resolve_gene']['gene'] == 'MAP2K1'),
            ('the recommended transcript is MANE Select',
             lambda r: r['list_transcripts']['transcripts'][0]['mane_select'] == 1),
            ('the run designed guides',
             lambda r: r['design_guides']['counts']['designs'] > 0),
            ('the run completed rather than stopping early',
             lambda r: r['design_guides']['counts']['errors'] == 0),
            ('a pathogenic classification was on offer to filter by',
             lambda r: any(s['value'] == 'Pathogenic'
                           for s in r['design_guides']['significance'])),
            ('the filter narrowed the result',
             lambda r: 0 < r['query_guides']['matched']
             < r['design_guides']['counts']['designs']),
            ('every guide returned is from the cytosine pass',
             lambda r: all(row[15] == 'C-T' for row in r['query_guides']['rows'])),
            ('every guide returned makes a missense edit',
             lambda r: all('Missense' in row[20] for row in r['query_guides']['rows'])),
        ],
    ),
    Task(
        name='isy1_minus_strand_exon',
        prompt=('ISY1 is on the minus strand. Show me the guides that sit in exon 3 '
                'and say how well the transcript is covered overall.'),
        golden='isy1',
        plan=[
            ('resolve_gene', lambda r: {'symbol': 'ISY1'}),
            ('list_transcripts', lambda r: {'gene': 'ISY1'}),
            ('design_guides', lambda r: dict(GOLDEN_PARAMS,
                                             transcript_id='ENST00000393295')),
            ('query_guides', lambda r: {'run_id': _run(r), 'exon': 3}),
        ],
        expect=[
            ('the gene is on the minus strand',
             lambda r: any(t['strand'] == -1
                           for t in r['list_transcripts']['transcripts'])),
            ('coverage is reported per exon',
             lambda r: len(r['design_guides']['exon_coverage']) > 1),
            ('every exon reports its length and whether it codes',
             lambda r: all('bp' in e and 'coding' in e
                           for e in r['design_guides']['exon_coverage'])),
            ('exon 3 holds the guides the map counted',
             lambda r: r['query_guides']['matched']
             == r['design_guides']['exon_coverage'][2]['guides']),
        ],
    ),
    Task(
        name='krtap25_1_single_exon',
        prompt=('KRTAP25-1 is a single-exon gene. Give me every guide that knocks '
                'the protein out.'),
        golden='krtap25_1',
        plan=[
            ('resolve_gene', lambda r: {'symbol': 'KRTAP25-1'}),
            ('design_guides', lambda r: dict(GOLDEN_PARAMS,
                                             transcript_id='ENST00000416044')),
            ('query_guides', lambda r: {'run_id': _run(r), 'consequence': 'lof'}),
        ],
        expect=[
            ('a hyphenated symbol resolves',
             lambda r: r['resolve_gene']['gene'] == 'KRTAP25-1'),
            ('the transcript has one exon',
             lambda r: len(r['design_guides']['exon_coverage']) == 1),
            ('the knockout filter is answerable',
             lambda r: r['query_guides']['matched'] >= 0),
        ],
    ),
    Task(
        name='hgs_partial_run_is_reported',
        prompt=("Design guides across HGS and tell me whether I can trust the "
                "result."),
        golden='hgs',
        plan=[
            ('design_guides', lambda r: dict(GOLDEN_PARAMS,
                                             transcript_id='ENST00000678176')),
            ('export_run', lambda r: {'run_id': _run(r), 'file': 'errors'}),
        ],
        expect=[
            # HGS has a CDS length that is not a multiple of 3, and the engine
            # abandons the transcript there. A summary that looked clean would be
            # the worst possible answer, so the error count has to carry it.
            ('the partial run is declared, not hidden',
             lambda r: r['design_guides']['counts']['errors'] > 0),
            ('the error file says what went wrong',
             lambda r: 'codon' in r['export_run']['data'].lower()),
        ],
    ),
    Task(
        name='malat1_non_coding',
        prompt='Can I base-edit MALAT1? Design over it and tell me what happened.',
        golden='malat1',
        plan=[
            ('design_guides', lambda r: dict(GOLDEN_PARAMS,
                                             transcript_id='ENST00000620902')),
            ('export_run', lambda r: {'run_id': _run(r), 'file': 'errors'}),
        ],
        expect=[
            # A non-coding transcript is a legitimate answer, not a crash: the tool
            # has to come back and say there is nothing to design.
            ('no guides were designed',
             lambda r: r['design_guides']['counts']['designs'] == 0),
            ('the reason is reported rather than raised',
             lambda r: r['design_guides']['counts']['errors'] == 1),
            ('the error names the transcript',
             lambda r: 'ENST00000620902' in r['export_run']['data']),
        ],
    ),
    Task(
        name='gfp_pasted_sequence',
        prompt=("Here is a GFP sequence from a plasmid that isn't in Ensembl. "
                'Design base-editor guides over it.'),
        golden='gfp',
        plan=[
            ('design_guides', lambda r: dict(GOLDEN_PARAMS,
                                             sequence=gfp_sequence()[1],
                                             sequence_name=gfp_sequence()[0])),
            ('query_guides', lambda r: {'run_id': _run(r), 'consequence': 'mis'}),
        ],
        expect=[
            ('guides were designed without a reference transcript',
             lambda r: r['design_guides']['counts']['designs'] > 0),
            ('no gene was invented for it',
             lambda r: r['design_guides']['target']['gene'] == ''),
            ('no ClinVar annotation is claimed',
             lambda r: r['design_guides']['clinvar'] == 'none'),
            ('it writes no annotations file',
             lambda r: r['design_guides']['files'] == ['designs', 'errors']),
        ],
    ),
]
