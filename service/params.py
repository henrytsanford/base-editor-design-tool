"""Turning a query string into a validated design request.

Everything a request can say is checked here, before a job is scheduled, so a
malformed or hostile request costs a regex rather than a CPU core. The rules are the
design doc's 5, plus the preset allowlist and the caps that bound the work of
checking.

A request carries two independent halves. The *design* parameters say what to compute
and are what the cache key is built from; the *view* parameters say how to show the
result and must never reach the key, or filtering a table would recompute it. The two
allowlists are disjoint, and asserted to be.

`ValidationError` carries a message written for the person who typed the URL. It
never quotes a filesystem path or an internal exception.
"""
import dataclasses
import re
from dataclasses import dataclass

from bedesign import DESIGN_COLUMNS
from bedesign.engine import ALL_EDITS, DesignParams, UnknownBaseEditor

from .cachekey import RESULT_FILES

# No value a request can carry is longer than this; anything longer is a probe.
MAX_VALUE_LEN = 64

TRANSCRIPT = re.compile(r'^ENST\d{11}$')
# IUPAC nucleotide codes, which is what revcom in the engine knows how to complement.
PAM = re.compile(r'^[ACGTRYSWKMBDHVN]{2,8}$')
WINDOW = re.compile(r'^(\d{1,2})-(\d{1,2})$')

SG_LEN_RANGE = (17, 24)
INTRON_BUFFER_RANGE = (0, 100)
EDITS = tuple(ALL_EDITS) + ('all',)

# Parameters that describe the editor itself. Giving these alongside a preset is
# ambiguous -- which one wins? -- so it is refused rather than silently resolved.
EDITOR_PARAMS = ('pam', 'window', 'sg_len', 'edit')
# Every design parameter is accepted by name, taken from the dataclass so a new one
# does not have to be added here as well to stop being rejected as unknown.
EDITOR_ALL = ('preset',) + tuple(f.name for f in dataclasses.fields(DesignParams))
DESIGN_PARAMS = ('transcript',) + EDITOR_ALL
# The gene search takes a symbol instead of a transcript, and carries the editor
# choice through so the links it renders arrive at /designs fully specified.
GENE_PARAMS = ('q',) + EDITOR_ALL
# A download names a stored file alongside the design that produced it.
DOWNLOAD_PARAMS = DESIGN_PARAMS + ('file',)

# How to show a result. None of these reach the cache key.
VIEW_PARAMS = ('mutation', 'significance', 'deaminase', 'strand', 'hide_bsmbi',
               'hide_4t', 'sort', 'dir', 'page')
# Deliberately 'deaminase' rather than 'edit': 'edit' is already a design parameter
# over the same C-T/A-G vocabulary, and one name for both would make a view setting
# change the cache key.
assert not set(VIEW_PARAMS) & set(DESIGN_PARAMS), 'a view parameter shadows a design one'

SORT_COLUMNS = frozenset(DESIGN_COLUMNS)
DIRECTIONS = ('asc', 'desc')
STRANDS = ('sense', 'antisense')
# Deep enough for TTN's 25,662 guides at 50 a page, and a bound on the arithmetic.
MAX_PAGE = 100000
# ClinVar classifications are free text from an external source: letters, digits and
# the punctuation its compound labels use, e.g. 'Benign/Likely benign' and
# 'Conflicting classifications of pathogenicity'. Shape only -- the value is checked
# against the result's own vocabulary once the frame is loaded.
FILTER_VALUE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 ,/()'._-]{0,63}$")
# A gene symbol, capped so a probe cannot make the lookup do work.
GENE_QUERY = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]{0,31}$')


class ValidationError(ValueError):
    """A request that will not be run, with a message safe to show the user."""


def _one(query, name):
    """The single value for a name, or None when it was not really given.

    Two values for one parameter is how a request tries to make the validator and
    the code that reads it disagree, so it is an error rather than last-wins.

    An empty value counts as absent. A GET form submits every field it has, including
    the ones left blank, so 'pam=' means 'I did not set a PAM' and must not be read as
    'my PAM is the empty string'.
    """
    values = query.getlist(name)
    if not values:
        return None
    if len(values) > 1:
        raise ValidationError('%s was given more than once.' % name)
    value = values[0].strip()
    if len(value) > MAX_VALUE_LEN:
        raise ValidationError('%s is too long.' % name)
    return value or None


def _int(value, name, low, high):
    try:
        number = int(value)
    except ValueError:
        raise ValidationError('%s must be a whole number.' % name)
    if not low <= number <= high:
        raise ValidationError('%s must be between %d and %d.' % (name, low, high))
    return number


def _flag(value, name):
    if value.lower() in ('true', '1', 'yes'):
        return True
    if value.lower() in ('false', '0', 'no'):
        return False
    raise ValidationError('%s must be true or false.' % name)


def _bool(query, name):
    """A checkbox: absent means unset, which is false."""
    value = _one(query, name)
    return _flag(value, name) if value else False


def _check_keys(query, allowed):
    """Refuses any parameter name not on the allowlist.

    Rejected rather than ignored: silently dropping a parameter would let a user
    believe a setting applied when the result ignored it.

    The allowlist is its own size bound -- no route can legitimately carry more names
    than it accepts -- so a longer request is refused before the names are compared.
    """
    keys = set(query.keys())
    if len(keys) > len(allowed):
        raise ValidationError('Too many parameters.')
    unknown = sorted(keys - set(allowed))
    if unknown:
        raise ValidationError('Unknown parameter: %s.' % ', '.join(unknown))


def parse_editor_params(query):
    """The editor half of a request, as DesignParams.

    Shared by /designs and /genes: the search page has to validate a preset before it
    knows which transcript the user will pick, and one implementation means the two
    routes cannot disagree about what a valid editor is.
    """
    overrides = {}
    intron_buffer = _one(query, 'intron_buffer')
    if intron_buffer is not None:
        overrides['intron_buffer'] = _int(intron_buffer, 'intron_buffer',
                                          *INTRON_BUFFER_RANGE)
    filter_gc = _one(query, 'filter_gc')
    if filter_gc is not None:
        overrides['filter_gc'] = _flag(filter_gc, 'filter_gc')

    preset = _one(query, 'preset')
    if preset:
        # On the values, not the names: the search form posts every Advanced field,
        # so a blank one is not a conflicting setting.
        conflicting = sorted(k for k in EDITOR_PARAMS if _one(query, k) is not None)
        if conflicting:
            raise ValidationError(
                'preset already sets %s; give one or the other, not both.'
                % ', '.join(conflicting))
        try:
            # The preset table in the engine is the allowlist, so the dropdown and
            # this check cannot drift apart.
            return DesignParams.from_preset(preset, **overrides)
        except UnknownBaseEditor:
            raise ValidationError('Unknown base editor %r.' % preset)

    fields = dict(overrides)
    pam = _one(query, 'pam')
    if pam is not None:
        pam = pam.upper()
        if not PAM.match(pam):
            raise ValidationError('pam must be 2-8 IUPAC nucleotide codes, e.g. NGG.')
        fields['pam'] = pam

    sg_len = _one(query, 'sg_len')
    if sg_len is not None:
        fields['sg_len'] = _int(sg_len, 'sg_len', *SG_LEN_RANGE)

    edit = _one(query, 'edit')
    if edit is not None:
        if edit not in EDITS:
            raise ValidationError('edit must be one of %s.' % ', '.join(EDITS))
        fields['edit'] = edit

    window = _one(query, 'window')
    if window is not None:
        match = WINDOW.match(window)
        if not match:
            raise ValidationError('window must look like 4-8.')
        start, end = int(match.group(1)), int(match.group(2))
        # Checked against the sg_len this request will actually use, not the default.
        limit = fields.get('sg_len', DesignParams().sg_len)
        if not 1 <= start <= end <= limit:
            raise ValidationError(
                'window must satisfy 1 <= start <= end <= sg_len (%d).' % limit)
        fields['window'] = window

    return DesignParams(**fields)


def parse_designs_query(query, transcript_exists=None):
    """Validates /designs parameters into (transcript_id, DesignParams).

    View parameters are accepted here so a filtered table is a legal URL, but they are
    not returned: only the design half reaches the cache key.

    `transcript_exists` is called only after the ID matches the pattern, so an
    unparseable ID never reaches the database.
    """
    _check_keys(query, DESIGN_PARAMS + VIEW_PARAMS)
    return _transcript(query, transcript_exists), parse_editor_params(query)


def _transcript(query, transcript_exists):
    transcript = _one(query, 'transcript')
    if not transcript:
        raise ValidationError('A transcript is required, e.g. ENST00000307102.')
    transcript = transcript.upper()
    if not TRANSCRIPT.match(transcript):
        raise ValidationError(
            'Transcript must look like ENST00000307102: ENST followed by 11 digits.')
    if transcript_exists is not None and not transcript_exists(transcript):
        raise ValidationError('Transcript %s is not in the reference bundle.' % transcript)
    return transcript


def parse_genes_query(query):
    """Validates /genes parameters into (symbol, DesignParams).

    The symbol may be empty: an empty search box is a blank results page, not an error.
    """
    _check_keys(query, GENE_PARAMS)
    symbol = _one(query, 'q') or ''
    if symbol and not GENE_QUERY.match(symbol):
        raise ValidationError(
            'A gene symbol is letters, digits, dot, dash or underscore, e.g. MAP2K1.')
    return symbol.upper(), parse_editor_params(query)


@dataclass(frozen=True)
class TableView:
    """How to show a result: filters, sort and page. Never part of the cache key."""
    mutation: str = ''
    significance: str = ''
    deaminase: str = ''
    strand: str = ''
    hide_bsmbi: bool = False
    hide_4t: bool = False
    sort: str = ''
    dir: str = 'asc'
    page: int = 1

    @property
    def filtered(self):
        return bool(self.mutation or self.significance or self.deaminase
                    or self.strand or self.hide_bsmbi or self.hide_4t)


def _choice(query, name, allowed):
    value = _one(query, name)
    if not value:
        return ''
    if value not in allowed:
        raise ValidationError('%s must be one of %s.' % (name, ', '.join(allowed)))
    return value


def _open_value(query, name):
    """A filter whose vocabulary comes from an external source.

    Only the shape is checked here. ClinVar's classification list grows without asking
    us, so the value is matched against the result's own vocabulary once the frame is
    loaded -- see results.UnknownFilterValue.
    """
    value = _one(query, name)
    if not value:
        return ''
    if not FILTER_VALUE.match(value):
        raise ValidationError('%s is not a value this table can filter on.' % name)
    return value


def parse_view_query(query):
    """Validates the view half of /designs into a TableView.

    Assumes _checked_keys has already refused unknown names; this only reads the ones
    it knows. Repeats are refused the same way design parameters are, through _one.
    """
    sort = _one(query, 'sort') or ''
    if sort and sort not in SORT_COLUMNS:
        raise ValidationError('sort must name a column of the table.')

    page = _one(query, 'page')
    page = _int(page, 'page', 1, MAX_PAGE) if page is not None else 1

    return TableView(
        mutation=_open_value(query, 'mutation'),
        significance=_open_value(query, 'significance'),
        deaminase=_choice(query, 'deaminase', ALL_EDITS),
        strand=_choice(query, 'strand', STRANDS),
        hide_bsmbi=_bool(query, 'hide_bsmbi'),
        hide_4t=_bool(query, 'hide_4t'),
        sort=sort,
        dir=_choice(query, 'dir', DIRECTIONS) or 'asc',
        page=page,
    )


def parse_download_query(query, transcript_exists=None):
    """Validates /designs/download into (transcript_id, DesignParams, file name).

    Filters are not accepted: a download is the complete result, byte-identical to
    what the CLI writes, so a filtered view cannot silently hand someone a partial
    file that looks like the whole one.

    `file` is matched against RESULT_FILES rather than turned into a name, which is
    what keeps it out of the storage path.
    """
    _check_keys(query, DOWNLOAD_PARAMS)
    transcript = _transcript(query, transcript_exists)
    name = _one(query, 'file') or 'designs'
    if name not in RESULT_FILES:
        raise ValidationError('file must be one of %s.'
                              % ', '.join(sorted(RESULT_FILES)))
    return transcript, parse_editor_params(query), name
