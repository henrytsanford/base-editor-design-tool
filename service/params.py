"""Turning a query string into a validated design request.

Everything a request can say is checked here, before a job is scheduled, so a
malformed or hostile request costs a regex rather than a CPU core. The rules are the
design doc's 5, plus the preset allowlist and the caps that bound the work of
checking.

`ValidationError` carries a message written for the person who typed the URL. It
never quotes a filesystem path or an internal exception.
"""
import dataclasses
import re

from bedesign.engine import ALL_EDITS, DesignParams, UnknownBaseEditor

# A request has at most a handful of parameters; anything longer is a probe.
MAX_PARAMS = 16
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
DESIGN_PARAMS = ('transcript', 'preset') + tuple(
    f.name for f in dataclasses.fields(DesignParams))


class ValidationError(ValueError):
    """A request that will not be run, with a message safe to show the user."""


def _one(query, name):
    """The single value for a name, refusing repeats.

    Two values for one parameter is how a request tries to make the validator and
    the code that reads it disagree, so it is an error rather than last-wins.
    """
    values = query.getlist(name)
    if not values:
        return None
    if len(values) > 1:
        raise ValidationError('%s was given more than once.' % name)
    value = values[0].strip()
    if len(value) > MAX_VALUE_LEN:
        raise ValidationError('%s is too long.' % name)
    return value


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


def parse_designs_query(query, transcript_exists=None):
    """Validates /designs parameters into (transcript_id, DesignParams).

    `transcript_exists` is called only after the ID matches the pattern, so an
    unparseable ID never reaches the database.
    """
    keys = set(query.keys())
    if len(keys) > MAX_PARAMS:
        raise ValidationError('Too many parameters.')
    unknown = sorted(keys - set(DESIGN_PARAMS))
    if unknown:
        # Rejected rather than ignored: silently dropping a parameter would let a
        # user believe a setting applied when the result ignored it.
        raise ValidationError('Unknown parameter: %s.' % ', '.join(unknown))

    transcript = _one(query, 'transcript')
    if not transcript:
        raise ValidationError('A transcript is required, e.g. ENST00000307102.')
    transcript = transcript.upper()
    if not TRANSCRIPT.match(transcript):
        raise ValidationError(
            'Transcript must look like ENST00000307102: ENST followed by 11 digits.')
    if transcript_exists is not None and not transcript_exists(transcript):
        raise ValidationError('Transcript %s is not in the reference bundle.' % transcript)

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
        conflicting = sorted(k for k in EDITOR_PARAMS if k in keys)
        if conflicting:
            raise ValidationError(
                'preset already sets %s; give one or the other, not both.'
                % ', '.join(conflicting))
        try:
            # The preset table in the engine is the allowlist, so the dropdown and
            # this check cannot drift apart.
            return transcript, DesignParams.from_preset(preset, **overrides)
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

    return transcript, DesignParams(**fields)
