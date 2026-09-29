"""Tool arguments, shaped so service.params can validate them.

The web app validates every request in service/params.py: the PAM pattern, the window
bounds, the sg_len and intron_buffer ranges, the preset conflict rule, the filter
vocabularies. Those rules are the contract the engine is called under, and there is no
version of them that is right for HTTP and wrong for MCP.

So rather than restating any of it, this wraps a dict of tool arguments in the part of
starlette's QueryParams that params.py actually uses -- `keys` and `getlist`, and
nothing else -- and the validators run unchanged. A rule added there reaches the tools
without anyone remembering to add it twice.
"""


def _text(value):
    """A tool argument as the string a query parameter would have carried.

    params._flag reads 'true' and 'false' case-insensitively, so a JSON boolean and a
    typed-in one arrive as the same answer.
    """
    if value is True:
        return 'true'
    if value is False:
        return 'false'
    return str(value)


class Query(object):
    """A read-only multidict over tool arguments.

    An argument left None was not given, and is dropped rather than passed as the
    string 'None'. That is what params._one reads as absent -- so an unset editor
    field is not a conflict with a preset, and an unset filter is not a filter.
    """

    def __init__(self, values):
        self._values = {name: _text(value) for name, value in values.items()
                        if value is not None}

    def keys(self):
        return self._values.keys()

    def getlist(self, name):
        value = self._values.get(name)
        return [] if value is None else [value]

