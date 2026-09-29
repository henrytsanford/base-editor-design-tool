"""The TSV that design rows are written as, in one place.

The CLI, the service worker and the MCP server all write the same three files. This
module owns the dialect so a change to it cannot reach one of them and not the
others; the golden fixtures would catch a drift, but only after the fact.

csv.writer terminates rows with CR LF, which is what the CLI leaves on disk and what
`as_text` translates back to a bare LF -- the form test/golden stores. Do not
"correct" the terminator here: the goldens decode universally and would keep passing
while every file on disk quietly changed.
"""
import csv
import gzip
import io

DELIMITER = '\t'


def writer(fileobj):
    """A csv.writer in the one dialect this repository writes."""
    return csv.writer(fileobj, delimiter=DELIMITER)


def tsv_bytes(columns, rows):
    """A header and its rows, as the bytes the CLI writes to disk.

    Rows are written as given. clinvar_annotations rows are legitimately ragged --
    28 fields on a SNP match, truncated after 'Mutation category' otherwise -- and
    are never padded out to the header.
    """
    buffer = io.StringIO(newline='')
    out = writer(buffer)
    out.writerow(columns)
    out.writerows(rows)
    return buffer.getvalue().encode('utf-8')


def gzipped(data):
    """`data` gzipped.

    mtime=0 for the same reason the golden fixtures use it: the bytes should depend
    on the designs, not on when they were computed.
    """
    raw = io.BytesIO()
    with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as gz:
        gz.write(data)
    return raw.getvalue()


def tsv_gz(columns, rows):
    """`tsv_bytes` gzipped."""
    return gzipped(tsv_bytes(columns, rows))


def as_text(data):
    """Bytes as text with universal newlines, the form test/golden stores.

    The same translation `open(path)` does, so a file read back from a CLI run and a
    payload built here compare equal.
    """
    return io.TextIOWrapper(io.BytesIO(data), encoding='utf-8', newline=None).read()
