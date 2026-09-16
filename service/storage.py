"""Where finished results live.

Two methods, so the same app runs against a directory locally and against S3 at M2
(design doc 6). Keys look like S3 keys -- 'results/<key>/designs.tsv.gz' -- and the
local backend maps them onto a directory.

Keys are built by the app from a sha256 digest and a fixed filename, never from user
input directly. `_safe_segments` enforces that anyway: a backend that silently
resolved '..' would turn any future key-handling mistake into an arbitrary write.
"""
import errno
import os
import re

# Deliberately narrow. Every key this app builds is hex digits, ASCII words, dots
# and dashes, so anything else means a bug upstream rather than a new use case.
SEGMENT = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]*$')


class BadKey(ValueError):
    """A storage key that is not safe to turn into a path."""


def _safe_segments(key):
    """Splits a key into path segments, refusing anything that could escape the root."""
    if not key or key != key.strip():
        raise BadKey('empty or padded storage key %r' % key)
    segments = key.split('/')
    for segment in segments:
        if not SEGMENT.match(segment) or segment in ('.', '..'):
            raise BadKey('unsafe segment %r in storage key %r' % (segment, key))
    return segments


class Storage(object):
    """The two things the app needs from a result store, plus a local-file escape hatch."""

    def exists(self, key):
        raise NotImplementedError

    def get(self, key):
        """The object's bytes. Raises KeyError if it is not there."""
        raise NotImplementedError

    def put(self, key, data):
        raise NotImplementedError


class LocalStorage(Storage):
    """A directory on disk. The M1 backend, and what `docker compose` uses at M2."""

    def __init__(self, root):
        os.makedirs(root, exist_ok=True)
        # realpath, not abspath: on macOS the temp directory is itself a symlink
        # (/var -> /private/var), and a root that is not already resolved would
        # make every key look like it escapes.
        self.root = os.path.realpath(root)

    def path(self, key):
        """The absolute path for a key, checked to be inside the root.

        The segment check already rules out traversal; resolving and re-checking
        catches the case the segment rules cannot see, a symlink inside the root
        pointing out of it.
        """
        path = os.path.join(self.root, *_safe_segments(key))
        resolved = os.path.realpath(path)
        if resolved != self.root and not resolved.startswith(self.root + os.sep):
            raise BadKey('storage key %r resolves outside the store' % key)
        return path

    def exists(self, key):
        return os.path.exists(self.path(key))

    def get(self, key):
        try:
            with open(self.path(key), 'rb') as fh:
                return fh.read()
        except IOError as e:
            if e.errno == errno.ENOENT:
                raise KeyError(key)
            raise

    def put(self, key, data):
        """Writes via a temporary file in the same directory, then renames.

        A reader that finds manifest.json must find a complete object behind it, and
        a job killed by the timeout mid-write must not leave a truncated one.
        """
        path = self.path(key)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = '%s.%d.tmp' % (path, os.getpid())
        with open(tmp, 'wb') as fh:
            fh.write(data)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
