"""Where finished results live.

A small interface, so the same app could run against a bucket as well as the
directory it uses now (design doc 6). Keys look like object-store keys --
'results/<key[:2]>/<key>/designs.tsv.gz' -- and the local backend maps them onto a
directory.

The results directory is in memory on Cloud Run, so `evict` keeps it under a size
budget: without it, a client could fill memory with distinct results.

Keys are built by the app from a sha256 digest and a fixed filename, never from user
input directly. `_safe_segments` enforces that anyway: a backend that silently
resolved '..' would turn any future key-handling mistake into an arbitrary write.
"""
import itertools
import os
import re
import shutil
import threading

# Appended to a group's directory name while `evict` deletes it. No key the app
# builds contains it: result directories are named by hex digest.
EVICTING = '.evicting-'

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
        # SEGMENT requires an alphanumeric first character, which is also what
        # rules out '.' and '..'.
        if not SEGMENT.match(segment):
            raise BadKey('unsafe segment %r in storage key %r' % (segment, key))
    return segments


class Storage(object):
    """What the app needs from a result store: read, write, and keep it bounded."""

    def exists(self, key):
        raise NotImplementedError

    def get(self, key):
        """The object's bytes. Raises KeyError if it is not there."""
        raise NotImplementedError

    def put(self, key, data):
        raise NotImplementedError

    def touch(self, key):
        """Marks an object as just used, for `evict`."""
        raise NotImplementedError

    def evict(self, prefix, marker, max_bytes):
        """Deletes least recently used groups under `prefix` down to `max_bytes`."""
        raise NotImplementedError


class LocalStorage(Storage):
    """A directory on disk. The M1 backend, and what `docker compose` uses at M2."""

    def __init__(self, root):
        os.makedirs(root, exist_ok=True)
        # realpath, not abspath: on macOS the temp directory is itself a symlink
        # (/var -> /private/var), and a root that is not already resolved would
        # make every key look like it escapes.
        self.root = os.path.realpath(root)
        self._evict_lock = threading.Lock()
        self._evictions = itertools.count()

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
        except FileNotFoundError:
            raise KeyError(key)

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

    def touch(self, key):
        """Marks an object as just used, which is what `evict` orders by. An object
        that is not there is ignored: it may have been evicted a moment ago."""
        try:
            os.utime(self.path(key))
        except FileNotFoundError:
            pass

    def evict(self, prefix, marker, max_bytes):
        """Deletes the least recently used groups under `prefix` until the rest fit
        in `max_bytes`. Returns how many it deleted.

        A group is a directory holding `marker` -- for results, one key's files and
        its manifest. A directory without the marker is a job still writing, so it
        is neither counted nor touched. Age is the marker's mtime, which `touch`
        refreshes on use.

        Each group is renamed out of the way before it is deleted. The rename is
        atomic, so a reader sees the whole group or none of it, and a job that
        starts recomputing the key meanwhile writes to a fresh directory rather than
        into one being deleted. A rename left behind by a crash is cleared next time.
        """
        base = self.path(prefix)
        with self._evict_lock:
            groups = []
            total = 0
            for directory, subdirs, files in os.walk(base):
                for name in [d for d in subdirs if EVICTING in d]:
                    shutil.rmtree(os.path.join(directory, name), ignore_errors=True)
                    subdirs.remove(name)
                if marker not in files:
                    continue
                size = 0
                for name in files:
                    try:
                        size += os.path.getsize(os.path.join(directory, name))
                    except FileNotFoundError:
                        pass
                try:
                    used = os.path.getmtime(os.path.join(directory, marker))
                except FileNotFoundError:
                    continue
                groups.append((used, size, directory))
                total += size

            evicted = 0
            for _, size, directory in sorted(groups):
                if total <= max_bytes:
                    break
                doomed = '%s%s%d-%d' % (directory, EVICTING, os.getpid(),
                                        next(self._evictions))
                try:
                    os.rename(directory, doomed)
                except FileNotFoundError:
                    continue
                shutil.rmtree(doomed, ignore_errors=True)
                total -= size
                evicted += 1
            return evicted
