"""The parent process's read-only view of the reference data.

The app needs two things from the bundle before it schedules anything: whether a
transcript exists, and the release and ClinVar versions that go into the cache key.
Worker processes open their own sources; nothing here is shared across processes.

Connections are per-thread. The ASGI server answers requests on a thread pool, and a
single sqlite3 connection used concurrently interleaves its cursors: one thread steps
another's statement, `fetchone` comes back empty, and a transcript that is in the
bundle is reported as missing. Read-only connections cost about 0.2 ms to open and
the pool is bounded, so one per thread is the cheap fix rather than a lock that would
serialise every lookup.
"""
import os
import re
import sqlite3
import threading

from bedesign.transcript_source import (TRANSCRIPTS_DB, find_bundle,
                                        find_clinvar_db)

CLINVAR_VERSION = re.compile(r'clinvar-(.+)\.db$')


class References(object):
    def __init__(self, refdata='refdata', clinvar_db=''):
        self.bundle = find_bundle(refdata)
        self.clinvar_db = find_clinvar_db(clinvar_db or None, refdata)
        self._uri = 'file:%s?mode=ro' % os.path.join(self.bundle, TRANSCRIPTS_DB)
        self._local = threading.local()
        meta = dict(self._db().execute('SELECT key, value FROM meta').fetchall())
        self.release = meta.get('release', 'unknown')
        match = CLINVAR_VERSION.search(os.path.basename(self.clinvar_db))
        self.clinvar_version = match.group(1) if match else 'unknown'

    def _db(self):
        """This thread's connection, opened on first use."""
        db = getattr(self._local, 'db', None)
        if db is None:
            db = sqlite3.connect(self._uri, uri=True)
            self._local.db = db
        return db

    def transcript_exists(self, transcript_id):
        """Parameterized, and only ever called with an ID that already matched
        the transcript pattern."""
        row = self._db().execute(
            'SELECT 1 FROM transcript WHERE transcript_id = ? LIMIT 1',
            (transcript_id,)).fetchone()
        return row is not None
