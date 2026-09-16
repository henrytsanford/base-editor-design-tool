"""The parent process's read-only view of the reference data.

The app needs two things from the bundle before it schedules anything: whether a
transcript exists, and the release and ClinVar versions that go into the cache key.
Worker processes open their own sources; this one is never shared across processes.
"""
import os
import re
import sqlite3

from bedesign.transcript_source import (TRANSCRIPTS_DB, find_bundle,
                                        find_clinvar_db)

CLINVAR_VERSION = re.compile(r'clinvar-(.+)\.db$')


class References(object):
    def __init__(self, refdata='refdata', clinvar_db=''):
        self.bundle = find_bundle(refdata)
        self.clinvar_db = find_clinvar_db(clinvar_db or None, refdata)
        # check_same_thread=False because the ASGI server answers requests on a
        # thread pool. Reads are serialised by SQLite and the connection is
        # read-only, so there is no writer to race with.
        self._db = sqlite3.connect(
            'file:%s?mode=ro' % os.path.join(self.bundle, TRANSCRIPTS_DB),
            uri=True, check_same_thread=False)
        meta = dict(self._db.execute('SELECT key, value FROM meta').fetchall())
        self.release = meta.get('release', 'unknown')
        match = CLINVAR_VERSION.search(os.path.basename(self.clinvar_db))
        self.clinvar_version = match.group(1) if match else 'unknown'

    def transcript_exists(self, transcript_id):
        """Parameterized, and only ever called with an ID that already matched
        the transcript pattern."""
        row = self._db.execute(
            'SELECT 1 FROM transcript WHERE transcript_id = ? LIMIT 1',
            (transcript_id,)).fetchone()
        return row is not None
