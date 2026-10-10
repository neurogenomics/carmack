"""Template identity and primary-mate validation shared by BAM deduplicators.

Coordinates describe duplicate groups, not the identity of a sequenced template.
QNAME is scoped by read group and cell; reuse within that namespace is ambiguous.
"""

import json
import logging
import sqlite3
import sys
from collections import OrderedDict
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import NamedTuple

import pysam

log = logging.getLogger(__name__)


class TemplateIdentity(NamedTuple):
    query_name: str
    read_group: str | None
    barcode: str


def template_identity(
    read: pysam.AlignedSegment, barcode: str | None = None, barcode_tag: str = "CB"
) -> TemplateIdentity:
    if not read.query_name:
        raise ValueError("Template identity requires a nonempty QNAME.")
    if barcode is None:
        if not read.has_tag(barcode_tag):
            raise ValueError(
                f"Read {read.query_name!r} does not have a barcode tag ({barcode_tag})."
            )
        barcode = read.get_tag(barcode_tag)
    if not isinstance(barcode, str) or not barcode:
        raise ValueError(f"Read {read.query_name!r} has an invalid cell barcode.")
    read_group = read.get_tag("RG") if read.has_tag("RG") else None
    if read_group is not None and not isinstance(read_group, str):
        raise ValueError(f"Read {read.query_name!r} has an invalid RG tag.")
    return TemplateIdentity(
        read.query_name,
        sys.intern(read_group) if read_group is not None else None,
        sys.intern(barcode),
    )


def _end(read: pysam.AlignedSegment, mate: bool = False) -> tuple | None:
    if read.mate_is_unmapped if mate else read.is_unmapped:
        # SAM permits unmapped records to carry a mapped mate's coordinates.
        return None
    return (
        (read.next_reference_id, read.next_reference_start, read.mate_is_reverse)
        if mate
        else (read.reference_id, read.reference_start, read.is_reverse)
    )


class PrimaryPairValidator:
    """Reject ambiguous/incomplete primary pairs without retaining BAM records.

    At most 100,000 outstanding mates are held in RAM; older entries spill to SQLite.
    A temporary UNIQUE index catches reused complete identities across chromosomes
    with an 8 MiB page cache. An unpaired third primary record fails at EOF.
    """

    def __init__(self, max_pending: int = 100_000) -> None:
        if max_pending < 1:
            raise ValueError("max_pending must be positive")
        self.max_pending = max_pending
        self.pending: OrderedDict[TemplateIdentity, tuple] = OrderedDict()
        self.spilled_count = 0
        self._temporary = TemporaryDirectory(prefix="carmack-template-index-")
        self.index_path = Path(self._temporary.name) / "completed.sqlite"
        self._index = sqlite3.connect(str(self.index_path))
        # Disposable validation scratch: no recovery journal or durable commit needed.
        self._index.execute("PRAGMA journal_mode=OFF")
        self._index.execute("PRAGMA synchronous=OFF")
        self._index.execute("PRAGMA cache_size=-8192")
        self._index.execute("CREATE TABLE completed (identity TEXT PRIMARY KEY) WITHOUT ROWID")
        self._index.execute(
            "CREATE TABLE pending (identity TEXT PRIMARY KEY, signature TEXT) WITHOUT ROWID"
        )
        self.completed_count = 0

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self._index.close()
        self._temporary.cleanup()

    def observe(self, read: pysam.AlignedSegment, identity: TemplateIdentity) -> TemplateIdentity:
        if not read.is_paired or read.is_read1 == read.is_read2:
            raise ValueError(f"Invalid primary mate flags for template {identity!r}.")
        previous = self.pending.pop(identity, None)
        if previous is None and self.spilled_count:
            encoded = json.dumps(identity)
            row = self._index.execute(
                "SELECT signature FROM pending WHERE identity = ?", (encoded,)
            ).fetchone()
            if row is not None:
                self._index.execute("DELETE FROM pending WHERE identity = ?", (encoded,))
                self.spilled_count -= 1
                role, own, mate = json.loads(row[0])
                previous = (
                    identity,
                    role,
                    tuple(own) if own else None,
                    tuple(mate) if mate else None,
                )
        current = (identity, read.is_read1, _end(read), _end(read, mate=True))
        if previous is None:
            if len(self.pending) >= self.max_pending:
                old_identity, old_signature = self.pending.popitem(last=False)
                self._index.execute(
                    "INSERT INTO pending VALUES (?, ?)",
                    (json.dumps(old_identity), json.dumps(old_signature[1:])),
                )
                self.spilled_count += 1
            self.pending[identity] = current
            return identity
        elif previous[1] == current[1]:
            raise ValueError(f"Ambiguous template identity {identity!r}: repeated primary mate.")
        elif previous[2] != current[3] or previous[3] != current[2]:
            raise ValueError(
                f"Non-reciprocal primary mate coordinates/flags for template {identity!r}."
            )
        else:
            # Full lossless identity, not a hash. One insert per complete pair;
            # there is no completed-index query for each incoming alignment.
            try:
                self._index.execute("INSERT INTO completed VALUES (?)", (json.dumps(identity),))
            except sqlite3.IntegrityError as error:
                raise ValueError(
                    f"Ambiguous template identity {identity!r}: repeated primary pair."
                ) from error
            self.completed_count += 1
            return previous[0]

    def finish(self) -> int:
        pending_count = len(self.pending) + self.spilled_count
        if pending_count:
            example = (
                next(iter(self.pending))
                if self.pending
                else TemplateIdentity(
                    *json.loads(
                        self._index.execute("SELECT identity FROM pending LIMIT 1").fetchone()[0]
                    )
                )
            )
            raise ValueError(
                f"Incomplete primary pairs: {pending_count} template identities lack a mate; "
                f"example {example!r}. Check barcode/RG consistency and input completeness."
            )
        self._index.commit()
        log.info(
            "Validated %s primary template identities; temporary collision index %s bytes "
            "(8 MiB page cache).",
            self.completed_count,
            self.index_path.stat().st_size,
        )
        return self.completed_count
