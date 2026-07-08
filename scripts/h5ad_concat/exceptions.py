from __future__ import annotations

from h5ad_concat.models import SkipReason


class FileRejected(Exception):
    """Raised when a single h5ad fails validation and must be excluded from concat."""

    def __init__(self, reason: SkipReason) -> None:
        self.reason = reason
        super().__init__(reason.value)
