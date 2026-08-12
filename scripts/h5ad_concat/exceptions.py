from h5ad_concat.models import QcStats, SkipReason


class FileRejected(Exception):
    """Raised when a single h5ad fails validation and must be excluded from concat."""

    def __init__(self, reason: SkipReason, qc: QcStats | None = None) -> None:
        self.reason = reason
        self.qc = qc
        super().__init__(reason.value)
