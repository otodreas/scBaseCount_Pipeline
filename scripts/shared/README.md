# shared

Internal utilities shared across all pipeline packages.

## Modules

### `repo.py`

Locates the repo root by walking up from `__file__` until a `.venv` or `.git` directory is found.

```python
from shared.repo import REPO_ROOT
```

`REPO_ROOT` is a `pathlib.Path` resolved once at import time.

### `logger.py`

Configures a file logger rooted at `REPO_ROOT / "logs"` and returns a named `logging.Logger`.

```python
from shared.logger import configure_file_logger

_log = configure_file_logger("my_package.log", __name__)
```

Writes append-only to `logs/<log_filename>` at `INFO` level with format `YYYY-MM-DD HH:MM:SS LEVEL message`.
