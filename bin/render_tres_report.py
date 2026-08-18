#!/usr/bin/env python3
"""Pipeline report entry point backed by the shared TrESFlow QC package."""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))

from tresflow_qc.cli import main  # noqa: E402
from tresflow_qc.core import read_duplicate_metrics  # noqa: E402,F401


if __name__ == "__main__":
    raise SystemExit(main())
