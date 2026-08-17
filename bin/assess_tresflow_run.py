#!/usr/bin/env python3
"""Standalone TrESFlow QC assessor backed by the shared report package."""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))

from tresflow_qc.cli import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
