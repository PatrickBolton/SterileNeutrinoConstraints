"""Regenerate plots/ by executing the notebooks headlessly.

The notebooks are executed in memory and are not modified on disk.

Usage:
    uv run python scripts/make_plots.py              # all notebooks
    uv run python scripts/make_plots.py 4_summary    # selected notebooks
"""

import os
import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient

REPO_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = [
    "0_electron_mixing",
    "1_muon_mixing",
    "2_tau_mixing",
    "4_summary",
    "5_dark_matter",
]


def run_notebook(name):
    path = REPO_ROOT / f"{name}.ipynb"
    nb = nbformat.read(path, as_version=4)
    client = NotebookClient(
        nb,
        timeout=1800,
        kernel_name="python3",
        resources={"metadata": {"path": str(REPO_ROOT)}},
    )
    start = time.perf_counter()
    client.execute()
    print(f"{name}: done in {time.perf_counter() - start:.1f} s")


def main(argv):
    # Fixed timestamp so regenerated PDFs are byte-stable
    os.environ.setdefault("SOURCE_DATE_EPOCH", "0")
    for name in argv or NOTEBOOKS:
        run_notebook(name.removesuffix(".ipynb"))


if __name__ == "__main__":
    main(sys.argv[1:])
