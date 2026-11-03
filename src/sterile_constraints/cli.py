"""Command-line entry points."""

import argparse
import os
import time

from matplotlib import pyplot as plt


def make_plots(argv=None):
    """Regenerate plots/ from the figure specifications (`uv run make-plots [name ...]`)."""
    parser = argparse.ArgumentParser(description=make_plots.__doc__)
    parser.add_argument("names", nargs="*", help="figure names (default: all), e.g. UeNsq_constraints")
    args = parser.parse_args(argv)

    # Fixed timestamp so regenerated PDFs are byte-stable
    os.environ.setdefault("SOURCE_DATE_EPOCH", "0")

    from .figures import FIGURES
    from .plotting import plot_figure, save_figure

    unknown = set(args.names) - set(FIGURES)
    if unknown:
        parser.error(f"unknown figure(s) {sorted(unknown)}; available: {sorted(FIGURES)}")
    for name in args.names or FIGURES:
        start = time.perf_counter()
        fig, _ = plot_figure(FIGURES[name])
        save_figure(fig, name)
        plt.close(fig)
        print(f"{name}: {time.perf_counter() - start:.1f} s")
