"""Pixel-compare plots/*.png against a baseline directory.

Usage:
    uv run python scripts/compare_plots.py --save      # copy plots/*.png into the baseline
    uv run python scripts/compare_plots.py             # compare plots/*.png with the baseline
    uv run python scripts/compare_plots.py --tol 0.5   # allow a small RMS difference

Exits non-zero if any image differs, is missing, or is new. Diff images are written
next to the regenerated plots as *-failed-diff.png (git-ignored).
"""

import argparse
import shutil
import sys
from pathlib import Path

from matplotlib.testing.compare import compare_images
from matplotlib.testing.exceptions import ImageComparisonFailure

REPO_ROOT = Path(__file__).resolve().parents[1]
PLOTS_DIR = REPO_ROOT / "plots"
BASELINE_DIR = REPO_ROOT / "tests" / "baseline"


def save_baseline(plots_dir, baseline_dir):
    baseline_dir.mkdir(parents=True, exist_ok=True)
    for png in sorted(plots_dir.glob("*.png")):
        if png.name.endswith("-failed-diff.png"):
            continue
        shutil.copy2(png, baseline_dir / png.name)
        print(f"saved {png.name}")


def compare(plots_dir, baseline_dir, tol):
    new = {p.name for p in plots_dir.glob("*.png") if not p.name.endswith("-failed-diff.png")}
    old = {p.name for p in baseline_dir.glob("*.png")}
    if not old:
        print(f"No baseline images in {baseline_dir}; run with --save first.")
        return 1
    failures = 0
    for name in sorted(old | new):
        if name not in new:
            print(f"MISSING  {name}")
            failures += 1
            continue
        if name not in old:
            print(f"NEW      {name}")
            failures += 1
            continue
        try:
            result = compare_images(str(baseline_dir / name), str(plots_dir / name), tol=tol, in_decorator=True)
        except ImageComparisonFailure as err:
            print(f"DIFFERS  {name}  ({err})")
            failures += 1
            continue
        if result is None:
            print(f"OK       {name}")
        else:
            print(f"DIFFERS  {name}  (rms={result['rms']:.4f})")
            failures += 1
    print(f"\n{len(old | new) - failures}/{len(old | new)} identical (tol={tol})")
    return 1 if failures else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--save", action="store_true", help="save plots/*.png as the new baseline")
    parser.add_argument("--tol", type=float, default=0.0, help="RMS tolerance (default 0: pixel-identical)")
    parser.add_argument("--plots", type=Path, default=PLOTS_DIR)
    parser.add_argument("--baseline", type=Path, default=BASELINE_DIR)
    args = parser.parse_args()
    if args.save:
        save_baseline(args.plots, args.baseline)
        return 0
    return compare(args.plots, args.baseline, args.tol)


if __name__ == "__main__":
    sys.exit(main())
