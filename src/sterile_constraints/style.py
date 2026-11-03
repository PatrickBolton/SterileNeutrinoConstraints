"""Matplotlib text rendering: LaTeX (as used for the published plots) or mathtext fallback.

The published plots are rendered with LaTeX (``text.usetex``), which needs ``latex`` and
``dvipng`` on the PATH. Control it with the ``STERILE_CONSTRAINTS_USETEX`` environment
variable, or by calling ``apply_style(use_tex=...)``:

    1 / true   always use LaTeX (error if it cannot be found)
    0 / false  never use LaTeX (matplotlib mathtext, Computer Modern math font)
    unset      use LaTeX if found, otherwise fall back to mathtext with a warning
"""

import glob
import os
import shutil
import warnings

from matplotlib import rc, rcParams

ENV_VAR = "STERILE_CONSTRAINTS_USETEX"
TEX_PROGRAMS = ("latex", "dvipng")

# Standard install locations, for when the PATH of a GUI-launched editor misses them
TEX_BIN_DIRS = (
    "/Library/TeX/texbin",
    *sorted(glob.glob("/usr/local/texlive/*/bin/*"), reverse=True),
)


def _tex_on_path():
    return all(shutil.which(p) for p in TEX_PROGRAMS)


def find_tex():
    """Return True if latex and dvipng are available, adding a standard TeX bin dir to PATH if needed."""
    if _tex_on_path():
        return True
    for d in TEX_BIN_DIRS:
        if all(os.access(os.path.join(d, p), os.X_OK) for p in TEX_PROGRAMS):
            os.environ["PATH"] += os.pathsep + d
            return True
    return False


def _env_setting():
    value = os.environ.get(ENV_VAR, "").strip().lower()
    if value in ("", "auto"):
        return None
    if value in ("1", "true", "yes", "on"):
        return True
    if value in ("0", "false", "no", "off"):
        return False
    raise ValueError(f"{ENV_VAR}={value!r}: expected 1/true, 0/false or unset")


def apply_style(use_tex=None):
    """Set the global matplotlib text rendering. Returns whether LaTeX is used."""
    if use_tex is None:
        use_tex = _env_setting()
    if use_tex is None:
        use_tex = find_tex()
        if not use_tex:
            warnings.warn(
                "LaTeX (latex + dvipng) not found: falling back to matplotlib mathtext. "
                f"Plots will differ slightly from the published ones. Set {ENV_VAR}=0 to silence.",
                stacklevel=2,
            )
    elif use_tex and not find_tex():
        raise RuntimeError(f"{ENV_VAR} requests LaTeX, but latex/dvipng were not found on PATH or in {TEX_BIN_DIRS}")

    if use_tex:
        rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
        rc('text', usetex=True)
    else:
        rcParams.update({
            'font.family': 'serif',
            'font.serif': ['cmr10', 'DejaVu Serif'],
            'text.usetex': False,
            'mathtext.fontset': 'cm',
            'axes.unicode_minus': False,
        })
    return use_tex
