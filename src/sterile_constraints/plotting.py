"""Draw constraint figures from a `Figure` specification (see figures.py).

A figure has four ordered groups of items, drawn in this order:

    xray     current bounds filled down to the τ_N = t_U line (capped at 1e-5)
    current  current bounds filled up to the top of the plot
    cosmo    like `current`, but still drawn in colour on 'future' figures
    future   future sensitivities, drawn as lines on 'future' figures only

On a 'future' figure, xray and current bounds are drawn as a single solid `block_color`
region without labels.
"""

from dataclasses import dataclass, field
from functools import cache

import matplotlib.cm as cm
import matplotlib.ticker as ticker
import numpy as np
from matplotlib import pyplot as plt

from .catalog import get, load_curve
from .paths import PLOTS_DIR
from .physics import N_DarkMatter_lifetime
from .style import apply_style

apply_style()

YEAR = 3600 * 24 * 365
T_UNIVERSE = YEAR * 13.8e9  # [s]
T_CMB = YEAR * 3.78e5  # [s]

# Mass grids [GeV] for computed (0vbb) curves
CURRENT_GRID = np.logspace(-15.1, 5.1, 1000)
FUTURE_GRID = np.logspace(-15.1, 10.1, 1000)

LABEL_FONTSIZE = 13


@dataclass(frozen=True)
class Curve:
    """A registry constraint drawn on a figure.

    label_at is the (log10 x, log10 y) position of the label in the figure's axes.
    label and the style keywords override the registry defaults for this figure.
    """

    id: str
    label_at: tuple | None = None
    rot: float = 0
    show_label: bool = True
    label: str | None = None
    color: str | None = None
    lw: float | None = None
    ls: str | None = None
    fill: float | None = None


@dataclass(frozen=True)
class Text:
    """A free-standing annotation at (log10 x, log10 y)."""

    text: str
    at: tuple
    rot: float = 0
    color: str = "k"
    show: bool = True


@dataclass(frozen=True)
class Figure:
    name: str  # output file stem in plots/
    mode: str = "current"  # 'current' or 'future'
    x_min: float = 1e-14
    x_max: float = 1e5
    y_min: float = 1e-15
    y_max: float = 1
    n_x_ticks: int = 100
    n_y_ticks: int = 100
    shape: str = "wide"  # 'wide' or 'tall'
    xlabel: str = r"$m_N~[\mathrm{GeV}]$"
    ylabel: str | None = None
    x_scale: float = 1  # tick labels show value * scale (e.g. 1e6 for keV)
    y_scale: float = 1
    x_rot: float = 0
    flip: bool = False  # mN on the y axis
    block_color: str = "lightgrey"
    seesaw: bool = True
    light_nu_masses: bool = True
    dm_window: bool = True
    dm_lifetime: bool = False
    xray: list = field(default_factory=list)
    current: list = field(default_factory=list)
    cosmo: list = field(default_factory=list)
    future: list = field(default_factory=list)


class Canvas:
    """Thin wrapper around an Axes taking (mN, |U|^2) data, whichever way round the axes are."""

    def __init__(self, ax, flip):
        self.ax = ax
        self.flip = flip

    def line(self, mN, Usq, **kwargs):
        return self.ax.loglog(Usq, mN, **kwargs) if self.flip else self.ax.loglog(mN, Usq, **kwargs)

    def fill(self, mN, U1, U2=0, **kwargs):
        if self.flip:
            return self.ax.fill_betweenx(mN, U1, U2, **kwargs)
        return self.ax.fill_between(mN, U1, U2, **kwargs)

    def text(self, log_x, log_y, s, rot=0, **kwargs):
        """Text at (10^log_x, 10^log_y) in display coordinates (not swapped by flip)."""
        return self.ax.text(10**log_x, 10**log_y, s, fontsize=LABEL_FONTSIZE, rotation=rot, **kwargs)


# Light neutrino masses (normal ordering, m1 = 0) [GeV]
M1 = 0
M2 = np.sqrt(M1**2 + 7.537e-5) * 1e-9
M3 = np.sqrt(M2**2 + 2.521e-3) * 1e-9


@cache
def dm_lifetime_bound(mN, time):
    """|U_eN|^2 at which the HNL lifetime equals `time` [s]."""
    return N_DarkMatter_lifetime(mN, 1, 0, 0, time)


def add_seesaw_band(c):
    mN = np.logspace(-16, 10, 100)
    seesaw_bound_d = M2 / (M2 + mN)
    seesaw_bound_u = 0.1e-9 / (0.1e-9 + mN)
    nlevels = 30
    cols = cm.Greys(np.linspace(0.1, 0.45, nlevels))
    levels = (np.linspace(1, np.sqrt(13), nlevels)) ** 2
    for i in range(nlevels - 1):
        c.fill(mN, np.maximum(seesaw_bound_u / levels[i], seesaw_bound_d),
               np.minimum(seesaw_bound_d * levels[i], seesaw_bound_u),
               alpha=0.025, color=cols[i, :], zorder=-1000, lw=0)
    at = (-13.6, 4.) if c.flip else (3.5, -13.2)
    c.ax.text(10**at[0], 10**at[1], r'$\mathrm{Seesaw}$', fontsize=16, rotation=0, color='black')


def add_light_nu_masses(c):
    c.line([M2, M2], [1e-40, 10], lw=1, ls='--', color='grey')
    c.line([M3, M3], [1e-40, 10], lw=1, ls='--', color='grey')
    m2_label = r'$m_2^{\mathrm{min}} = \sqrt{\Delta m_\odot^2}$'
    m3_label = r'$m_3^{\mathrm{min}} = \sqrt{\Delta m_\odot^2 + \Delta m_{\mathrm{atm}}^2}$'
    if c.flip:
        c.text(-27.6, -10.85, m2_label, rot=0, color='k')
        c.text(-27.6, -10., m3_label, rot=0, color='k')
    else:
        c.text(-11.5, -14.6, m2_label, rot=90, color='k')
        c.text(-10.73, -14.6, m3_label, rot=90, color='k')


def add_dm_window(c):
    c.line([1e-8, 1e-8, 2e-3, 2e-3], [1e-17, 1e-5, 1e-5, 1e-17], lw=0.5, ls='--', color='grey')
    mvals = np.logspace(-4.695, np.log10(2e-3), 100)
    c.line(mvals, [dm_lifetime_bound(m, T_UNIVERSE) for m in mvals], lw=1.0, ls='--', color='k')
    c.text(-4.2, -8.4, r'$\tau_N < t_U$', rot=-70, color='k')
    c.text(-8.35, -14.6, r'$N$ Dark Matter', rot=90, color='k')


def add_dm_lifetime_lines(c):
    """τ_N = t_U, t_CMB and 1 s."""
    mvals = np.logspace(-10, 5, 2000)
    for time in (T_UNIVERSE, T_CMB, 1):
        c.line(mvals, [dm_lifetime_bound(m, time) for m in mvals], lw=1.0, ls='--', color='k')


def _style(item, constraint):
    style = dict(constraint.style)
    for key in ("color", "lw", "ls", "fill"):
        if getattr(item, key) is not None:
            style[key] = getattr(item, key)
    return style


def _draw_label(c, item, constraint, color):
    if item.label_at is not None and item.show_label:
        label = item.label if item.label is not None else constraint.label
        c.text(*item.label_at, label, rot=item.rot, color=color)


def draw_item(c, spec, group, item):
    """Draw one Curve or Text of a figure group."""
    current = spec.mode == "current"
    if isinstance(item, Text):
        visible = group == "cosmo" or (group == "future") != current
        if visible and item.show:
            c.text(*item.at, item.text, rot=item.rot, color=item.color)
        return

    if group == "future" and current:
        return
    constraint = get(item.id)
    mN, Usq = load_curve(constraint, mN=FUTURE_GRID if group == "future" else CURRENT_GRID)
    style = _style(item, constraint)
    color = style["color"]

    if group == "future":
        c.line(mN, Usq, lw=style["lw"], ls=style["ls"], color=color)
        _draw_label(c, item, constraint, color)
        return

    if group == "xray":
        lower = np.minimum([dm_lifetime_bound(m, T_UNIVERSE) for m in mN], 1e-5)
        fill_args = (mN, lower, Usq)
    else:
        fill_args = (mN, 10, Usq)

    if current or group == "cosmo":
        c.line(mN, Usq, lw=style["lw"], ls=style["ls"], color=color)
        c.fill(*fill_args, facecolor=color, alpha=style["fill"])
        _draw_label(c, item, constraint, color)
    else:
        c.fill(*fill_args, facecolor=spec.block_color, alpha=1)


class ScaledLogFormatterMathtext(ticker.LogFormatterMathtext):
    """Log tick labels showing value * scale (e.g. GeV data labelled in keV)."""

    def __init__(self, scale=1.0, **kwargs):
        super().__init__(**kwargs)
        self.scale = scale

    def __call__(self, x, pos=None):
        return super().__call__(x * self.scale, pos)


def format_axes(fig, ax, spec):
    ax.tick_params(axis='x', which='major', labelsize=26, pad=7.5, labelrotation=spec.x_rot)
    ax.tick_params(axis='y', which='major', labelsize=26)
    ax.set_xlim(spec.x_min, spec.x_max)
    ax.set_ylim(spec.y_min, spec.y_max)

    for axis, n_ticks, scale, set_label, label, rot in (
            (ax.xaxis, spec.n_x_ticks, spec.x_scale, ax.set_xlabel, spec.xlabel, 0),
            (ax.yaxis, spec.n_y_ticks, spec.y_scale, ax.set_ylabel, spec.ylabel, 90)):
        axis.set_major_locator(ticker.LogLocator(base=10.0, numticks=n_ticks))
        axis.set_major_formatter(ScaledLogFormatterMathtext(scale=scale, base=10.0))
        axis.set_minor_locator(ticker.LogLocator(base=10.0, subs=[i/10 for i in range(1, 10)], numticks=100))
        axis.set_minor_formatter(ticker.NullFormatter())
        set_label(label, fontsize=30, rotation=rot)

    # Ticks on all four sides, pointing inward
    ax.tick_params(axis='both', which='major', top=True, right=True, direction='in', length=8)
    ax.tick_params(axis='both', which='minor', top=True, right=True, direction='in', length=4)

    x0, x1 = np.log10(ax.get_xlim())
    y0, y1 = np.log10(ax.get_ylim())
    # Axes box aspect: golden ratio, landscape ('wide') or portrait ('tall')
    if spec.shape == 'wide':
        ratio = (1 + np.sqrt(5)) / 2
        ax.set_aspect(((x1 - x0) / (y1 - y0)) / ratio)
        fig.set_size_inches(15, 15)
    elif spec.shape == 'tall':
        ratio = 2 / (1 + np.sqrt(5))
        ax.set_aspect(((x1 - x0) / (y1 - y0)) / ratio)
        fig.set_size_inches(10, 10)


def plot_figure(spec):
    """Draw a Figure spec; returns (fig, ax). Saving is left to the caller (see save_figure)."""
    fig, ax = plt.subplots()
    c = Canvas(ax, spec.flip)

    if spec.seesaw:
        add_seesaw_band(c)
    if spec.light_nu_masses:
        add_light_nu_masses(c)
    if spec.dm_window:
        add_dm_window(c)
    if spec.dm_lifetime:
        add_dm_lifetime_lines(c)

    for group in ("xray", "current", "cosmo", "future"):
        for item in getattr(spec, group):
            draw_item(c, spec, group, item)

    format_axes(fig, ax, spec)
    return fig, ax


def save_figure(fig, name, plots_dir=PLOTS_DIR, formats=("pdf", "png")):
    for fmt in formats:
        fig.savefig(plots_dir / f"{name}.{fmt}", bbox_inches='tight')
