"""Constraint registry (catalog/constraints.toml) and loading of curves in canonical units.

Each constraint is identified by its path under data/ without the .csv suffix, e.g.
'UeNsq/current/KATRIN'. The CSV header is the source of truth for the description,
references and units; the registry adds what the header cannot say (category, default
plot label and style, unit-conversion choices). Computed constraints (0vbb) have a
`computed` table instead of a CSV.

    >>> from sterile_constraints.catalog import load_curve
    >>> mN, Usq = load_curve("UeNsq/current/KATRIN")   # mN in GeV, |U_eN|^2
"""

import re
import tomllib
from dataclasses import dataclass, field
from functools import cache

import numpy as np
import pandas as pd

from .paths import DATA_DIR, REPO_ROOT
from .physics import N_0vbb_bound
from .units import parse_units

REGISTRY_PATH = REPO_ROOT / "catalog" / "constraints.toml"

FLAVOURS = {"UeNsq": "e", "UmuNsq": "mu", "UtauNsq": "tau", "UalphaNsq": "alpha"}

# CSVs still waiting for a header (description, source, units line); skipped with a warning
PENDING_HEADERS = frozenset({
    "UtauNsq/future/MesonDecays_LNV_tau",
    "UtauNsq/future/HL_LHC_tau",
    "UtauNsq/future/CODEX_b_tau",
    "UmuNsq/future/ATLAS_1_mu",
    "UmuNsq/future/ATLAS_2_mu",
    "UmuNsq/future/ATLAS_3_mu",
    "UmuNsq/future/CODEX_b_mu",
})


class PendingHeaderError(ValueError):
    pass


@dataclass(frozen=True)
class Constraint:
    id: str
    category: str
    label: str | None = None
    style: dict = field(default_factory=dict)
    y_conversion: str = "exact"
    computed: dict | None = None
    note: str = ""

    @property
    def flavour(self):
        return FLAVOURS[self.id.split("/")[0]]

    @property
    def status(self):
        return self.id.split("/")[1]

    @property
    def path(self):
        return None if self.computed else DATA_DIR / f"{self.id}.csv"

    @property
    def header(self):
        """The CSV's '#' comment lines, without the leading '# '."""
        if self.computed:
            return []
        with open(self.path, encoding="utf-8") as f:
            return [line[1:].strip() for line in f if line.startswith("#")]

    @property
    def description(self):
        return self.header[0] if self.header else ""

    @property
    def references(self):
        """Header lines between the description and the units line (figure, arXiv URL, ...)."""
        return self.header[1:-1]

    @property
    def arxiv_ids(self):
        return re.findall(r"arxiv\.org/(?:pdf|abs)/([\w.\-/]+?)(?:v\d+)?(?:\.pdf)?(?=\s|$)", " ".join(self.references))

    @property
    def units(self):
        if self.id in PENDING_HEADERS:
            raise PendingHeaderError(f"{self.id}: CSV has no header yet")
        return parse_units(self.header[-1], self.y_conversion)


@cache
def load_registry(path=REGISTRY_PATH):
    """Return {id: Constraint} for every entry of the registry."""
    with open(path, "rb") as f:
        raw = tomllib.load(f)
    return {cid: Constraint(id=cid, **entry) for cid, entry in raw.items()}


def get(cid):
    try:
        return load_registry()[cid]
    except KeyError:
        raise KeyError(f"{cid!r} is not in {REGISTRY_PATH.relative_to(REPO_ROOT)}") from None


def read_columns(path):
    """Read the two data columns of a constraint CSV as float arrays."""
    df = pd.read_csv(path, sep=",", names=["X", "Y"], comment="#")
    return df["X"].to_numpy(dtype=float), df["Y"].to_numpy(dtype=float)


def load_curve(constraint, mN=None):
    """Return (mN [GeV], |U|^2) arrays for a constraint id or Constraint.

    For computed constraints, mN is the mass grid to evaluate on (required).
    """
    c = constraint if isinstance(constraint, Constraint) else get(constraint)
    if c.computed:
        if mN is None:
            raise ValueError(f"{c.id} is computed: pass a mass grid mN")
        if c.computed["model"] != "0vbb":
            raise ValueError(f"{c.id}: unknown computed model {c.computed['model']!r}")
        Usq = np.array([N_0vbb_bound(m, c.computed["thalf"], c.computed["isotope"]) for m in mN])
        return np.asarray(mN), Usq
    return c.units.convert(read_columns(c.path))
