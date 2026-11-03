import warnings

import numpy as np
import pytest

from sterile_constraints.catalog import PENDING_HEADERS, load_curve, load_registry
from sterile_constraints.paths import DATA_DIR
from sterile_constraints.plotting import CURRENT_GRID
from sterile_constraints.units import Y_CONVERSIONS

CATEGORIES = {"oscillation", "beta_decay", "0vbb", "decay_in_flight", "peak_search", "fixed_target",
              "collider", "precision", "low_energy", "xray", "supernovae", "cosmology"}
CSV_IDS = sorted(str(p.relative_to(DATA_DIR).with_suffix("")) for p in DATA_DIR.glob("*/*/*.csv"))
REGISTRY = load_registry()

# Known inconsistencies awaiting a decision (units line vs directory)
KNOWN_FLAVOUR_MISMATCH = {
    "UalphaNsq/future/FCC_ee": "units line says |UeN|^2 but the file is used as flavour-universal",
}


def test_every_csv_is_registered():
    assert sorted(set(CSV_IDS) - set(REGISTRY)) == []


def test_every_registry_entry_has_data():
    missing = [cid for cid, c in REGISTRY.items() if not c.computed and not c.path.exists()]
    assert missing == []


def test_pending_headers_are_known_files():
    assert PENDING_HEADERS <= set(CSV_IDS)


@pytest.mark.parametrize("cid", sorted(REGISTRY))
def test_registry_entry(cid):
    c = REGISTRY[cid]
    assert c.category in CATEGORIES
    assert c.y_conversion in Y_CONVERSIONS
    assert c.flavour in {"e", "mu", "tau", "alpha"} and c.status in {"current", "future"}
    if c.label is not None:
        assert set(c.style) == {"color", "lw", "ls", "fill"}


@pytest.mark.parametrize("cid", CSV_IDS)
def test_csv_loads(cid):
    c = REGISTRY[cid]
    if cid in PENDING_HEADERS:
        warnings.warn(f"{cid}: CSV header pending, not checked")
        pytest.skip("header pending")
    assert c.description and c.references, "header needs a description and a source line"
    if cid in KNOWN_FLAVOUR_MISMATCH:
        pytest.xfail(KNOWN_FLAVOUR_MISMATCH[cid])
    if c.units.flavour is not None:
        assert c.units.flavour == c.flavour, "units line flavour does not match directory"
    mN, Usq = load_curve(c)
    assert len(mN) > 1
    assert np.all(np.isfinite(mN)) and np.all(mN > 0)
    assert np.all(np.isfinite(Usq)) and np.all(Usq >= 0)


@pytest.mark.parametrize("cid", [cid for cid, c in REGISTRY.items() if c.computed])
def test_computed_curve(cid):
    mN, Usq = load_curve(cid, mN=CURRENT_GRID[::50])
    assert np.all(np.isfinite(Usq)) and np.all(Usq > 0)
