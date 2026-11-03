"""Parse the units line of a constraint CSV and convert its columns to (mN [GeV], |U|^2).

Every CSV header ends with a units line naming its two columns, e.g.

    # log10(mN/GeV), log10(|UeN|^2)
    # log10(sin^2(2θ)), log10(Δm41^2/eV^2)
    # mN/keV, |UeN|^2/10^-3

One column is a mass (mN, mN^2 or Δm41^2 in eV/keV/MeV/GeV) and the other a mixing quantity
(|U|^2, |U|, sin^2(2θ), tan^2(θ) or Omega_DM sin^2(2θ)), either of them optionally log10'd.
The conversion is expressed as a chain of elementary operations (see `transform`), so that the
result is bit-for-bit what the original per-notebook transforms produced.
"""

import re
from dataclasses import dataclass

import numpy as np

# Lightest active neutrino mass squared [eV^2] assumed when converting Δm41^2 to mN^2
M1SQ = 0.0**2
# 1/Omega_DM, to undo the Omega_DM weighting of Omega_DM sin^2(2θ) bounds
OMEGA_DM_INV = 4.083333333061112

MASS_UNITS = {"eV": 1e-9, "keV": 1e-6, "MeV": 1e-3, "GeV": 1}

# How a sin^2(2θ) or tan^2(θ) column is converted to |U|^2
Y_CONVERSIONS = ("exact", "small_angle", "none")

OPERATIONS = {
    "pow10": lambda x: 10**x,
    "sqrt": np.sqrt,
    "square": np.square,
    "comp": lambda x: 1 - x,
    "add_m1sq": lambda x: M1SQ + x,
    "tan2_to_sin2": lambda x: x / (1 + x),
}


def transform(value, operations):
    """Apply a chain of operations: names from OPERATIONS, or numbers to multiply by."""
    result = value
    for op in operations:
        if isinstance(op, (int, float)):
            result = result * op
        elif op in OPERATIONS:
            result = OPERATIONS[op](result)
        else:
            raise ValueError(f"Unknown operation {op!r}; available: {list(OPERATIONS)}")
    return result


class UnitsError(ValueError):
    pass


@dataclass(frozen=True)
class Units:
    """Parsed units line: which column holds the mass and how to convert each column."""

    line: str
    mass_column: int  # 0 or 1
    mass_ops: tuple
    mixing_ops: tuple
    flavour: str | None  # 'e', 'mu', 'tau', 'alpha', or None for sin^2(2θ) etc.

    def convert(self, columns):
        """Convert the two raw columns to (mN [GeV], |U|^2) arrays."""
        mass = columns[self.mass_column]
        mixing = columns[1 - self.mass_column]
        return transform(mass, self.mass_ops), transform(mixing, self.mixing_ops)


def _strip_log10(term):
    m = re.fullmatch(r"log10\((.*)\)", term)
    return (m.group(1), True) if m else (term, False)


_MASS_RE = re.compile(r"(?P<q>mN|mN\^2|Δm41\^2)/(?P<unit>eV|keV|MeV|GeV)(?P<sq>\^2)?")
_U_RE = re.compile(r"\|U(?P<flav>e|mu|tau|alpha)N\|(?P<sq>\^2)?(?P<scale>/10\^-3)?")


def _mass_ops(term):
    inner, log = _strip_log10(term)
    m = _MASS_RE.fullmatch(inner)
    if not m:
        return None
    squared = m["q"] != "mN"
    if squared != bool(m["sq"]):
        raise UnitsError(f"Inconsistent powers in mass term {term!r}")
    ops = ["pow10"] if log else []
    if m["q"] == "Δm41^2":
        ops.append("add_m1sq")
    if squared:
        ops.append("sqrt")
    if MASS_UNITS[m["unit"]] != 1:
        ops.append(MASS_UNITS[m["unit"]])
    return ops


def _mixing_ops(term, y_conversion):
    inner, log = _strip_log10(term)
    ops = ["pow10"] if log else []
    m = _U_RE.fullmatch(inner)
    if m:
        if m["scale"]:
            ops.append(1e-3)
        if not m["sq"]:
            ops.append("square")
        return ops, m["flav"]
    if inner.startswith("Omega_DM "):
        ops.append(OMEGA_DM_INV)
        inner = inner.removeprefix("Omega_DM ")
    if inner == "sin^2(2θ)":
        ops += {"exact": ["comp", "sqrt", "comp", 0.5], "small_angle": [0.25], "none": []}[y_conversion]
        return ops, None
    if inner == "tan^2(θ)":
        ops += {"exact": ["tan2_to_sin2"], "small_angle": [], "none": []}[y_conversion]
        return ops, None
    return None, None


def parse_units(line, y_conversion="exact"):
    """Parse a units line such as '# log10(mN/GeV), log10(|UeN|^2)'.

    y_conversion selects how sin^2(2θ) / tan^2(θ) become |U|^2: 'exact', 'small_angle'
    (sin^2(2θ)/4, tan^2(θ)) or 'none' (taken as |U|^2 directly).
    """
    if y_conversion not in Y_CONVERSIONS:
        raise UnitsError(f"y_conversion must be one of {Y_CONVERSIONS}, got {y_conversion!r}")
    terms = [t.strip() for t in line.lstrip("#").split(",")]
    if len(terms) != 2:
        raise UnitsError(f"Expected two comma-separated columns in units line {line!r}")
    for mass_column in (0, 1):
        mass_ops = _mass_ops(terms[mass_column])
        if mass_ops is None:
            continue
        mixing_ops, flavour = _mixing_ops(terms[1 - mass_column], y_conversion)
        if mixing_ops is None:
            raise UnitsError(f"Unrecognised mixing term {terms[1 - mass_column]!r} in {line!r}")
        return Units(line.strip(), mass_column, tuple(mass_ops), tuple(mixing_ops), flavour)
    raise UnitsError(f"No recognised mass column in units line {line!r}")
