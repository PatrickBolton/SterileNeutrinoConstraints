import numpy as np
import pytest

from sterile_constraints.units import UnitsError, parse_units


def convert(line, x, y, **kwargs):
    return parse_units(line, **kwargs).convert((np.array([x], dtype=float), np.array([y], dtype=float)))


@pytest.mark.parametrize("line, x, y, mN, Usq", [
    ("# log10(mN/GeV), log10(|UeN|^2)", 0.0, -3.0, 1.0, 1e-3),
    ("# mN/MeV, log10(|UmuN|^2)", 500.0, -6.0, 0.5, 1e-6),
    ("# mN/keV, |UeN|^2/10^-3", 10.0, 2.0, 1e-5, 2e-3),
    ("# log10(mN/eV), log10(|UeN|^2)", 3.0, -2.0, 1e-6, 1e-2),
    ("# log10(mN/GeV), log10(|UtauN|)", 1.0, -3.0, 10.0, 1e-6),
    ("# log10(mN^2/eV^2), |UeN|^2", 2.0, 0.1, 1e-8, 0.1),
    ("# log10(|UeN|^2), log10(Δm41^2/eV^2)", -2.0, 2.0, 1e-8, 1e-2),  # mass in second column
])
def test_conversions(line, x, y, mN, Usq):
    m, u = convert(line, x, y)
    assert m[0] == pytest.approx(mN, rel=1e-12)
    assert u[0] == pytest.approx(Usq, rel=1e-12)


def test_sin2_2theta_conversions():
    line = "# log10(mN/keV), log10(sin^2(2θ))"
    s = 0.36  # sin^2(2θ) = 4|U|^2(1 - |U|^2) for |U|^2 = 0.1
    _, exact = convert(line, 0.0, np.log10(s))
    _, small = convert(line, 0.0, np.log10(s), y_conversion="small_angle")
    _, none = convert(line, 0.0, np.log10(s), y_conversion="none")
    assert exact[0] == pytest.approx(0.1)
    assert small[0] == pytest.approx(s / 4)
    assert none[0] == pytest.approx(s)


def test_tan2_theta_exact():
    _, u = convert("# log10(tan^2(θ)), log10(Δm41^2/eV^2)", np.log10(0.25), 0.0)
    assert u[0] == pytest.approx(0.2)  # sin^2 = tan^2 / (1 + tan^2)


@pytest.mark.parametrize("line", [
    "# log10(mN/GeV)",
    "# log10(mN/TeV), log10(|UeN|^2)",
    "# log10(mN/GeV), log10(|UxN|^2)",
    "# log10(mN^2/GeV), log10(|UeN|^2)",
])
def test_bad_units_raise(line):
    with pytest.raises(UnitsError):
        parse_units(line)
