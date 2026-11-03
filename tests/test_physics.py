import pytest

from sterile_constraints.physics import N_lifetime


@pytest.mark.parametrize("mN", [5.0, 10.0, 30.0, 70.0])
def test_electron_and_muon_widths_agree_above_tau_threshold(mN):
    """Lepton universality: e and mu mixing differ only by small charged-lepton mass effects."""
    assert N_lifetime(mN, 1, 0, 0) == pytest.approx(N_lifetime(mN, 0, 1, 0), rel=1e-2)


@pytest.mark.parametrize("mN", [30.0, 70.0, 200.0])
def test_tau_width_agrees_well_above_tau_threshold(mN):
    assert N_lifetime(mN, 0, 0, 1) == pytest.approx(N_lifetime(mN, 1, 0, 0), rel=3e-2)
