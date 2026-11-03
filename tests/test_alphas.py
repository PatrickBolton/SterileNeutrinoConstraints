import numpy as np
import pytest

from sterile_constraints.physics import N_lifetime
from sterile_constraints.physics.alphas import ALPHA_s, ALPHA_s_f


@pytest.mark.parametrize("mN", [0.0, 0.5, 4.2, 91.1876, 999.0, 1e3 - 1e-8, 1e3, 1e3 + 1e-8, 1e4, 1e5])
def test_alpha_s_defined_everywhere(mN):
    a = ALPHA_s_f(mN)
    assert np.isfinite(a) and 0.05 < a < 0.5


def test_table_used_up_to_its_end():
    assert ALPHA_s.x[-1] == 1e3
    assert ALPHA_s_f(999.0) == ALPHA_s(999.0 + 1e-7)


def test_lifetime_at_1_TeV():
    assert np.isfinite(N_lifetime(1e3, 1, 0, 0))
