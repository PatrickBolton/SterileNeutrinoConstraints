import pytest
from matplotlib import pyplot as plt

from sterile_constraints.catalog import load_registry
from sterile_constraints.figures import FIGURES
from sterile_constraints.plotting import Curve, plot_figure


@pytest.mark.parametrize("name", sorted(FIGURES))
def test_figure_references_registered_constraints(name):
    spec = FIGURES[name]
    assert spec.name == name and spec.mode in {"current", "future"}
    registry = load_registry()
    for group in ("xray", "current", "cosmo", "future"):
        for item in getattr(spec, group):
            if isinstance(item, Curve):
                assert item.id in registry
                assert registry[item.id].label is not None, f"{item.id} has no default label/style"


@pytest.mark.parametrize("name", sorted(FIGURES))
def test_figure_builds(name):
    fig, ax = plot_figure(FIGURES[name])
    assert ax.lines or ax.collections
    plt.close(fig)
