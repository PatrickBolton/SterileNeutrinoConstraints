"""HNL physics: decay widths/lifetimes, 0vbb bounds and alpha_s."""

from .lifetime import N_DarkMatter_lifetime, N_lifetime
from .zeronubb import N_0vbb_bound

__all__ = ["N_lifetime", "N_DarkMatter_lifetime", "N_0vbb_bound"]
