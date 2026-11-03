"""Kinematic helper functions shared by the decay-rate calculations."""

import numpy as np

def lam(x, y, z):
    "Square root of Källén function."
    kallen = (x - y - z)**2 - 4 * y * z
    return np.sqrt(kallen * np.heaviside(kallen,0))

def sqrt_pos(x):
    "Restrict square root to positive arguments." 
    return np.sqrt(x * np.heaviside(x, 0))

def kin(a, b, c):
    return np.heaviside(1 - b/a - c/a, 0)
