"""Strong coupling alpha_s(mN) used in the hadronic HNL decay widths and the 0vbb QCD running."""

from importlib.resources import files

import numpy as np
from scipy.interpolate import interp1d

with files(__package__).joinpath("data", "alpha_data.txt").open() as _f:
    alpha_dat = np.loadtxt(_f)
ALPHA_s = interp1d(np.append(1e-7,alpha_dat[::,0]), np.append(alpha_dat[1,1],alpha_dat[::,1]), kind="linear")
def ALPHA_s_f(mN):
    # Interpolated table up to its last point (1 TeV), one-loop running with 6 flavours above
    mu = mN * np.heaviside(mN, 0) + 1e-7
    if mu <= ALPHA_s.x[-1]:
        a = ALPHA_s(mu)
    else:
        a = 2*np.pi/((11 - 2/3 * 6) * np.log(mN/0.043))
    return a

def ALPHA_s_ff(mN):
    mbottom = 4.1859
    mtop = 172.76
    if mN < mbottom:
        b0 = 11 - 2/3 * 4
        a = 2 * np.pi/(b0 * np.log(mN/0.119))
    elif mN >= mbottom and mN < mtop:
        b0 = 11 - 2/3 * 5
        a = 2 * np.pi/(b0 * np.log(mN/0.087))
    elif mN >= mtop:
        b0 = 11 - 2/3 * 6
        a = 2 * np.pi/(b0 * np.log(mN/0.043))
    return a
