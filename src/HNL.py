import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import interp1d
from scipy.integrate import fixed_quad
from scipy.special import spence
import os
os.environ['PATH'] += os.pathsep + "/Library/TeX/texbin"
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def transform(value, *operations):
    functions = {
        'log10': np.log10,
        'sqrt': np.sqrt,
        'square': np.square,
        'pow10': lambda x: 10 ** x,
        'comp': lambda x: 1 - x,
        'split': lambda x: 0.0**2 + x,
    }

    # Flatten: allow operations as separate args or as a list/tuple
    flat_ops = []
    for op in operations:
        if isinstance(op, (list, tuple)):
            flat_ops.extend(op)
        else:
            flat_ops.append(op)

    result = value
    for operation in flat_ops:
        if isinstance(operation, (int, float)):
            result = result * operation
        elif isinstance(operation, str):
            if operation not in functions:
                raise ValueError(f"Unknown operation: '{operation}'. "
                                  f"Available: {list(functions.keys())}")
            result = functions[operation](result)
        else:
            raise TypeError(f"Each operation must be a number or string, got {type(operation)}")

    return result

def lam(x, y, z):
    "Square root of Källén function."
    kallen = (x - y - z)**2 - 4 * y * z
    return np.sqrt(kallen * np.heaviside(kallen,0))

def sqrt_pos(x):
    "Restrict square root to positive arguments." 
    return np.sqrt(x * np.heaviside(x, 0))

def kin(a, b, c):
    return np.heaviside(1 - b/a - c/a, 0)

alpha_dat = np.loadtxt('src/alpha_data.txt')
ALPHA_s = interp1d(np.append(1e-7,alpha_dat[::,0]), np.append(alpha_dat[1,1],alpha_dat[::,1]), kind="cubic")
def R_QCD(mN):
    "Correction term for hadronic HNL decays."
    a = ALPHA_s(mN * np.heaviside(mN, 0) + 1e-7)
    return 1 + a/np.pi + 5.2 * (a/np.pi)**2 + 26.4 * (a/np.pi)**3

def R(a, b, c):
    Rsq = 1 - (a + b + c)**2
    return np.sqrt(Rsq * np.heaviside(Rsq, 0))

def t1(ma, mb, st=10):
    arg = np.minimum(- st * (1 + ma - mb), 500)
    return 1/(1 + np.exp(arg))

def t2(ma, mb, st=10):
    return 1/(1 + np.exp(st * (1 + ma - mb)))

def I1(x, y, z):
    "I1 integral for decay rates."
    if x + y + z > 1:
        int = 0
    elif x == 0 and y == 0 and z == 0:
        int = 1
    elif y == 0 and z == 0:
        int = 1 - 8*x**2 + 8*x**6 - x**8 - 24*x**4*np.log(x)
    elif x == 0 and y == z:
        int = (np.sqrt(1 - 4*z**2)*(1 - 2*z**2*(7 + z**2 + 6*z**4)) 
               + 12*z**4*(-1 + z**4)*np.log((2*z**2*(1 - np.sqrt(1 - 4*z**2)))/((1 + np.sqrt(1 - 4*z**2))*(1 - 2*z**2 + np.sqrt(1 - 4*z**2)))))
    else:
        def intg(s):
            return (12 * 1/s * (1 + z**2 - s) * (s - x**2 - y**2) * lam(1, s, z**2) * lam(s, x**2, y**2) * np.heaviside(1 - x - y - z, 0))
        int, _ = fixed_quad(intg, (x + y)**2, (1 - z)**2, n=50)
    return int

def I2(x, y, z):
    "I2 integral for decay rates."
    if x + y + z > 1:
        int = 0
    elif x == 0 or y == 0:
        int = 0
    elif z == 0 and x == y:
        int = (8*x**2*np.sqrt(1 - 4*x**2)*(1 + 5*x**2 - 6*x**4) 
               - 48*x**4*(1 - 2*x**2 + 2*x**4)*np.log((1 - 2*x**2 + np.sqrt(1 - 4*x**2))/(2.*x**2)))
    else:
        def intg(s):
            return (24 * x * y * 1/s * (1 + z**2 - s) * lam(1, s, z**2) * lam(s, x**2, y**2) * np.heaviside(1 - x - y - z, 0))
        int, _ = fixed_quad(intg, (x + y)**2, (1 - z)**2, n=50)
    return int

def N_lifetime(mN, UeN, UmuN, UtauN):
    GF = 1.1663788e-5 # GeV^-2
    v = 1/np.sqrt(np.sqrt(2) * GF)
    ALPHA = 1/137.035999084
    SWsq = 0.23129
    MW = 80.379
    MZ = 91.1876
    MH = 125.25
    
    # Fermion masses in GeV
    Me = 0.510998946151099894e-3
    Mmu = 105.6583745e-3
    Mtau = 1776.86e-3
    Mu = 2.16e-3
    Md = 4.67e-3
    Ms = 93.4e-3
    Mc = 1.27
    Mb = 4.18
    Mt = 172.9

    # Meson Decay Factors
    Vud = 0.97370
    Vus = 0.2245
    Vub = 3.82e-3
    Vcd = 0.221
    Vcs = 0.987
    Vcb = 41.0e-3
    Mpi0 = 134.9768e-3
    Meta = 547.862e-3
    Metap = 957.78e-3
    Mrho0 = 775.26e-3
    Momega = 782.65e-3
    Mphi = 1019e-3
    Mpi = 139.57039e-3
    MK = 493.677e-3
    MD = 1.86962
    MDs = 1.96835
    MB = 5.27934
    Mrho = 892e-3
    MKs = 892e-3

    fpi0 = 130.3e-3
    feta = 78.4e-3
    fetap = -95.7e-3
    frho0 = 220e-3
    fomega = 195e-3
    fphi = 229e-3
    fpi = 130.3e-3
    fK = 156.4e-3
    fD = 222.6e-3
    fDs = 280.1e-3
    fB = 190e-3
    frho = 220e-3
    fKs = 204e-3
    krho0 = 1 - 2*SWsq
    komega = -2/3*SWsq
    kphi = -np.sqrt(2)*(1/2 - 2/3*SWsq)
    
    # Fermion Z couplings
    gLe = -1/2 + SWsq
    gRe = SWsq
    gLu = 1/2 - 2/3*SWsq
    gRu = - 2/3*SWsq
    gLd = -1/2 + 1/3*SWsq
    gRd = 1/3*SWsq

    w = t1(1.05*MW, mN, 1)
    z = t1(1.05*MZ, mN, 1)
    h = t1(MH, mN, 1)
    
    Uesq = np.abs(UeN)**2
    Umsq = np.abs(UmuN)**2
    Utsq = np.abs(UtauN)**2
    Usq = np.abs(UeN)**2 + np.abs(UmuN)**2 + np.abs(UtauN)**2

    two = GF**2 * mN**3/(16 * np.pi)
    thr = GF**2 * mN**5/(96 * np.pi**3)

    def f(r):
        if r <= 1e-3:
            return 3/4
        elif r > 1e-3 and r <= 1:
            return 1/r - (1-r)/r**2 * spence(1-r)
        else:
            return 1/r - (1-r)/r**2 * (np.pi**2/3 - 1/2 * np.log(r)**2 - spence(1-1/r) + 1j*np.pi*np.log(r))
    
    Ga_3v = z * thr * Usq
    Ga_vA = two * mN**2 * ALPHA/np.pi**3 * Usq * np.abs(f(mN**2/MW**2))**2
    
    Ga_vee = thr * Uesq * (((w + z*gLe)**2 + z*gRe**2) * I1(0, Me/mN, Me/mN) + (w + z*gLe) * z*gRe * I2(0, Me/mN, Me/mN))
    Ga_vmm = thr * Umsq * (((w + z*gLe)**2 + z*gRe**2) * I1(0, Mmu/mN, Mmu/mN) + (w + z*gLe) * z*gRe * I2(0, Mmu/mN, Mmu/mN))
    Ga_vtt = thr * Utsq * (((w + z*gLe)**2 + z*gRe**2) * I1(0, Mtau/mN, Mtau/mN) + (w + z*gLe) * z*gRe * I2(0, Mtau/mN, Mtau/mN))

    Ga_vem = w * thr/2 * (Uesq * I1(0, Me/mN, Mmu/mN) + Umsq * I1(0, Mmu/mN, Me/mN))
    Ga_vet = w * thr/2 * (Uesq * I1(0, Me/mN, Mtau/mN) + Utsq * I1(0, Mtau/mN, Me/mN))
    Ga_vmt = w * thr/2 * (Umsq * I1(0, Mmu/mN, Mtau/mN) + Utsq * I1(0, Mtau/mN, Mmu/mN))

    Ga_vpi0 = z * two * Usq * fpi0**2 * (1 - (Mpi0/mN)**2)**2 * kin(mN, 0, Mpi0)
    Ga_veta = z * two * Usq * feta**2 * (1 - (Meta/mN)**2)**2 * kin(mN, 0, Meta)
    Ga_vetap = z * two * Usq * fetap**2 * (1 - (Metap/mN)**2)**2 * kin(mN, 0, Metap)
    Ga_vrho0 = z * two * Usq * krho0**2 * frho0**2 * (1 - (Mrho0/mN)**2)**2 * (1 + 2*(Mrho0/mN)**2) * kin(mN, 0, Mrho0)
    Ga_vomega = z * two * Usq * komega**2 * fomega**2 * (1 - (Momega/mN)**2)**2 * (1 + 2*(Momega/mN)**2) * kin(mN, 0, Momega)
    Ga_vphi = z * two * Usq * kphi**2 * fphi**2 * (1 + 2*(Mphi/mN)**2) * (1 - (Mphi/mN)**2)**2 * kin(mN, 0, Mphi)
    Ga_vmeson = Ga_vpi0 + Ga_veta + Ga_vetap + Ga_vrho0 + Ga_vomega + Ga_vphi

    Ga_epi = w * two * Uesq * Vud**2 * fpi**2 * lam(1, (Me/mN)**2, (Mpi/mN)**2) * ((1 - (Me/mN)**2)**2 - (1 + (Me/mN)**2) * (Mpi/mN)**2) * kin(mN, Me, Mpi)
    Ga_eK = w * two * Uesq * Vus**2 * fK**2 * lam(1, (Me/mN)**2, (MK/mN)**2) * ((1 - (Me/mN)**2)**2 - (1 + (Me/mN)**2) * (MK/mN)**2) * kin(mN, Me, MK)
    Ga_eD = w * two * Uesq * Vcd**2 * fD**2 * lam(1, (Me/mN)**2, (MD/mN)**2) * ((1 - (Me/mN)**2)**2 - (1 + (Me/mN)**2) * (MD/mN)**2) * kin(mN, Me, MD)
    Ga_eDs = w * two * Uesq * Vcs**2 * fDs**2 * lam(1, (Me/mN)**2, (MDs/mN)**2) * ((1 - (Me/mN)**2)**2 - (1 + (Me/mN)**2) * (MDs/mN)**2) * kin(mN, Me, MDs)
    Ga_erho = w * two * Uesq * Vud**2 * frho**2 * lam(1, (Me/mN)**2, (Mrho/mN)**2) * ((1 - (Me/mN)**2)**2 + (1 + (Me/mN)**2) * (Mrho/mN)**2 - 2*(Mrho/mN)**4) * kin(mN, Me, Mrho)
    Ga_eKs = w * two * Uesq * Vus**2 * fKs**2 * lam(1, (Me/mN)**2, (MKs/mN)**2) * ((1 - (Me/mN)**2)**2 + (1 + (Me/mN)**2) * (MKs/mN)**2 - 2*(MKs/mN)**4) * kin(mN, Me, MKs)
    Ga_emeson = Ga_epi + Ga_eK + Ga_eD + Ga_eDs + Ga_erho + Ga_eKs

    Ga_mpi = w * two * Umsq * Vud**2 * fpi**2 * lam(1, (Mmu/mN)**2, (Mpi/mN)**2) * ((1 - (Mmu/mN)**2)**2 - (1 + (Mmu/mN)**2) * (Mpi/mN)**2) * kin(mN, Mmu, Mpi)
    Ga_mK = w * two * Umsq * Vus**2 * fK**2 * lam(1, (Mmu/mN)**2, (MK/mN)**2) * ((1 - (Mmu/mN)**2)**2 - (1 + (Mmu/mN)**2) * (MK/mN)**2) * kin(mN, Mmu, MK)
    Ga_mD = w * two * Umsq * Vcd**2 * fD**2 * lam(1, (Mmu/mN)**2, (MD/mN)**2) * ((1 - (Mmu/mN)**2)**2 - (1 + (Mmu/mN)**2) * (MD/mN)**2) * kin(mN, Mmu, MD)
    Ga_mDs = w * two * Umsq * Vcs**2 * fDs**2 * lam(1, (Mmu/mN)**2, (MDs/mN)**2) * ((1 - (Mmu/mN)**2)**2 - (1 + (Mmu/mN)**2) * (MDs/mN)**2) * kin(mN, Mmu, MDs)
    Ga_mrho = w * two * Umsq * Vud**2 * frho**2 * lam(1, (Mmu/mN)**2, (Mrho/mN)**2) * ((1 - (Mmu/mN)**2)**2 + (1 + (Mmu/mN)**2) * (Mrho/mN)**2 - 2*(Mrho/mN)**4) * kin(mN, Mmu, Mrho)
    Ga_mKs = w * two * Umsq * Vus**2 * fKs**2 * lam(1, (Mmu/mN)**2, (MKs/mN)**2) * ((1 - (Mmu/mN)**2)**2 + (1 + (Mmu/mN)**2) * (MKs/mN)**2 - 2*(MKs/mN)**4) * kin(mN, Mmu, MKs)
    Ga_mmeson = Ga_mpi + Ga_mK + Ga_mD + Ga_mDs + Ga_mrho + Ga_mKs

    Ga_tpi = w * two * Utsq * Vud**2 * fpi**2 * lam(1, (Mtau/mN)**2, (Mpi/mN)**2) * ((1 - (Mtau/mN)**2)**2 - (1 + (Mtau/mN)**2) * (Mpi/mN)**2) * kin(mN, Mtau, Mpi)
    Ga_tK = w * two * Utsq * Vus**2 * fK**2 * lam(1, (Mtau/mN)**2, (MK/mN)**2) * ((1 - (Mtau/mN)**2)**2 - (1 + (Mtau/mN)**2) * (MK/mN)**2) * kin(mN, Mtau, MK)
    Ga_tD = w * two * Utsq * Vcd**2 * fD**2 * lam(1, (Mtau/mN)**2, (MD/mN)**2) * ((1 - (Mtau/mN)**2)**2 - (1 + (Mtau/mN)**2) * (MD/mN)**2) * kin(mN, Mtau, MD)
    Ga_tDs = w * two * Utsq * Vcs**2 * fDs**2 * lam(1, (Mtau/mN)**2, (MDs/mN)**2) * ((1 - (Mtau/mN)**2)**2 - (1 + (Mtau/mN)**2) * (MDs/mN)**2) * kin(mN, Mtau, MDs)
    Ga_trho = w * two * Utsq * Vud**2 * frho**2 * lam(1, (Mtau/mN)**2, (Mrho/mN)**2) * ((1 - (Mtau/mN)**2)**2 + (1 + (Mtau/mN)**2) * (Mrho/mN)**2 - 2*(Mrho/mN)**4) * kin(mN, Mtau, Mrho)
    Ga_tKs = w * two * Utsq * Vus**2 * fKs**2 * lam(1, (Mtau/mN)**2, (MKs/mN)**2) * ((1 - (Mtau/mN)**2)**2 + (1 + (Mtau/mN)**2) * (MKs/mN)**2 - 2*(MKs/mN)**4) * kin(mN, Mtau, MKs)
    Ga_tmeson = Ga_tpi + Ga_tK + Ga_tD + Ga_tDs + Ga_trho + Ga_tKs

    Ga_vuu = z * 3 * thr * Usq * ((gLu**2 + gRu**2) * I1(0, Mu/mN, Mu/mN) + gLu * gRu * I2(0, Mu/mN, Mu/mN))
    Ga_vcc = z * 3 * thr * Usq * ((gLu**2 + gRu**2) * I1(0, Mc/mN, Mc/mN) + gLu * gRu * I2(0, Mc/mN, Mc/mN))
    Ga_vdd = z * 3 * thr * Usq * ((gLd**2 + gRd**2) * I1(0, Md/mN, Md/mN) + gLd * gRd * I2(0, Md/mN, Md/mN))
    Ga_vss = z * 3 * thr * Usq * ((gLd**2 + gRd**2) * I1(0, Ms/mN, Ms/mN) + gLd * gRd * I2(0, Ms/mN, Ms/mN))
    Ga_vbb = z * 3 * thr * Usq * ((gLd**2 + gRd**2) * I1(0, Mb/mN, Mb/mN) + gLd * gRd * I2(0, Mb/mN, Mb/mN))
    Ga_vhad = (R_QCD(mN) * t2(0, mN) 
               * (R(0, Mpi/mN, Mpi/mN) * Ga_vuu
                  + R(0, MD/mN, MD/mN) * Ga_vcc
                  + R(0, Mpi/mN, Mpi/mN) * Ga_vdd 
                  + R(0, MK/mN, MK/mN) * Ga_vss
                  + R(0, MB/mN, MB/mN) * Ga_vbb))

    Ga_eud = w * 3 * thr/2 * Uesq * Vud**2 * I1(Me/mN, Mu/mN, Md/mN)
    Ga_eus = w * 3 * thr/2 * Uesq * Vus**2 * I1(Me/mN, Mu/mN, Ms/mN)
    Ga_eub = w * 3 * thr/2 * Uesq * Vub**2 * I1(Me/mN, Mu/mN, Mb/mN)
    Ga_ecd = w * 3 * thr/2 * Uesq * Vcd**2 * I1(Me/mN, Mc/mN, Md/mN)
    Ga_ecs = w * 3 * thr/2 * Uesq * Vcs**2 * I1(Me/mN, Mc/mN, Ms/mN)
    Ga_ecb = w * 3 * thr/2 * Uesq * Vcb**2 * I1(Me/mN, Mc/mN, Mb/mN)
    Ga_ehad = (R_QCD(mN - Me) * t2(Me, mN)
               * (R(Me/mN, Mpi/mN, Mpi/mN) * Ga_eud 
                  + R(Me/mN, Mpi/mN, MK/mN) * Ga_eus
                  + R(Me/mN, Mpi/mN, MB/mN) * Ga_eub
                  + R(Me/mN, MD/mN, Mpi/mN) * Ga_ecd
                  + R(Me/mN, MD/mN, MK/mN) * Ga_ecs
                  + R(Me/mN, MD/mN, MB/mN) * Ga_ecb))

    Ga_mud = w * 3 * thr/2 * Umsq * Vud**2 * I1(Mmu/mN, Mu/mN, Md/mN)
    Ga_mus = w * 3 * thr/2 * Umsq * Vus**2 * I1(Mmu/mN, Mu/mN, Ms/mN)
    Ga_mub = w * 3 * thr/2 * Umsq * Vub**2 * I1(Mmu/mN, Mu/mN, Mb/mN)
    Ga_mcd = w * 3 * thr/2 * Umsq * Vcd**2 * I1(Mmu/mN, Mc/mN, Md/mN)
    Ga_mcs = w * 3 * thr/2 * Umsq * Vcs**2 * I1(Mmu/mN, Mc/mN, Ms/mN)
    Ga_mcb = w * 3 * thr/2 * Umsq * Vcb**2 * I1(Mmu/mN, Mc/mN, Mb/mN)
    Ga_mhad = (R_QCD(mN - Mmu) * t2(Mmu, mN)
               * (R(Mmu/mN, Mpi/mN, Mpi/mN) * Ga_mud 
                  + R(Mmu/mN, Mpi/mN, MK/mN) * Ga_mus
                  + R(Mmu/mN, Mpi/mN, MB/mN) * Ga_mub
                  + R(Mmu/mN, MD/mN, Mpi/mN) * Ga_mcd
                  + R(Mmu/mN, MD/mN, MK/mN) * Ga_mcs
                  + R(Mmu/mN, MD/mN, MB/mN) * Ga_mcb))

    Ga_tud = w * 3 * thr/2 * Uesq * Vud**2 * I1(Mtau/mN, Mu/mN, Md/mN)
    Ga_tus = w * 3 * thr/2 * Uesq * Vus**2 * I1(Mtau/mN, Mu/mN, Ms/mN)
    Ga_tub = w * 3 * thr/2 * Uesq * Vub**2 * I1(Mtau/mN, Mu/mN, Mb/mN)
    Ga_tcd = w * 3 * thr/2 * Uesq * Vcd**2 * I1(Mtau/mN, Mc/mN, Md/mN)
    Ga_tcs = w * 3 * thr/2 * Uesq * Vcs**2 * I1(Mtau/mN, Mc/mN, Ms/mN)
    Ga_tcb = w * 3 * thr/2 * Uesq * Vcb**2 * I1(Mtau/mN, Mc/mN, Mb/mN)
    Ga_thad = (R_QCD(mN - Mtau) * t2(Mtau, mN)
               * (R(Mtau/mN, Mpi/mN, Mpi/mN) * Ga_tud 
                  + R(Mtau/mN, Mpi/mN, MK/mN) * Ga_tus
                  + R(Mtau/mN, Mpi/mN, MB/mN) * Ga_tub
                  + R(Mtau/mN, MD/mN, Mpi/mN) * Ga_tcd
                  + R(Mtau/mN, MD/mN, MK/mN) * Ga_tcs
                  + R(Mtau/mN, MD/mN, MB/mN) * Ga_tcb))

    Ga_vZ = (1 - z) * 2 * two * Usq * v**2 * (1 - (MZ/mN)**2)**2 * (1 + 2 * (MZ/mN)**2) * kin(mN, 0, MZ)

    Ga_eW = (1 - w) * 2 * two * Uesq * v**2 * lam(1, (Me/mN)**2, (MW/mN)**2) * ((1 - (Me/mN)**2)**2 + (1 + (Me/mN)**2) * (MW/mN)**2 - 2 * (MW/mN)**4) * kin(mN, Me, MW)
    Ga_mW = (1 - w) * 2 * two * Umsq * v**2 * lam(1, (Mmu/mN)**2, (MW/mN)**2) * ((1 - (Mmu/mN)**2)**2 + (1 + (Mmu/mN)**2) * (MW/mN)**2 - 2 * (MW/mN)**4) * kin(mN, Mmu, MW)
    Ga_tW = (1 - w) * 2 * two * Utsq * v**2 * lam(1, (Mtau/mN)**2, (MW/mN)**2) * ((1 - (Mtau/mN)**2)**2 + (1 + (Mtau/mN)**2) * (MW/mN)**2 - 2 * (MW/mN)**4) * kin(mN, Mtau, MW)

    Ga_vh = (1 - h) * 2 * two * Usq * v**2 * (1 - (MH/mN)**2)**2 * kin(mN, 0, MH)
    
    Ga_tot = (Ga_vA
              + Ga_3v
              + Ga_vee + Ga_vmm + Ga_vtt + 2*(Ga_vem + Ga_vet + Ga_vmt)
              + t1(0, mN) * Ga_vmeson + Ga_vhad
              + t1(Me, mN) * 2*Ga_emeson + 2*Ga_ehad
              + t1(Mmu, mN) * 2*Ga_mmeson + 2*Ga_mhad
              + t1(Mtau, mN) * 2*Ga_tmeson + 2*Ga_thad
              + Ga_vZ + 2*(Ga_eW + Ga_mW + Ga_tW) + Ga_vh)

    return Ga_tot

def N_DarkMatter_lifetime(mN, UeN, UmuN, UtauN, time):
    GeVtoinvsec = 1/6.5821e-25
    Ga_tot = N_lifetime(mN, UeN, UmuN, UtauN)
    return 1/(Ga_tot * GeVtoinvsec) * 1/time

def N_0vbb_bound(mN, Thalf):
    fmtoinvGeV = 5.0677
    G01 = 1.4e-14 #yrs^-1
    Vud = 0.97
    Q = 2.5e-3 #GeV
    me = 0.511e-3 #GeV
    mpi = 0.1349768 #GeV
    gA = 1.27
    gvNN0 = - 1.01 * fmtoinvGeV**2
    g1pipi = 0.36
    g1piN = 1
    g1NN = (1 + 3 * gA**2)/4
    M0 = 2.7
    MFsd = -1.94
    MGTPPsd = 0.25
    MGTAPsd = -0.79
    MTPPsd = -0.01
    MTAPsd = 0.01
    RA = 1.2 * (136)**(1/3) * fmtoinvGeV

    def M(mN):
        ma = 0.157
        mb = 0.221
        return M0 * 1/(1 + mN/ma + (mN/mb)**2)

    def dMdmN(mN):
        ma = 0.157
        mb = 0.221
        return - M0 * (1/ma + 2*mN/mb**2)/(1 + mN/ma + (mN/mb)**2)**2

    def gvNN(mN):
        mc = 1.0
        md = 0.146
        return gvNN0 * (1 + (mN/mc)**2)/(1 + (mN/mc)**2 * (mN/md)**2)

    def Avfunc(mN):
        return - M(mN) - 2 * gvNN(mN) * mpi**2 * MFsd/gA**2

    def Avsoft(mN):
        def f(m, E):
            if m >= E:
                return (- 2 * (E * (1 + np.log(mpi/m)) 
                            + np.sqrt(m**2 - E**2) * (np.pi/2 - np.arctan2(E, np.sqrt(m**2 - E**2)))))
            else:
                return (- 2 * (E * (1 + np.log(mpi/m)) 
                            - np.sqrt(E**2 - m**2) * np.log((E + np.sqrt(E**2 - m**2))/m)))

        # E_n - E_i [MeV]
        energies = [0.17, 0.63, 0.89, 1.02, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5,
        1.6, 1.7, 1.9, 2.0, 2.2, 2.3, 2.6, 2.7, 2.9, 3.1,
        3.3, 3.6, 3.8, 4.0, 4.3, 4.6, 4.8, 5.1, 5.4, 5.7,
        6.1, 6.3, 6.7, 7.0, 7.3, 7.5, 7.7, 8.1, 8.4, 8.8,
        9.1, 9.4, 9.8, 10.1, 10.5, 10.9, 11.3, 11.7, 12.0, 12.4,
        12.8, 13.1, 13.5, 13.9, 14.2, 14.6, 15.1, 16.2, 17.3, 17.7]

        #<n| sigma tau+ |0_i^+>
        me_initial = [1.0, -0.19, -0.25, 0.30, 0.23, -0.13, 0.12, 0.16, -0.23, 0.20,
        -0.36, -0.24, 0.22, 0.34, 0.35, -0.49, 0.62, -0.91, 0.37, 0.30,
        0.39, 0.39, 0.45, -0.44, -0.35, -0.36, 0.44, 0.44, -0.55, 0.63,
        0.85, -1.2, -1.3, -1.9, 3.1, -4.0, 2.6, 1.4, -1.0, -0.93,
        0.80, 0.59, -0.50, 0.35, 0.26, -0.22, 0.17, -0.16, -0.16, 0.14,
        0.12, 0.092, -0.079, 0.071, -0.070, -0.035, -0.051, -0.039, -0.043, 0.11]

        #<0_f^+| sigma tau+ |n>
        me_final = [0.13, -0.0063, -0.016, 0.036, 0.025, -0.00076, -0.0052, -0.0028, -0.0098, -0.012,
        0.0084, 0.00058, 0.011, 0.0070, 0.0060, -0.0086, 0.021, -0.024, 0.0064, 0.0013,
        -0.0013, 0.0021, -0.013, -0.0032, -0.0038, -0.0067, 0.0083, 0.0066, -0.0093, 0.012,
        0.013, -0.016, -0.014, -0.016, 0.023, -0.028, 0.017, 0.0091, -0.0057, -0.0064,
        0.0038, 0.0014, 0.0027, -0.0027, -0.00053, -0.00021, -0.00037, -0.00054, -0.0010, 0.00092,
        -0.00014, -0.00040, -0.00019, -0.00026, 0.000031, 0.00021, -0.00015, 0.00011, -0.000091, -0.000029]

        list = []
        for i in range(len(energies)):
            DelE = Q/2 + me + energies[i] * 1e-3
            list.append(me_initial[i] * me_final[i] * 2 * f(mN, DelE))

        return (- RA/(2 * np.pi) * sum(list))

    if mN < 0.1:
        Av = Avfunc(mN) + mN * dMdmN(mN) + Avsoft(mN)
    elif mN >= 0.1 and mN < 2:
        Av = Avfunc(mN)
    elif mN >= 2:
        Av = - 2 * (mpi/mN)**2 * (5/6 * g1pipi * (MGTPPsd + MTPPsd)
                                  + 1/2 * g1piN * (MGTAPsd + MTAPsd) 
                                  - 2 * g1NN/gA**2 * MFsd)

    return (me/(Vud**2 * gA**2 * mN * np.abs(Av - 0*Avfunc(0)) * np.sqrt(G01 * Thalf)))
    
def plot(x_min=1e-14, x_max=1e5, y_min=1e-15, y_max=1, n_x_ticks = 100, n_y_ticks = 100,
         lim_type = 'current', seesaw = True, v_masses = True, DM_window = True, DM_lifetime = False, svfg = True,
         current_constraints = {}, 
         current_xray_constraints = {}, 
         current_cosmo_constraints = {},
         block_color = 'lightgrey',
         future_constraints = {}, title = None, xlabel = r'$m_N~[\mathrm{GeV}]$', ylabel = None, x_scale = 1, y_scale = 1, x_rot = 0, xy_flip = False):

    fig, ax = plt.subplots()

    if xy_flip == True:
        def ll(x, y, ax=None, **kwargs):
            if ax is None:
                ax = plt.gca()
            return ax.loglog(y, x, **kwargs)
        def ff(y, x1, x2=0, ax=None, **kwargs):
            if ax is None:
                ax = plt.gca()
            return ax.fill_betweenx(y, x1, x2, **kwargs) 
    else:
        def ll(x, y, ax=None, **kwargs):
            if ax is None:
                ax = plt.gca()
            return ax.loglog(x, y, **kwargs)
        def ff(x, y1, y2=0, ax=None, **kwargs):
            if ax is None:
                ax = plt.gca()
                return ax.fill_between(x, y1, y2, **kwargs)

    if seesaw == True:
        seesaw_bound = 0.05e-9/(0.05e-9 + 10**np.arange(-16, 10, 0.2))
        ll(10**np.arange(-16, 10, 0.2), seesaw_bound, lw=1, ls='dotted', color='black') # Seesaw line
        if xy_flip == True:
            plt.text(10**(-14.), 10**(4.), r'$\mathrm{Seesaw}$',fontsize=16,rotation=0,color='black')
        else:
            plt.text(10**(3.5), 10**(-13.5), r'$\mathrm{Seesaw}$',fontsize=16,rotation=0,color='black')

    if v_masses == True:
        m1NO = 0
        m2NO = np.sqrt(m1NO**2 + 7.537e-5) * 1e-9
        m3NO = np.sqrt(m2NO**2 + 2.521e-3) * 1e-9
        ll([m2NO,m2NO], [1e-40,10], lw=1, ls='--', color='grey')
        ll([m3NO,m3NO], [1e-40,10], lw=1, ls='--', color='grey')
        if xy_flip == True:
            plt.text(10**(-27.6), 10**(-10.85), r'$m_2^{\mathrm{min}} = \sqrt{\Delta m_\odot^2}$',fontsize=13,rotation=0,color='k')
            plt.text(10**(-27.6), 10**(-10.), r'$m_3^{\mathrm{min}} = \sqrt{\Delta m_\odot^2 + \Delta m_{\mathrm{atm}}^2}$',fontsize=13,rotation=0,color='k')
        else:
            plt.text(10**(-11.5), 10**(-14.6), r'$m_2^{\mathrm{min}} = \sqrt{\Delta m_\odot^2}$',fontsize=13,rotation=90,color='k')
            plt.text(10**(-10.73), 10**(-14.6), r'$m_3^{\mathrm{min}} = \sqrt{\Delta m_\odot^2 + \Delta m_{\mathrm{atm}}^2}$',fontsize=13,rotation=90,color='k')

    if DM_window == True:
        ll([1e-8,1e-8,2e-3,2e-3], [1e-17,1e-5,1e-5,1e-17], lw=0.5, ls='--', color='grey')
        mvals = np.logspace(-4.695,np.log10(2e-3),100)
        lvals = [N_DarkMatter_lifetime(m, 1, 0, 0, 3600*24*365*13.8e9) for m in mvals]
        ll(mvals, lvals, lw=1.0, ls='--', color='k')
        plt.text(10**(-4.2), 10**(-8.4), r'$\tau_N < t_U$',fontsize=13,rotation=-70,color='k')
        plt.text(10**(-8.35), 10**(-14.6), r'$N$ Dark Matter',fontsize=13,rotation=90,color='k')

    if DM_lifetime == True:
        mvals = np.logspace(-10,5,2000)
        lval_u = [N_DarkMatter_lifetime(m, 1, 0, 0, 60*60*24*365*13.8e9) for m in mvals]
        lvals_r = [N_DarkMatter_lifetime(m, 1, 0, 0, 60*60*24*365*3.78e5) for m in mvals]
        lvals_1sec = [N_DarkMatter_lifetime(m, 1, 0, 0, 1) for m in mvals]

        ll(mvals, lval_u, lw=1.0, ls='--', color='k')
        ll(mvals, lvals_r, lw=1.0, ls='--', color='k')
        ll(mvals, lvals_1sec, lw=1.0, ls='--', color='k')

    for name, con in current_xray_constraints.items():
        label = con["label"]
        label_on = con["label_on"]
        label_loc = con["label_loc"]
        label_rot = con["label_rot"]
        color = con['color']
        lw = con['lw']
        ls = con['ls']
        fill = con['fill']
        xr = con['xr']
        yr = con['yr']
        xy = con['xy']

        df = pd.read_csv(con["data"], sep=",", names = ["X", "Y"], comment = "#")
        data = np.zeros((len(df),2))
        for i in range(len(df.index)):
            if xy == True:
                data[i,0] = transform(df.iloc[i]['X'], xr)
                data[i,1] = transform(df.iloc[i]['Y'], yr)
            else:
                data[i,0] = transform(df.iloc[i]['Y'], xr)
                data[i,1] = transform(df.iloc[i]['X'], yr)

        up_lim = [N_DarkMatter_lifetime(m, 1, 0, 0, 60*60*24*365*13.8e9) for m in data[:,0]]
            
        if lim_type == 'current':
            ll(data[:,0], data[:,1], lw=lw, ls=ls, color=color)
            ff(data[:,0], np.minimum(up_lim, 1e-5), data[:,1], facecolor=color, alpha=fill)
            if label_on == True:
                plt.text(10**label_loc[0], 10**label_loc[1], label,fontsize=13,rotation=label_rot,color=color)
    
        if lim_type == 'future':
            if lim_type == 'future':
                ll(data[:,0], data[:,1], lw=lw, ls=ls, color=color, alpha = 0)
                ff(data[:,0], np.minimum(up_lim, 1e-5), data[:,1], facecolor=block_color, alpha=1)

    for name, con in current_constraints.items():

        label = con["label"]
        label_on = con["label_on"]
        label_loc = con["label_loc"]
        label_rot = con["label_rot"]
        color = con['color']
        lw = con['lw']
        ls = con['ls']
        fill = con['fill']
        xr = con['xr']
        yr = con['yr']
        xy = con['xy']

        if name == '0vbb':
            mvals = np.logspace(-15.1,10.1,1000)
            Uvals = [N_0vbb_bound(m, 3.8e26) for m in mvals]
            data = np.zeros((len(mvals),2))
            data[:,0] = mvals
            data[:,1] = Uvals

        else:
            if con["data"] == False:
                data = np.array([[0,0],[0,0]])
            else:
                df = pd.read_csv(con["data"], sep=",", names = ["X", "Y"], comment = "#")
                data = np.zeros((len(df),2))
                for i in range(len(df.index)):
                    if xy == True:
                        data[i,0] = transform(df.iloc[i]['X'], xr)
                        data[i,1] = transform(df.iloc[i]['Y'], yr)
                    else:
                        data[i,0] = transform(df.iloc[i]['Y'], xr)
                        data[i,1] = transform(df.iloc[i]['X'], yr)

        if lim_type == 'current':
            ll(data[:,0], data[:,1], lw=lw, ls=ls, color=color)
            ff(data[:,0], 10, data[:,1], facecolor=color, alpha=fill)
            if label_on == True:
                plt.text(10**label_loc[0], 10**label_loc[1], label,fontsize=13,rotation=label_rot,color=color)

        elif lim_type == 'future':
            ll(data[:,0], data[:,1], lw=lw, ls=ls, color=color, alpha = 0)
            ff(data[:,0], 10, data[:,1], facecolor=block_color, alpha=1)

    for name, con in current_cosmo_constraints.items():
            label = con["label"]
            label_on = con["label_on"]
            label_loc = con["label_loc"]
            label_rot = con["label_rot"]
            color = con['color']
            lw = con['lw']
            ls = con['ls']
            fill = con['fill']
            xr = con['xr']
            yr = con['yr']
            xy = con['xy']
    
            df = pd.read_csv(con["data"], sep=",", names = ["X", "Y"], comment = "#")
            data = np.zeros((len(df),2))
            for i in range(len(df.index)):
                if xy == True:
                    data[i,0] = transform(df.iloc[i]['X'], xr)
                    data[i,1] = transform(df.iloc[i]['Y'], yr)
                else:
                    data[i,0] = transform(df.iloc[i]['Y'], xr)
                    data[i,1] = transform(df.iloc[i]['X'], yr)

            ll(data[:,0], data[:,1], lw=lw, ls=ls, color=color)
            ff(data[:,0], 10, data[:,1], facecolor=color, alpha=fill)
            if label_on == True:
                plt.text(10**label_loc[0], 10**label_loc[1], label,fontsize=13,rotation=label_rot,color=color)

    for name, con in future_constraints.items():
        label = con["label"]
        label_on = con["label_on"]
        label_loc = con["label_loc"]
        label_rot = con["label_rot"]
        color = con['color']
        lw = con['lw']
        ls = con['ls']
        fill = con['fill']
        xr = con['xr']
        yr = con['yr']
        xy = con['xy']
            
        df = pd.read_csv(con["data"], sep=",", names = ["X", "Y"], comment = "#")
        data = np.zeros((len(df),2))
        for i in range(len(df.index)):
            if xy == True:
                data[i,0] = transform(df.iloc[i]['X'], xr)
                data[i,1] = transform(df.iloc[i]['Y'], yr)
            else:
                data[i,0] = transform(df.iloc[i]['Y'], xr)
                data[i,1] = transform(df.iloc[i]['X'], yr)

        if lim_type == 'future':
            ll(data[:,0], data[:,1], lw=lw, ls=ls, color=color)
            if label_on == True:
                plt.text(10**label_loc[0], 10**label_loc[1], label,fontsize=13,rotation=label_rot,color=color)

    ax.tick_params(axis='x', which='major', labelsize=26, pad=7.5, labelrotation=x_rot)
    ax.tick_params(axis='y', which='major', labelsize=26)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    # Major ticks: one per decade, labeled as 10^n
    class ScaledLogFormatterMathtext(ticker.LogFormatterMathtext):
        def __init__(self, scale=1.0, **kwargs):
            super().__init__(**kwargs)
            self.scale = scale
            
        def __call__(self, x, pos=None):
            return super().__call__(x * self.scale, pos)
    ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=n_x_ticks))
    #ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{x * 1e6:g}"))
    #ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10.0))
    ax.xaxis.set_major_formatter(ScaledLogFormatterMathtext(scale=x_scale, base=10.0))
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=[i/10 for i in range(1,10)], numticks=100))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xlabel(xlabel,fontsize=30,rotation=0)

    # same for y-axis
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=n_y_ticks))
    #ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10.0))
    ax.yaxis.set_major_formatter(ScaledLogFormatterMathtext(scale=y_scale, base=10.0))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=[i/10 for i in range(1,10)], numticks=100))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_ylabel(ylabel,fontsize=30,rotation=90)

    # Ticks on all four sides, mirrored, pointing inward
    ax.tick_params(axis='both', which='major', top=True, right=True,
                direction='in', length=8)
    ax.tick_params(axis='both', which='minor', top=True, right=True,
                direction='in', length=4)

    x0,x1 = np.log10(ax.get_xlim())
    y0,y1 = np.log10(ax.get_ylim())
    ax.set_aspect(2*(x1-x0)/(1+np.sqrt(5))/(y1-y0))
    fig.set_size_inches(15,15)

    if svfg == True:
        if lim_type == 'current':
            plt.savefig(f"plots/{title}_constraints.pdf", bbox_inches='tight')
            plt.savefig(f"plots/{title}_constraints.png", bbox_inches='tight')
        elif lim_type == 'future':
            plt.savefig(f"plots/{title}_future.pdf", bbox_inches='tight')
            plt.savefig(f"plots/{title}_future.png", bbox_inches='tight')

    return fig, ax