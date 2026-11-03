"""Total HNL decay width and lifetime."""

import numpy as np
from scipy.integrate import fixed_quad
from scipy.special import spence

from .alphas import ALPHA_s_f
from .common import kin, lam

def R_QCD(mN):
    "Correction term for hadronic HNL decays."
    a = ALPHA_s_f(mN)
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
