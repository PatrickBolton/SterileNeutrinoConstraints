import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
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

def N_DarkMatter_lifetime(mN):
    GeVtoinvsec = 1/6.5821e-25
    GF = 1.1663788e-5 # GeV^-2
    ALPHA = 1/137.035999084
    Gamma_3nu = GF**2 * mN**5/(96 * np.pi**3)
    Gamma_nugamma = 9 * ALPHA * GF**2 * mN**5/(256 * np.pi**4)
    return 1/((Gamma_3nu + Gamma_nugamma) * GeVtoinvsec) * 1/(3600*24*365*13.8e9)

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
    
def plot(m_min=10**(-14.1), m_max=10**(4.1), U_min=10**(-16.1), U_max=10**(0.1),
         lim_type = 'current', seesaw = True, v_masses = True, DM_window = True, svfg = True,
         current_constraints = {}, 
         current_xray_constraints = {}, 
         current_cosmo_constraints = {},
         future_constraints = {}, title = None, ylabel = None):

    fig, ax = plt.subplots()

    if seesaw == True:
        seesaw_bound = 0.05e-9/(0.05e-9 + 10**np.arange(-16, 6.2, 0.2))
        ax.plot(10**np.arange(-16, 6.2, 0.2), seesaw_bound, lw=1, ls='dotted', color='black') # Seesaw line
        plt.text(10**(3), 10**(-13.2), r'$\mathrm{Seesaw}$',fontsize=16,rotation=0,color='black')

    if v_masses == True:
        m1NO = 0
        m2NO = np.sqrt(m1NO**2 + 7.537e-5) * 1e-9
        m3NO = np.sqrt(m2NO**2 + 2.521e-3) * 1e-9
        plt.loglog([m2NO,m2NO], [1e-17,10], lw=1, ls='--', color='grey')
        plt.loglog([m3NO,m3NO], [1e-17,10], lw=1, ls='--', color='grey')
        plt.text(10**(-11.47), 10**(-15.7), r'$m_2^{\mathrm{min}} = \sqrt{\Delta m_\odot^2}$',fontsize=13,rotation=90,color='k')
        plt.text(10**(-10.7), 10**(-15.7), r'$m_3^{\mathrm{min}} = \sqrt{\Delta m_\odot^2 + \Delta m_{\mathrm{atm}}^2}$',fontsize=13,rotation=90,color='k')

    if DM_window == True:
        plt.loglog([1e-8,1e-8,7.4e-4,7.4e-4], [1e-17,1e-5,1e-5,1e-17], lw=0.5, ls='--', color='grey')
        plt.loglog(np.logspace(-4.695,-3.1307682803,100), 
                   N_DarkMatter_lifetime(np.logspace(-4.695,-3.1307682803,100)), lw=0.5, ls='--', color='blue')
        plt.text(10**(-4.2), 10**(-8.4), r'$\tau_N < t_U$',fontsize=13,rotation=-70,color='blue')
        plt.text(10**(-8.35), 10**(-15.7), r'$N$ Dark Matter',fontsize=13,rotation=90,color='k')

    for name, con in current_constraints.items():

        label = con["label"]
        label_on = con["label_on"]
        label_loc = con["label_loc"]
        label_rot = con["label_rot"]
        color = con['color']
        lw = con['lw']
        ls = con['ls']
        fill = con['fill']
        fill_band = con['fill_band']
        xr = con['xr']
        yr = con['yr']
        xy = con['xy']

        if name == '0vbb':
            mvals = np.logspace(-15.1,4.1,100)
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
            plt.loglog(data[:,0], data[:,1], lw=lw, ls=ls, color=color)
            if fill_band == 'up':
                plt.fill_between(data[:,0], data[:,1]*10, data[:,1], facecolor=color, alpha=fill)
            elif fill_band == 'down':
                plt.fill_between(data[:,0], data[:,1]/10, data[:,1], facecolor=color, alpha=fill)
            else:
                plt.fill_between(data[:,0], 10, data[:,1], facecolor=color, alpha=fill)
            if label_on == True:
                plt.text(10**label_loc[0], 10**label_loc[1], label,fontsize=13,rotation=label_rot,color=color)

        elif lim_type == 'future':
            plt.loglog(data[:,0], data[:,1], lw=lw, ls=ls, color=color, alpha = 0)
            if fill_band == 'up':
                plt.fill_between(data[:,0], data[:,1]*10, data[:,1], facecolor='lightgray', alpha=1)
            elif fill_band == 'down':
                plt.fill_between(data[:,0], data[:,1]/10, data[:,1], facecolor='lightgray', alpha=1)
            else:
                plt.fill_between(data[:,0], 10, data[:,1], facecolor='lightgray', alpha=1)

    for name, con in current_xray_constraints.items():
            label = con["label"]
            label_on = con["label_on"]
            label_loc = con["label_loc"]
            label_rot = con["label_rot"]
            color = con['color']
            lw = con['lw']
            ls = con['ls']
            fill = con['fill']
            fill_band = con['fill_band']
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

            if lim_type == 'current':
                plt.loglog(data[:,0], data[:,1], lw=lw, ls=ls, color=color)
                if fill_band == 'up':
                    up_lim = np.zeros(len(data[:,1])) + 1e-5
                    plt.fill_between(data[:,0], np.minimum(data[:,1]*500,up_lim), data[:,1], facecolor=color, alpha=fill)
                elif fill_band == 'down':
                    plt.fill_between(data[:,0], data[:,1]/10, data[:,1], facecolor=color, alpha=fill)
                else:
                    plt.fill_between(data[:,0], 10, data[:,1], facecolor=color, alpha=fill)
                if label_on == True:
                    plt.text(10**label_loc[0], 10**label_loc[1], label,fontsize=13,rotation=label_rot,color=color)
    
            if lim_type == 'future':
                if lim_type == 'future':
                    up_lim = np.zeros(len(data[:,1])) + 1e-5
                    plt.loglog(data[:,0], data[:,1], lw=lw, ls=ls, color=color, alpha = 0)
                    plt.fill_between(data[:,0], np.minimum(data[:,1]*1000,up_lim), data[:,1], facecolor='lightgray', alpha=1)

    for name, con in current_cosmo_constraints.items():
            label = con["label"]
            label_on = con["label_on"]
            label_loc = con["label_loc"]
            label_rot = con["label_rot"]
            color = con['color']
            lw = con['lw']
            ls = con['ls']
            fill = con['fill']
            fill_band = con['fill_band']
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

            plt.loglog(data[:,0], data[:,1], lw=lw, ls=ls, color=color)
            if fill_band == 'up':
                plt.fill_between(data[:,0], data[:,1]*10, data[:,1], facecolor=color, alpha=fill)
            elif fill_band == 'down':
                plt.fill_between(data[:,0], data[:,1]/10, data[:,1], facecolor=color, alpha=fill)
            else:
                plt.fill_between(data[:,0], 10, data[:,1], facecolor=color, alpha=fill)
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
        fill_band = con['fill_band']
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
            plt.loglog(data[:,0], data[:,1], lw=lw, ls=ls, color=color)
            if label_on == True:
                plt.text(10**label_loc[0], 10**label_loc[1], label,fontsize=13,rotation=label_rot,color=color)


    ax.set_xlabel(r'$m_N~[\mathrm{GeV}]$',fontsize=30,rotation=0)
    ax.set_ylabel(ylabel,fontsize=30,rotation=90)
    ax.tick_params(axis='x', which='major', labelsize=26, pad=7.5)
    ax.tick_params(axis='y', which='major', labelsize=26)
    ax.set_xlim(m_min, m_max)
    ax.set_ylim(U_min, U_max)
    # Major ticks: one per decade, labeled as 10^n
    ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=19))
    ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10.0))
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=[i/10 for i in range(1,10)], numticks=19))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())

    # same for y-axis
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=17))
    ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10.0))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=[i/10 for i in range(1,10)], numticks=17))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())

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