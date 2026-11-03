import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import os
os.environ['PATH'] += os.pathsep + "/Library/TeX/texbin"
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

from .paths import PLOTS_DIR, REPO_ROOT
from .physics import N_0vbb_bound, N_DarkMatter_lifetime

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


def plot(x_min=1e-14, x_max=1e5, y_min=1e-15, y_max=1, n_x_ticks = 100, n_y_ticks = 100, plot_type = 'wide',
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

    m1NO = 0
    m2NO = np.sqrt(m1NO**2 + 7.537e-5) * 1e-9
    m3NO = np.sqrt(m2NO**2 + 2.521e-3) * 1e-9

    if seesaw == True:
        mN = np.logspace(-16,10,100)
        seesaw_bound_d = m2NO/(m2NO + mN)
        seesaw_bound_u = 0.1e-9/(0.1e-9 + mN)
        nlevels = 30
        cmap=cm.Greys
        cols = cmap(np.linspace(0.1,0.45,nlevels))
        levels = (np.linspace(1,np.sqrt(13),nlevels))**2
        for i in range(nlevels-1):
            ff(mN,np.maximum(seesaw_bound_u/levels[i],seesaw_bound_d),np.minimum(seesaw_bound_d*levels[i],seesaw_bound_u),alpha=0.025,color=cols[i,:],zorder=-1000,lw=0)
        
        if xy_flip == True:
            plt.text(10**(-13.6), 10**(4.), r'$\mathrm{Seesaw}$',fontsize=16,rotation=0,color='black')
        else:
            plt.text(10**(3.5), 10**(-13.2), r'$\mathrm{Seesaw}$',fontsize=16,rotation=0,color='black')

    if v_masses == True:
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

        df = pd.read_csv(REPO_ROOT / con["data"], sep=",", names = ["X", "Y"], comment = "#")
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

        if name.startswith('0vbb'):
            mvals = np.logspace(-15.1,5.1,1000)
            Uvals = [N_0vbb_bound(m, con["data"][0], con["data"][1]) for m in mvals]
            data = np.zeros((len(mvals),2))
            data[:,0] = mvals
            data[:,1] = Uvals

        else:
            if con["data"] == False:
                data = np.array([[0,0],[0,0]])
            else:
                df = pd.read_csv(REPO_ROOT / con["data"], sep=",", names = ["X", "Y"], comment = "#")
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
    
            df = pd.read_csv(REPO_ROOT / con["data"], sep=",", names = ["X", "Y"], comment = "#")
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

        if name.startswith('0vbb'):
            mvals = np.logspace(-15.1,10.1,1000)
            Uvals = [N_0vbb_bound(m, con["data"][0], con["data"][1]) for m in mvals]
            data = np.zeros((len(mvals),2))
            data[:,0] = mvals
            data[:,1] = Uvals

        else:
            if con["data"] == False:
                data = np.array([[0,0],[0,0]])
            
            else:
                df = pd.read_csv(REPO_ROOT / con["data"], sep=",", names = ["X", "Y"], comment = "#")
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

    if plot_type == 'wide':
        ratio = (1 + np.sqrt(5))/2
        ax.set_aspect(((x1-x0)/(y1-y0))/ratio)
        fig.set_size_inches(15,15)

    if plot_type == 'tall':
        ratio = 2/(1 + np.sqrt(5))
        ax.set_aspect(((x1-x0)/(y1-y0))/ratio)
        fig.set_size_inches(10,10)

    if svfg == True:
        if lim_type == 'current':
            plt.savefig(PLOTS_DIR / f"{title}_constraints.pdf", bbox_inches='tight')
            plt.savefig(PLOTS_DIR / f"{title}_constraints.png", bbox_inches='tight')
        elif lim_type == 'future':
            plt.savefig(PLOTS_DIR / f"{title}_future.pdf", bbox_inches='tight')
            plt.savefig(PLOTS_DIR / f"{title}_future.png", bbox_inches='tight')

    return fig, ax