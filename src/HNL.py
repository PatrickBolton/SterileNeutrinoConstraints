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
        'split': lambda x: 0.001**2 + x,
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
    

def plot(m_min=10**(-9.1), m_max=10**(4.1), U_min=10**(-12.1), U_max=10**(0.1),
         lim_type = 'current', seesaw = True, svfg = True,
         current_constraints = {}, future_constraints = {}, title = None, ylabel = None):

    fig, ax = plt.subplots()
    if seesaw == True:
        seesaw_bound = 0.05*1e-9/10**np.arange(-12, 6.2, 0.2)
        ax.plot(10**np.arange(-12, 6.2, 0.2), seesaw_bound, lw=1, ls='dotted', color='black') # Seesaw line
        plt.text(10**(1.4), 10**(-11.5), r'$\mathrm{Seesaw}$',fontsize=16,rotation=0,color='black')

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
    ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=14))
    ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10.0))
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=[i/10 for i in range(1,10)], numticks=15))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())

    # same for y-axis
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=14))
    ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10.0))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=[i/10 for i in range(1,10)], numticks=15))
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