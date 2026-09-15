import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os
os.environ['PATH'] += os.pathsep + "/Library/TeX/texbin"
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

import numpy as np

import numpy as np

import numpy as np

def transform(value, *operations):
    functions = {
        'log10': np.log10,
        'sqrt': np.sqrt,
        'square': np.square,
        'pow10': lambda x: 10 ** x,
        'complement': lambda x: 1 - x
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
    

def plot(m_min=10**(-12.1), m_max=10**(0.1), U_min=10**(-9.1), U_max=10**(4.1),
         lim_type = 'current', seesaw = True, labels_on = True, svfg = True,
         constraints = {}, title = None, ylabel = None):

    fig, ax = plt.subplots()

    for name, con in constraints.items():

        label = con["label"]
        label_loc = con["loc"]
        color = con['color']
        ls = con['ls']
        fill = con['fill']
        xr = con['xr']
        yr = con['yr']

        df = pd.read_csv(con["data"],header=None, sep=",", names = ["X", "Y"])
        data = np.zeros((len(df),2))
        for i in range(len(df.index)):
            data[i,0] = transform(df.iloc[i]['X'], xr)
            data[i,1] = transform(df.iloc[i]['Y'], yr)

        plt.loglog(data[:,0], data[:,1], lw=1.5, ls=ls, color=color)

        if lim_type == 'current':
            plt.fill_between(data[:,0], 10, data[:,1], facecolor=color, alpha=fill)
        if labels_on == True:
            plt.text(10**label_loc[0], 10**label_loc[1], label,fontsize=13,rotation=0,color=color)

    if seesaw == True:
        seesaw_bound = 0.05*1e-9/10**np.arange(-12, 6.2, 0.2)
        ax.plot(10**np.arange(-12, 6.2, 0.2), seesaw_bound, lw=1, ls='dotted', color='black') # Seesaw line
        plt.text(10**(1.4), 10**(-11.5), r'$\mathrm{Seesaw}$',fontsize=16,rotation=0,color='black')

    ax.set_xlabel(r'$m_N~[\mathrm{GeV}]$',fontsize=30,rotation=0)
    ax.set_ylabel(ylabel,fontsize=30,rotation=90)
    ax.tick_params(axis='x', which='major', labelsize=26, pad=7.5)
    ax.tick_params(axis='y', which='major', labelsize=26)
    ax.set_xlim(m_min, m_max)
    ax.set_ylim(U_min, U_max)
    x0,x1 = np.log10(ax.get_xlim())
    y0,y1 = np.log10(ax.get_ylim())
    ax.set_aspect(2*(x1-x0)/(1+np.sqrt(5))/(y1-y0))
    fig.set_size_inches(15,15)

    if svfg == True:
        plt.savefig(f"plots/{title}_constraints.pdf", bbox_inches='tight')
        plt.savefig(f"plots/{title}_constraints.png", bbox_inches='tight')

    return fig, ax