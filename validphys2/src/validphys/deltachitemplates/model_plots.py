"""
Fit the data in the 'results' directory and plot the coefficients n dependent.
In 'results' are stored the files for every 'numreplicas' set passed to 'writechi2',
as (sigma fraction, 1 / sigma fraction, cv, std).

Name of the files: 'pdfname_hessian_neig_nrep_numreplicas.txt'

Use optional argument --show to plt.show() all plots. 
"""
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

from validphys import mplstyles
plt.style.use(str(mplstyles.smallstyle))
from validphys import plotutils

from lib import fit_delta_chi2, model_delta_chi2, coeff_n_dep, coeff_n_indep, model_std

# path to store the plots
path_save = '../../reports/results/'

coeff_names = {0: '$a_N$', 1: '$b_N$', 2: '$c_N$', 3: '$d_N$'}

def plot_chi2(axis, one_over_k, delta_chi2, std, coeff, nrep, title):

    xmin, xmax = np.min(one_over_k), np.max(one_over_k)
    x = np.linspace(xmin, xmax)
    y = model_delta_chi2(x, *coeff)
    # model results
    axis.plot(x, y, label='fit')
    # sort the results
    s = np.argsort(one_over_k)
    x = [one_over_k[i] for i in s]
    y = [delta_chi2[i] for i in s]
    err = [std[i] for i in s]
    # data results
    axis.errorbar(x, y, yerr=err, fmt='.', label=title)

    nbatch = np.floor(1000/nrep)
    axis.set_title('%d batches - %d replicas per batch' % (nbatch, nrep), fontsize=12)
    axis.grid(linewidth=.5)
    #axis.legend()

    return

def plot_coeff(which_method):
    
    cname = {'a': 0, 'b': 1, 'c': 3, 'd': 4}
    # Plot for single coefficients
    font = {'size': 26}

    import matplotlib
    matplotlib.rc('font', **font)
    coeff_name = input('Which coefficient:')
    i = cname[coeff_name]
    stdup = coeff_n[:,i] + coeff_std[:,i]
    stddown = coeff_n[:,i] - coeff_std[:,i]

    fig, ax = plt.subplots(1, figsize=[9.4, 7.8])
    # same colors as validphys plots
    color = next(ax._get_lines.prop_cycler)['color']

    ax.fill_between(fit_x, y1=stdup, y2=stddown, label='model',
            hatch=hatch, alpha=.5, edgecolor=color, zorder=1)
    ax.plot(fit_x, coeff_n[:,i], color=color)

    color = next(ax._get_lines.prop_cycler)['color']
    ax.errorbar(fit_x, fit_y[:,i], yerr=fit_err[:,i], fmt='o', elinewidth=2, capthick=2.5, capsize=3, label='fit', color=color)

    neig = which_method.split('_')[-1]
    if which_method.split('_')[0] == 'NNPDF31':
        ax.set_title('NNPDF3.1 - Coefficient $%s_N$ - %d eigenv' % (coeff_name, neig), fontsize=24)
    else:
        ax.set_title('n3fit - Coefficient $%s_N$ - %d eigenv' % (coeff_name, neig), fontsize=24)
    #ax.set_ylabel(coeff_names[str(i)])
    ax.set_xlabel('$N_{rep}$')
    #ax.set_xscale('log')
    ax.grid(linewidth=.5, which='both')
    ax.legend(prop={'size': 18})

    if args.save:
        plt.tight_layout()
        plt.savefig(path_save + args.pdf_results + '_coeff_' + coeff_name + '.pdf')
    
    return


parser = argparse.ArgumentParser()

parser.add_argument('pdf_results', help='name of the directory with the results to fit')
parser.add_argument('--plot_coeff', action='store_true', help='use to plot a single coefficient')
parser.add_argument('--show', action='store_true',
                        help='show delta chi2 and coefficients plots')
parser.add_argument('--save', action='store_true', help='use to save plots')

args = parser.parse_args()

# directory containing txt files to fit 
txt_dir = '../dat/' + args.pdf_results

# list all txt files with results
# elements are os.DirEntry objects 
files = [f for f in os.scandir(txt_dir) if f.is_file()]
# order files in ascending order
files = sorted(files, key=lambda e: int(e.name.split('_')[-1][:-4]))

# prepare array for plots of fit results
num_of_sets = len(files)

fit_x = np.zeros(num_of_sets, dtype=int) # number of replicas
# 4 coefficients nrep dependent for every 'numreplicas' set
fit_y = np.zeros((num_of_sets, 4))
# 4 uncertainties for every 'numreplicas' set
fit_err = np.zeros((num_of_sets, 4))
# 4 coefficients nrep independent for every 'numreplicas' set
param = np.zeros((num_of_sets, 4))
param_err = np.zeros((num_of_sets, 4))


# set delta chi2 subplots
n_rows = int(np.ceil(len(files)/2))
fig, ax = plt.subplots(n_rows, 2, figsize=[14.4, 9.8], sharex='col')
if args.pdf_results.split('_')[0] == 'NNPDF31':
    title = 'NNPDF3.1'
elif args.pdf_results.split('_')[0] == 'PN3':
    title = 'n3fit'
else:
    title = 'feature scaling'
fig.suptitle('$\Delta\chi^2$ fit - %s - %s eigenvectors' % (title, args.pdf_results.split('_')[-1]))

# set y axis label only for left hand plots
for axis in ax[:,0]:
    axis.set_ylabel('$<\Delta\chi^2>$')

# set x label only for low plots
ax[n_rows-1,0].set_xlabel('$1/k$')
ax[n_rows-1,1].set_xlabel('$1/k$')
    
    
for i, f in enumerate(files):

    one_over_sigfrac, cv, std = np.loadtxt(txt_dir + '/' + f.name, unpack=True,
                                            usecols=(1, 2, 3))

    # number of replicas used in the current set
    nrep = int(f.name.split('_')[-1][:-4])

    # results of the fit -> coefficients nrep dependent to plot
    fit_x[i] = nrep
    # get parameters nrep dependent: fit <delta_chi2> as a function of 1/k
    fit_y[i,:], fit_err[i,:] = fit_delta_chi2(one_over_sigfrac, cv, std)

    # get parameters nrep independent
    param[i,:] = coeff_n_indep(*fit_y[i,:], nrep)
    param_err[i,:] = coeff_n_indep(*fit_err[i,:], nrep)

    plot_chi2(ax.reshape(-1)[i], one_over_sigfrac, cv, std, fit_y[i,:], nrep, title)

    if f == files[-1]:
        handles, labels = axis.get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper right', fontsize='small')   


plt.tight_layout()
plt.subplots_adjust(top=0.9)
# save delta chi2 plot
if args.save:
    delta_chi2_fig = path_save + args.pdf_results + '_chi2.pdf'
    plt.savefig(delta_chi2_fig)

# array with <a>, <b>, <c>, <d>
param_name = ['a', 'b', 'c', 'd']
param_av = np.mean(param, axis=0)
param_std = 1/len(param) * np.sqrt(np.sum(param_err**2, axis=0))
#param_std = np.std(param, axis=0)
print('Results of the fit for the coefficients independent from n:')
for name, av, std in zip(param_name, param_av, param_std):
    print('%s: %.6f +- %.6f' % (name, av, std))

# plot coeffs as a function of nrep as predicted by the model
lowlim = np.min(fit_x)
suplim = np.max(fit_x)
#nrep = np.linspace(lowlim, suplim)
#nrep = np.linspace(50, 1000)

# compute values predicted by the model
#coeff_n = np.zeros((len(nrep), 4))
#coeff_std = np.zeros((len(nrep), 4))
coeff_n = np.zeros((len(fit_x), 4))
coeff_std = np.zeros((len(fit_x), 4))
#for i, n in enumerate(nrep):
for i, n in enumerate(fit_x):
    coeff_n[i,:] = coeff_n_dep(*param_av, n)
    coeff_std[i,:] = model_std(*param_av, n)

# same hatch patterns as validphys
hatchit = plotutils.hatch_iter()
hatch = next(hatchit)

# define chi2 to estimate goodness fit vs model
chi2 = np.zeros(4)

# prepare axes for coefficients plots
fig, ax = plt.subplots(2, 2, figsize=[10.4, 5.8])
fig.suptitle('Coefficients: $a_N,b_N,c_N,d_N$ - %s' % title)#, args.pdf_results.split('_')[-1]))

for i, axis in enumerate(np.ndarray.flatten(ax)):

    for j in range(len(fit_x)):
        chi2[i] += (fit_y[j,i] - coeff_n[j,i])**2 / fit_err[j,i]**2
        
    # same colors as validphys plots
    color = next(axis._get_lines.prop_cycler)['color']

    # define error bands
    stdup = coeff_n[:,i] + coeff_std[:,i]
    stddown = coeff_n[:,i] - coeff_std[:,i]
    #axis.fill_between(nrep, y1=stdup, y2=stddown, label='model',
    #                    hatch=hatch, alpha=.5, edgecolor=color, zorder=1)
    axis.fill_between(fit_x, y1=stdup, y2=stddown, label='model',
                        hatch=hatch, alpha=.5, edgecolor=color, zorder=1)
    axis.plot(fit_x, coeff_n[:,i], color=color)

    color = next(axis._get_lines.prop_cycler)['color']
    axis.errorbar(fit_x, fit_y[:,i], yerr=fit_err[:,i], fmt='.-', label='fit', color=color)

    #axis.set_title('$\chi^2=%.0f$' % chi2[i])
    axis.set_ylabel(coeff_names[i])
    axis.set_xlabel('$N_{rep}$')
    #axis.set_xlabel('Number of replicas')
    axis.set_xscale('log')
    axis.grid(linewidth=.2, which='both')
    if axis == np.ndarray.flatten(ax)[-1]:
        handles, labels = axis.get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper right', fontsize='small')   
    #axis.legend()
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    # save coefficients plot
    if args.save:
        coeffs_fig = path_save + args.pdf_results + '_coeffs.pdf'
        plt.savefig(coeffs_fig)


print(chi2)

if args.plot_coeff:
    plot_coeff(args.pdf_results)

if args.show:
    plt.show()
