"""
Plot the coefficients a_N, b_N, c_N, d_N for all the subsets of chi2
considered
"""
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

from validphys import mplstyles
plt.style.use(str(mplstyles.smallstyle))
from validphys import plotutils

coeff_names = {0: '$a_N$', 1: '$b_N$', 2: '$c_N$', 3: '$d_N$'}

def name_label(ns):
    if 's' in ns and '1' in ns:
        return '20 eigenv', '.-'
    elif 's' in ns and '2' in ns:
        return '21-40 eigenv', '--'
    elif 's' in ns and '3' in ns:
        return '41-60 eigenv', '--'
    elif 's' in ns and '4' in ns:
        return '61-80 eigenv', '--'
    elif 's' in ns and '5' in ns:
        return '81-100 eigenv', '--'
    elif 'c' in ns and '1' in ns:
        return '40 eigenv', '.-'
    elif 'c' in ns and '2' in ns:
        return '60 eigenv', '.-'
    elif 'c' in ns and '3' in ns:
        return '80 eigenv', '.-'
    else:
        return '100 eigen', '.-'
    
def get_out_dir(fn):
    if fn == 1:
        out_dir = '../subdata_NNPDF/coeffs/'
        return out_dir, os.listdir(out_dir), 'NNPDF3.1'
    elif fn == 2:
        out_dir = '../subdata_n3fit/coeffs/'
        return out_dir, os.listdir(out_dir), 'n3fit'
    elif fn == 3:
        out_dir = '../subdata_feature_scaling/coeffs/'
        return out_dir, os.listdir(out_dir), 'feature scaling'
    else:
        out_dir = '../subdata_new_featsc/coeffs/'
        return out_dir, os.listdir(out_dir), 'feature scaling'

parser = argparse.ArgumentParser()

parser.add_argument('fitname', help='name of the fit. 1: NNPDF, 2: n3fit,'
                                    ' 3: feature scaling, 4: new feat. sc.', type=int,
                                    choices={1,2,3,4}) 
parser.add_argument('--save', help='save plots', action='store_true')

args = parser.parse_args()

out_dir, outfiles, title = get_out_dir(args.fitname)

# prepare axes for coefficients plots
fig, ax = plt.subplots(2, 2, figsize=[11.4, 7.4])
#fig.suptitle('$a_N,b_N,c_N,d_N$ with %s eigenvectors from %s MC sets of $N_{rep}$ replicas' % (outfiles[-1].split('_')[-2], title))
fig.suptitle('Coefficients: $a_N,b_N,c_N,d_N$ - %s' % title)

# elements are os.DirEntry objects 
files = [f for f in os.scandir(out_dir) if f.is_file()]
# order files in ascending order
sorted_outfiles = sorted(files, key=lambda e: e.name.split('_')[-1][:-4])
sorted_outfiles[0], sorted_outfiles[4] = sorted_outfiles[4], sorted_outfiles[0]

for data in sorted_outfiles:
    
    nrep, a_N, b_N, c_N, d_N, err_a, err_b, err_c, err_d = np.loadtxt(out_dir + data.name, unpack=True) 

    fit_coeff = np.asarray([a_N, b_N, c_N, d_N])
    fit_err = np.asarray([err_a, err_b, err_c, err_c])

    name_subset = data.name.split('_')[-1] 
    label, fmt = name_label(name_subset)

    for i, axis in enumerate(np.ndarray.flatten(ax)):
            
        # same colors as validphys plots
        color = next(axis._get_lines.prop_cycler)['color']

        axis.errorbar(nrep, fit_coeff[i], yerr=fit_err[i], fmt=fmt, label=label, color=color)

        #axis.set_title('$\chi^2=%.0f$' % chi2[i])
        axis.set_ylabel(coeff_names[i])
        axis.set_xlabel('$N_{rep}$')
        #axis.set_xlabel('Number of replicas')
        axis.set_xscale('log')
        axis.grid(linewidth=.2, which='both')
        #if i == 1:
        #    axis.set_title('Coefficient $b_N$ - %s' % title)

        #if data == sorted_outfiles[-1]:
        #    handles, labels = axis.get_legend_handles_labels()
        #    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=[.98,.93])   

plt.tight_layout()
plt.subplots_adjust(top=0.9)

# path to store the plots
path_save = '../../reports/results/'
if args.fitname == 1:
    name = 'NNPDF31'
elif args.fitname == 2:
    name = 'n3fit'
elif args.fitname == 3:
    name = 'feature_scaling'
else:
    name = 'newfeatsc'

if args.save:
    coeffs_fig = path_save + name + '_coeffs_subsets.pdf'
    plt.savefig(coeffs_fig)
    
plt.show()
