""" Plot the average deltachi2 as a function of x=1/k for all
the Nrep used
"""
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

from validphys import mplstyles

plt.style.use(str(mplstyles.smallstyle))
from validphys import plotutils

from lib import fit_delta_chi2, model_delta_chi2, coeff_n_indep

# custom legend
from matplotlib.lines import Line2D


def plot_chi2(axis, one_over_k, delta_chi2, std, coeff, nrep, txt_dir):

    xmin, xmax = np.min(one_over_k), np.max(one_over_k)
    x = np.linspace(xmin, xmax)
    y = model_delta_chi2(x, *coeff)
    # model results
    axis.plot(x, y, color="#B2ACAC")

    # sort the results
    s = np.argsort(one_over_k)
    x = [one_over_k[i] for i in s]
    y = [delta_chi2[i] for i in s]
    yerr = [std[i] for i in s]
    # data results
    if txt_dir.name.split("_")[0] == "NNPDF31":
        label = "NNPDF3.1"
        color = "#66c2a5"
    elif txt_dir.name.split("_")[0] == "PN3":
        label = "n3fit"
        color = "#fc8d62"
    elif txt_dir.name.split("_")[0] == "feature":
        label = "feature scaling"
        color = "#daa79a"
    else:
        label = "feature scaling"
        color = "#ffd92f"

    axis.errorbar(
        x, y, yerr=yerr, fmt=".--", label=label, color=color
    )  # , elinewidth=1.5, capsize=2., capthick=1.5)

    axis.set_title("$<\Delta\chi^2>$ computation vs. model - $N_{\t{rep}}$=%d" % nrep, fontsize=12)
    axis.grid(linewidth=0.2)
    axis.set_ylabel("$<\Delta\chi^2>$")
    axis.set_xlabel("$x$")

    return


def build_plots(files):
    n_rows = int(np.ceil(len(files) / 3))
    fig, ax = plt.subplots(n_rows, 3, figsize=[16, 8])  # , sharex='col')

    return fig, ax


parser = argparse.ArgumentParser()

parser.add_argument("neigs", help="number of eigenvectors used in the conversions")
parser.add_argument(
    "--num_plot", help="1: 100-125-150-175 nrep, 2: 200-250-300-350 nrep", choices={1, 2}, type=int
)
parser.add_argument("--save", help="use to save plots", action="store_true")

args = parser.parse_args()

# ugly.. I want the results ordered as NNPDF, n3fit, feature sc.
temp = [d for d in os.scandir("../dat/")]
temp = sorted(temp, key=lambda e: e.name)
# sort_dir = [temp[1], temp[2], temp[3], temp[0]]
sort_dir = [temp[1], temp[2], temp[0]]

for j, txt_dir in enumerate(sort_dir):

    # list all txt files with results
    # elements are os.DirEntry objects
    files = [f for f in os.scandir("../dat/" + txt_dir.name) if f.is_file()]
    # order files in ascending order
    files = sorted(files, key=lambda e: int(e.name.split("_")[-1][:-4]))

    if args.num_plot == 1:
        files = files[:6]
    elif args.num_plot == 2:
        files = files[6:]

    if j == 0:
        fig, ax = build_plots(files)

    if not txt_dir.name.endswith(args.neigs):
        continue

    # prepare array for plots of fit results
    num_of_sets = len(files)

    fit_x = np.zeros(num_of_sets, dtype=int)  # number of replicas
    # 4 coefficients nrep dependent for every 'numreplicas' set
    fit_y = np.zeros((num_of_sets, 4))
    # 4 uncertainties for every 'numreplicas' set
    fit_err = np.zeros((num_of_sets, 4))
    # 4 coefficients nrep independent for every 'numreplicas' set
    param = np.zeros((num_of_sets, 4))

    for i, f in enumerate(files):

        one_over_sigfrac, cv, std = np.loadtxt(
            "../dat/" + txt_dir.name + "/" + f.name, unpack=True, usecols=(1, 2, 3)
        )

        # number of replicas used in the current set
        nrep = int(f.name.split("_")[-1][:-4])

        # results of the fit -> coefficients nrep dependent to plot
        fit_x[i] = nrep
        # get parameters nrep dependent: fit <delta_chi2> as a function of 1/k
        fit_y[i, :], fit_err[i, :] = fit_delta_chi2(one_over_sigfrac, cv, std)

        # get parameters nrep independent
        param[i, :] = coeff_n_indep(*fit_y[i, :], nrep)

        plot_chi2(ax.reshape(-1)[i], one_over_sigfrac, cv, std, fit_y[i, :], nrep, txt_dir)


handles, labels = ax.reshape(-1)[i].get_legend_handles_labels()
handles.append(Line2D([0], [0], color="#B2ACAC", lw=1))
labels.append("Fit")
fig.legend(handles, labels, loc="upper center", fontsize="x-small")
# fig.suptitle('NNPDF3.1 - n3fit - feature scaling')# - %s eigenvectors' % args.neigs)

plt.tight_layout()
plt.subplots_adjust(top=0.85)

if args.save:
    path_save = "../../reports/results/"
    name_save = "chi2_modelcomparison_" + str(args.num_plot) + ".pdf"
    plt.savefig(path_save + name_save)

plt.show()
