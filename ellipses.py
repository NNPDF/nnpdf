#!/usr/bin/env python3
import posterior_mc_alphas
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
from matplotlib import rc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
from pathlib import PosixPath
plt.rcParams['text.usetex'] = True
# rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
# rc("text", **{"usetex": True, "latex.preamble": r"\usepackage{amssymb}"})

def plot_ellipses(fitnames, fitlabels, save_directory: PosixPath):
    fig, axs = plt.subplots(nrows=1, ncols=len(fitnames), sharex=True, sharey=True, figsize=(10, 5))
    if len(fitnames)==1:
        axs = [axs]
    #fitnames_and_axs = list(zip(fitnames, axs))
    for i, ax in enumerate(axs):
        fitname=fitnames[i]
        fitlabel=fitlabels[i]
        print(f"results for {fitname}: \n")
        mean, cov, str_result = posterior_mc_alphas.theory_cov_method(fitname=fitname)
        posterior_mc_alphas.confidence_ellipse(ax, cov, mean, confidence_level=68)
        ax.set_xlabel(r'$m_\mathrm{charm}\mathrm{ [GeV]}$')
        ax.set_ylabel(r'$\alpha_\mathrm{s}$')
        ax.set_title(fitlabel.replace("_", " "))
        with open(save_directory / str("TCM"+fitname+".txt"), "w") as f:
            f.write(str_result)
    fig.set_label(r"TCM 68% confidence ellipses")
    fig.savefig(save_directory / str("".join(fitnames)+".pdf"))

if __name__=="__main__":
    fitnames = sys.argv[1].split(" ")
    fitlabels = sys.argv[2].split(" ")
    save_directory = PosixPath(sys.argv[3])
    if len(fitnames) == len(fitlabels):
        plot_ellipses(fitnames, fitlabels, save_directory)
    else:
        print("input lengths do not match")
        