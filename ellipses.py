#!/usr/bin/env python3
import posterior_mc_alphas
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
from matplotlib import rc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
plt.rcParams['text.usetex'] = True
# rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
# rc("text", **{"usetex": True, "latex.preamble": r"\usepackage{amssymb}"})

def plot_ellipses(fitnames, fitlabels):
    fig, axs = plt.subplots(nrows=1, ncols=len(fitnames), sharex=True, sharey=True, figsize=(10, 5))
    if len(fitnames)==1:
        axs = [axs]
    #fitnames_and_axs = list(zip(fitnames, axs))
    for i, ax in enumerate(axs):
        fitname=fitnames[i]
        fitlabel=fitlabels[i]
        print(f"results for {fitname}: \n")
        mean, cov = posterior_mc_alphas.theory_cov_method(fitname=fitname)
        posterior_mc_alphas.confidence_ellipse(ax, cov, mean, confidence_level=68)
        ax.set_xlabel(r'$m_\mathrm{charm}\mathrm{ [GeV]}$')
        ax.set_ylabel(r'$\alpha_\mathrm{s}$')
        ax.set_title(fitlabel.replace("_", " "))
    fig.set_label(r"TCM 68% confidence ellipses")
    fig.savefig("".join(fitnames)+".pdf")

if __name__=="__main__":
    fitnames = sys.argv[1].split(" ")
    fitlabels = sys.argv[2].split(" ")
    if len(fitnames) == len(fitlabels):
        plot_ellipses(fitnames, fitlabels)
    else:
        print("input lengths do not match")
        