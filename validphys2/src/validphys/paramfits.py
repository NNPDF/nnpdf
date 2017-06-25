"""
paramfits.py

Functionality to determine parameters from a scan over
PDF fits. αs is so far the only example.

The functions here are high level and specialized, and rely on the more low
level modules (e.g. fitdata.py and results.py) for most of the functionality.

They also need to work around the limitations in libnnpdf, and so the
performance may not be optimal.
"""

import numpy as np
import matplotlib.pyplot as plt
from reportengine.figure import figure

@figure
def plot_fits_as_profile(fits_pdfs, fits_total_chi2):
    """Plot the total central chi² as a function of the value of α_s.
    Note that this plots as a function of the key "AlphaS_MZ" in the LHAPDF
    file, which is annoyingly *not* α_s(MZ) for Nf<5."""
    fig, ax = plt.subplots()
    alphas = [pdf.AlphaS_MZ for pdf in fits_pdfs]
    ax.plot(alphas, fits_total_chi2)
    ax.set_xlabel(r'$\alpha_S$')
    ax.set_ylabel(r'$\chi²/N_{dat}$')
    return fig

@figure
def plot_fitted_replicas_as_profiles_matched(fits_pdfs,
        fits_replica_data_correlated, max_ndiscarded:int=4):
    """Plot chi²(as) keeping the replica nnfit index matched.

    The ``max_ndiscarded`` parameter defines th number of points
    discarded by postfit from which we discard the curve.
    """

    alphas = [pdf.AlphaS_MZ for pdf in fits_pdfs]
    df = fits_replica_data_correlated
    def ap(x):
        x.columns = x.columns.droplevel(0)
        return (x['training'] + x['validation'])/2
    table = df.groupby(axis=1, level=0).apply(ap)
    filt = table.isnull().sum(axis=1) < max_ndiscarded
    table = table[filt]
    table = table.as_matrix()

    minimums = []
    asarr = np.asarray(alphas)
    for row in table:
        filt =  np.isfinite(row)
        a,b,c = np.polyfit(asarr[filt], row[filt], 2)
        minimums.append(-b/2/a)
    minimums = np.asarray(minimums)

    fig, ax = plt.subplots()

    from matplotlib.collections import LineCollection

    lc = LineCollection([list(x for x in zip(alphas, t) if np.isfinite(x[1])) for t in table])
    lc.set_array(minimums)
    lc.set_clim(*np.percentile(minimums, (5,95)))
    ax.add_collection(lc)
    ax.set_xlim(min(alphas), max(alphas))
    ax.set_ylim(np.nanmin(table), np.nanmax(table))
    fig.colorbar(lc, label=r"Preferred $\alpha_S$")
    plt.title(rf"$\alpha_S$ from quadratic fit = ${np.mean(minimums):.4f} \pm {np.std(minimums):.4f}$ N={len(table)}")



    #ax.plot(alphas, np.array(table).T, color='#ddddcc')
    ax.set_xlabel(r'$\alpha_S$')
    ax.set_ylabel(r'$\chi²/N_{dat}$')
    return fig


@figure
def plot_poly_as_fit(fits_pdfs,
        fits_replica_data_correlated, max_ndiscarded:int=4, polorder:int=2):
    """Plot a polynomial fit of chi²(as) of `degree polorder`, keeping the
    replica index matched.

    The ``max_ndiscarded`` parameter defines th number of points
    discarded by postfit from which we discard the curve.
    """

    alphas = [pdf.AlphaS_MZ for pdf in fits_pdfs]
    df = fits_replica_data_correlated
    def ap(x):
        x.columns = x.columns.droplevel(0)
        return (x['chi2'])
    table = df.groupby(axis=1, level=0).apply(ap)
    filt = table.isnull().sum(axis=1) < max_ndiscarded
    table = table[filt]
    table = table.as_matrix()
    fig, ax = plt.subplots()

    minimums = []
    asarr = np.asarray(alphas)
    for i, row in enumerate(table):
        filt =  np.isfinite(row)
        fit = np.polyfit(asarr[filt], row[filt], polorder)
        extrema = np.roots(np.polyder(fit))

        #Cast away the zero complex part
        candidates = np.real(extrema[np.isreal(extrema)])
        candidates = candidates[(candidates > alphas[0]) & (candidates<alphas[-1])]
        if len(candidates):
            minimums.append(candidates[np.argmin(np.polyval(fit, candidates))])
            color = f'C{i%10}'
            ax.plot(alphas, np.polyval(fit, alphas), color=color)
            ax.plot(asarr[filt], row[filt], 'o', color=color)
    minimums = np.asarray(minimums)


    ax.set_xlim(min(alphas), max(alphas))
    ax.set_ylim(np.nanmin(table), np.nanmax(table))
    plt.title(rf"$\alpha_S$ from order {polorder} polynomial fit "
              rf"= ${np.mean(minimums):.4f} \pm {np.std(minimums):.4f}$"
              rf"N={len(table)}")
    ax.set_xlabel(r'$\alpha_S$')
    ax.set_ylabel(r'$\chi²/N_{dat}$')
    return fig
