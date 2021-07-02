"""
    Target functions to minimize during hyperparameter scan

    All functions in this module have the same signature:
        fold_losses: list of loss-per-fold
        **kwargs

    New loss functions can be added directly in this module
    the name in the runcard must match the name in the module

    Example
    -------
    >>> import n3fit.hyper_optimization.rewards
    >>> f = ["average", "best_worst", "std"]
    >>> losses = [2.34, 1.234, 3.42]
    >>> for fname in f:
    >>>    fun = getattr(n3fit.hyper_optimization.rewards, fname) 
    >>>    print(f"{fname}: {fun(losses, None):2.4f}")
    average: 2.3313
    best_worst: 3.4200
    std: 0.8925

"""
import numpy as np
from validphys.pdfgrids import xplotting_grid, distance_grids


def average(fold_losses, n3pdf_objects, **kwargs):
    """ Returns the average of fold losses """
    return np.average(fold_losses)


def best_worst(fold_losses, n3pdf_objects, **kwargs):
    """ Returns the maximum loss of all k folds """
    return np.max(fold_losses)


def std(fold_losses, n3pdf_objects, **kwargs):
    """ Return the standard dev of the losses of the folds """
    return np.std(fold_losses)

def fit_distance(fold_losses, n3pdf_objects):
    """ Loss function for hyperoptimization based on the distance of
    the fits of all folds to the first fold
    """
    xgrid = np.concatenate([np.logspace(-6, -1, 20), np.linspace(0.11, 0.9, 30)])
    plotting_grids = [xplotting_grid(pdf, 1.6, xgrid) for pdf in n3pdf_objects]
    distances = distance_grids(n3pdf_objects, plotting_grids, 0)
    # The first distance will obviously be 0
    # TODO: define this more sensibly, for now it is just a template
    max_distance = 0
    for distance in distances:
        max_distance = max(max_distance, distance.grid_values.max())
    return max_distance
