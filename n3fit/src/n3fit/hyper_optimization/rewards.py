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
    >>>    print(f"{fname}: {fun(losses):2.4f}")
    average: 2.3313
    best_worst: 3.4200
    std: 0.8925

"""
import numpy as np


def average(fold_losses, **kwargs):
    """ Returns the average of fold losses """
    return np.average(fold_losses)


def best_worst(fold_losses, **kwargs):
    """ Returns the maximum loss of all k folds """
    return np.max(fold_losses)


def std(fold_losses, **kwargs):
    """ Return the standard dev of the losses of the folds """
    return np.std(fold_losses)
