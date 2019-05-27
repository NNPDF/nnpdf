"""
actions.py

Basic tools to study the IRIS dataset.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, hamming_loss
from sklearn.model_selection import train_test_split

from reportengine import collect
from reportengine.figure import figure
from reportengine.table import table
from reportengine.checks import make_argcheck, CheckError

def fit_result(algorithm, dataset):
    """Fit to a sample of the dataset where some labels have been
    shuffled  using the default algorithm
    parameters"""
    y_train = dataset.target.copy()
    np.random.seed(1520)
    mask = np.random.randint(len(y_train), size=len(y_train)//2)
    to_shuffle = y_train[mask]
    np.random.shuffle(to_shuffle)
    y_train[mask] = to_shuffle
    return algorithm().fit(dataset.data, y_train)

@figure
def plot_2d(scatterdata):
    """Generate scatter plot of the values of xaxis vs yaxis"""
    x,y,category = scatterdata
    fig, ax = plt.subplots()
    ax.scatter(x,y, c=category, cmap=plt.cm.coolwarm)
    return fig

@make_argcheck
def _check_can_predict_probabilities(algorithm):
    res = hasattr(algorithm().fit([[0],[1]], [0,1]), 'predict_proba')
    if not res:
        raise CheckError(f"Algorithm {algorithm.__name__} doesn't support "
                "'predict_proba'")


@make_argcheck
def _check_fpr_threshold(fpr_threshold):
    """Check that it's in (0,1) if given"""
    if fpr_threshold is None:
        return
    if not 0<fpr_threshold<1:
        raise CheckError('fpr_threshold must be contained in (0,1)')

@figure
@_check_fpr_threshold
@_check_can_predict_probabilities
def plot_roc(fit_result, algorithm, dataset,
        fpr_threshold:(float, type(None))=None):
    """Plot the ROC curve for each category. Mark the true positive
    rate at the ``fpr_threshold`` if given"""
    probs = fit_result.predict_proba(dataset.data)
    fig, ax = plt.subplots()
    for i, label in enumerate(dataset.target_names):
        y_pred = probs[:,i]
        y_true = (dataset.target == i)
        fpr, tpr, _ = roc_curve(y_true, y_pred)
        color = f'C{i}'
        ax.plot(fpr, tpr, label=label, color=color)
        if fpr_threshold is not  None:
            pos = np.searchsorted(fpr, fpr_threshold)
            ax.axhline(tpr[pos], linestyle='--', lw=0.5, color=color)

    ax.set_xlabel("False Positive rate")
    ax.set_ylabel("True Positive rate")
    ax.set_title("ROC curve")
    ax.legend()

    return fig

fit_results = collect(fit_result, ['algorithms'])

@table
def hamming_loss_table(algorithms, fit_results, dataset):
    records = []
    for algorithm, res in zip(algorithms, fit_results):
        records.append({'algorithm': algorithm.__name__, "Hamming "
                "loss": hamming_loss(dataset.target,
                    res.predict(dataset.data))})
    return pd.DataFrame(records)

