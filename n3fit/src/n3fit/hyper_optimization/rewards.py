"""
    Target functions to minimize during hyperparameter scan

    Not all functions will use all arguments.
    Keyword arguments that model_trainer.py will pass to this file are:

    - fold_losses: a list with the loss of each fold
    - n3pdfs: a list of N3PDF objects for each fit (which can contain more than 1 replica)
    - experimental_models: a reference to the model that contains the cv for all data (no masks)

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
import logging

import numpy as np

from n3fit.vpinterface import N3PDF, compute_phi2
from validphys.pdfgrids import distance_grids, xplotting_grid

log = logging.getLogger(__name__)


class HyperLoss:
    """
    Class to compute the hyper_loss based on the individual replica losses.

    Computes the statistic over the replicas and then over the folds, both
    statistics default to the average.

    Args:
        loss (str): the loss over the replicas to use
        fold_statistic (str): the statistic over the folds to use
        replica_statistic (str): the statistic over the replicas to use, for per replica losses
    """

    def __init__(self, loss: str = None, fold_statistic: str = None):
        self.implemented_stats = {
            "average": self._average,
            "best_worst": self._best_worst,
            "std": self._std,
        }
        self._default_statistic = "average"

        self.loss = loss
        self.implemented_losses = {
            "chi2": self._chi2,
            "phi2": self._phi2,
        }
        self._default_losss = "chi2"

        self.reduce_over_folds = self._parse_statistic(fold_statistic, "fold_statistic")
        self.reduce_over_replicas = self._parse_statistic(replica_statistic, "replica_statistic")

    def compute_loss(self, penalties, experimental_loss, pdf_models, experimental_data):
        """
        Compute the loss, including added penalties, for a single fold.

        Args:
            penalties: List(NDArray(replicas))
            experimental_loss: NDArray(replicas)
            pdf_models: List(MetaModel)
            experimental_data: List(dict)
                values:
                covariances:

        Returns:
            float
        """
        total_penalties = sum(np.mean(penalty) for penalty in penalties)

        if self.loss == "chi2":
            loss = self.reduce_over_replicas(experimental_loss)
        elif self.loss == "phi2":
            loss = compute_phi2(N3PDF(pdf_models), experimental_data)

        return total_penalties + loss

    def _parse_statistic(self, statistic: str, name) -> str:
        """
        Parse the statistic and return the default if None.

        Args:
            statistic (str): the statistic to parse
            name (str): the name of the statistic

        Returns:
            str: the parsed statistic
        """
        if statistic is None:
            statistic = self._default_statistic
            log.warning(f"No {name} selected in HyperLoss, defaulting to {statistic}")
        log.info(f"Using '{statistic}' as the {name} for hyperoptimization")
        return self.implemented_stats[statistic]

    @staticmethod
    def _average(fold_losses: np.ndarray, axis: int = 0) -> np.ndarray:
        return np.average(fold_losses, axis=axis).item()

    @staticmethod
    def _best_worst(fold_losses: np.ndarray, axis: int = 0) -> np.ndarray:
        return np.max(fold_losses, axis=axis).item()

    @staticmethod
    def _std(fold_losses: np.ndarray, axis: int = 0) -> np.ndarray:
        return np.std(fold_losses, axis=axis).item()


def fit_distance(n3pdfs=None, **_kwargs):
    """Loss function for hyperoptimization based on the distance of
    the fits of all folds to the first fold
    """
    if n3pdfs is None:
        raise ValueError("fit_distance needs n3pdf models to act upon")
    xgrid = np.concatenate([np.logspace(-6, -1, 20), np.linspace(0.11, 0.9, 30)])
    plotting_grids = [xplotting_grid(pdf, 1.6, xgrid) for pdf in n3pdfs]
    distances = distance_grids(n3pdfs, plotting_grids, 0)
    # The first distance will obviously be 0
    # TODO: define this more sensibly, for now it is just a template
    max_distance = 0
    for distance in distances:
        max_distance = max(max_distance, distance.grid_values.max())
    return max_distance


def _set_central_value(n3pdf, model):
    """Given an N3PDF object and a MetaModel, set the PDF_0 layer
    to be the central value of the n3pdf object"""
    from n3fit.backends import operations as op

    # Get the input x
    for key, grid in model.x_in.items():
        if key != "integration_grid":
            input_x = grid.numpy()
            break

    # Compute the central value of the PDF
    full_pdf = n3pdf(input_x)
    cv_pdf = op.numpy_to_tensor(np.mean(full_pdf, axis=0, keepdims=True))

    def central_value(x, training=None):  # pylint: disable=unused-argument
        return cv_pdf

    model.get_layer("PDF_0").call = central_value
    # This model won't be trainable ever again
    model.trainable = False
    model.compile()


def fit_future_tests(n3pdfs=None, experimental_models=None, **_kwargs):
    """Use the future tests as a metric for hyperopt

    NOTE: this function should only be called once at the end of
    every hyperopt iteration, as it is destructive for the models
    """
    if n3pdfs is None:
        raise ValueError("fit_future_test needs n3pdf models to act upon")
    if experimental_models is None:
        raise ValueError("fit_future_test needs experimental_models to compute chi2")

    from n3fit.backends import MetaModel

    compatibility_mode = False
    try:
        import tensorflow as tf

        from n3fit.backends import set_eager

        tf_version = tf.__version__.split(".")
        if int(tf_version[0]) == 2 and int(tf_version[1]) < 4:
            set_eager(True)
            compatibility_mode = True
    except ImportError:
        pass

    # For the last model the PDF covmat doesn't need to be computed.
    # This is because the last model corresponds to an empty k-fold partition,
    # meaning that all datasets were used during training.
    # but the mask needs to be flipped in the folding for the appropiate datasets
    last_model = experimental_models[-1]
    _set_central_value(n3pdfs[-1], last_model)

    # Loop over all models but the last (our reference!)
    total_loss = 0.0
    for n3pdf, exp_model in zip(n3pdfs[:-1], experimental_models[:-1]):
        _set_central_value(n3pdf, exp_model)

        # Get the full input and the total chi2
        full_input = exp_model.input

        # Now update the loss with the PDF covmat
        for layer in exp_model.get_layer_re(".*_exp$"):
            # Get the input to the loss layer
            model_output = MetaModel(full_input, layer.input)
            # Get the full predictions and generate the PDF covmat
            y = model_output.predict()[0]  # The first is a dummy-dim
            pdf_covmat = np.cov(y, rowvar=False)
            # Update the covmat of the loss
            layer.add_covmat(pdf_covmat)
            # Update the mask of the last_model so that its synced with this layer
            last_model.get_layer(layer.name).update_mask(layer.mask)

        # Compute the loss with pdf errors
        pdf_chi2 = exp_model.compute_losses()["loss"][0]

        # And the loss of the best (most complete) fit
        best_chi2 = last_model.compute_losses()["loss"][0]

        # Now make this into a measure of the total loss
        # for instance, any deviation from the "best" value is bad
        total_loss += np.abs(best_chi2 - pdf_chi2)

    if compatibility_mode:
        set_eager(False)

    return total_loss
