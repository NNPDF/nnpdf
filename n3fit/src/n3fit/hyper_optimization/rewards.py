"""
Target functions to minimize during hyperparameter scan

These are implemented in the HyperLoss class which incorporates
various statistics (average, standard deviation, best/worst case)
both across multiple replicas of a model and across different folds.

Key functionalities include:
- Support for different loss types such as Chi-square (chi2) and phi-square (phi2).
- Calculation of statistical measures (average, best_worst, std) over replicas and folds.
- Incorporation of penalties into the loss computation.
- Detailed tracking and storage of loss metrics for further analysis.

New statistics can be added directly in this class as staticmethods and
via `IMPLEMENTED_STATS`; their name in the runcard must
match the name in the module

Example
-------
>>> import numpy as np
>>> from n3fit.hyper_optimization.rewards import HyperLoss
>>> losses = np.array([1.0, 2.0, 3.0])
>>> loss_average = HyperLoss(fold_statistic="average")
>>> loss_best_worst = HyperLoss(fold_statistic="best_worst")
>>> loss_std = HyperLoss(fold_statistic="std")
>>> print(f"{loss_average.reduce_over_folds.__name__} {loss_average.reduce_over_folds(losses)}")
>>> print(f"{loss_best_worst.reduce_over_folds.__name__} {loss_best_worst.reduce_over_folds(losses)}")
>>> print(f"{loss_std.reduce_over_folds.__name__} {loss_std.reduce_over_folds(losses)}")
_average 2.0
_best_worst 3.0
_std 0.816496580927726
"""

import logging
from typing import Callable

import numpy as np

from n3fit.vpinterface import N3PDF, HyperoptMetrics, compute_hyperopt_metrics
from validphys.core import DataGroupSpec
from validphys.pdfgrids import distance_grids, xplotting_grid

log = logging.getLogger(__name__)


def _average(fold_losses: np.ndarray, axis: int = 0, **kwargs) -> float:
    """
    Compute the average of the input array along the specified axis.

    Parameters
    ----------
        fold_losses: np.ndarray
            Input array.
        axis: int, optional
            Axis along which the mean is computed. Default is 0.

    Returns
    -------
        float: The average along the specified axis.
    """
    return np.average(fold_losses, axis=axis).item()


def _best_worst(fold_losses: np.ndarray, axis: int = 0, **kwargs) -> float:
    """
    Compute the maximum value of the input array along the specified axis.

    Parameters
    ----------
        fold_losses: np.ndarray
            Input array.
        axis: int, optional
            Axis along which the maximum is computed. Default is 0.

    Returns
    -------
        float: The maximum value along the specified axis.
    """
    return np.max(fold_losses, axis=axis).item()


def _std(fold_losses: np.ndarray, axis: int = 0, **kwargs) -> float:
    """
    Compute the standard deviation of the input array along the specified axis.

    Parameters
    ----------
        fold_losses: np.ndarray
            Input array.
        axis: int, optional
            Axis along which the standard deviation is computed. Default is 0.

    Returns
    -------
        float: The standard deviation along the specified axis.
    """
    return np.std(fold_losses, axis=axis).item()


IMPLEMENTED_STATS = {
    "average_best": _average,  # NB: all averages are average best now
    "average": _average,
    "best_worst": _best_worst,
    "std": _std,
}
IMPLEMENTED_LOSSES = ["chi2", "phi2", "logp", "chi2p"]


def _pdfs_to_n3pdfs(pdfs_per_fold):
    """Convert a list of multi-replica PDFs to a list of N3PDFs"""
    return [N3PDF(pdf.split_replicas(), name=f"fold_{k}") for k, pdf in enumerate(pdfs_per_fold)]


class HyperLoss:
    """
    Class to compute the hyper_loss based on the individual replica losses.

    Computes the statistic over the replicas and then over the folds, both
    statistics default to the average.

    The ``compute_loss`` method saves intermediate metrics such as the
    chi2 of the folds or the phi regardless of the loss type that has been selected.
    These metrics are saved in the properties
        ``phi2_vector``: list of phi per fold
        ``chi2_matrix``: list of chi2 per fold, per replica


    Parameters
    ----------
        loss_type: str
            the type of loss over the replicas to use.
            Options are "chi2" and "phi2".
        replica_statistic: str
            the statistic over the replicas to use, for per replica losses.
            Options are ``average``, ``best_worst`` and ``_std``.
        fold_statistic: str
            the statistic over the folds to use.
            Options are ``average``, ``best_worst`` and ``_std``.
        reduce_proportion: float (default 0.85)
            Proportion of replicas to select when computing statistics.
        penalties_in_loss: bool
            whether the penalties should be included in the output of ``compute_loss``
    """

    def __init__(
        self,
        loss_type: str = None,
        replica_statistic: str = None,
        fold_statistic: str = None,
        reduce_proportion: float = 0.85,
        penalties_in_loss: bool = False,
    ):
        self._default_loss = "chi2"
        self._penalties_in_loss = penalties_in_loss
        self._proportion = reduce_proportion

        self.loss_type = self._parse_loss(loss_type)

        self.hyper_chi2_vector = []
        self.hyper_phi2_vector = []
        self.hyper_logp_vector = []
        self.exp_chi2_matrix = []
        self.penalties = {}

        self.reduce_over_replicas = self._parse_statistic(replica_statistic, "replica")
        self.reduce_over_folds = self._parse_statistic(fold_statistic, "fold")

    def compute_loss(
        self,
        penalties: dict[str, np.ndarray],
        validation_loss: np.ndarray,
        experimental_loss: np.ndarray,
        pdf_object: N3PDF,
        experimental_data: list[DataGroupSpec],
        fold_idx: int = 0,
    ) -> float:
        """
        Compute the loss, including added penalties, for a single fold.

        Save the phi of the assemble and the chi2 of the separate replicas,
        and the penalties into the ``phi2_vector``, ``chi2_matrix`` and ``penalties`` attributes.

        Parameters
        ----------
            penalties: Dict[str, NDArray(replicas)]
                Dict of penalties for each replica.
                Possible keys are 'saturation', 'patience' and 'integrability'
                as defined in 'penalties.py' and instantiated within :class:`~n3fit.model_trainer.ModelTrainer`.
            experimental_loss: NDArray(replicas)
                Experimental loss for each replica.
            pdf_object: :class:`n3fit.vpinterface.N3PDF`
                N3fitted PDF
            experimental_data: List[validphys.core.DataGroupSpec]
                List of tuples containing `validphys.core.DataGroupSpec` instances for each group data set
            fold_idx: int
                k-fold index. Defaults to 0.

        Returns
        -------
            loss: float
                The computed loss over the replicas.

        Example
        -------
        >>> import numpy as np
        >>> from n3fit.hyper_optimization.rewards import HyperLoss
        >>> from n3fit.model_gen import generate_pdf_model
        >>> from n3fit.vpinterface import N3PDF
        >>> from validphys.loader import Loader
        >>> hyper = HyperLoss(loss_type="chi2", replica_statistic="average", fold_statistic="average")
        >>> penalties = {'saturation': np.array([1.0, 2.0]), 'patience': np.array([3.0, 4.0]), 'integrability': np.array([5.0, 6.0]),}
        >>> experimental_loss = np.array([0.1, 0.2])
        >>> ds = Loader().check_dataset("NMC_NC_NOTFIXED_P_EM-SIGMARED", variant="legacy", theoryid=399, cuts="internal")
        >>> experimental_data = [Loader().check_experiment("My DataGroupSpec", [ds])]
        >>> fake_fl = [{'fl' : i, 'largex' : [0,1], 'smallx': [1,2]} for i in ['u', 'ubar', 'd', 'dbar', 'c', 'g', 's', 'sbar']]
        >>> pdf_model = generate_pdf_model(nodes=[8], activations=['linear'], seed=[0,2], flav_info=fake_fl, fitbasis="FLAVOUR")
        >>> pdf = N3PDF(pdf_model.split_replicas())
        >>> loss = hyper.compute_loss(penalties, experimental_loss, pdf, experimental_data)
        """
        if np.isnan(validation_loss).any():
            log.warning(f"{np.isnan(validation_loss).sum()} replicas have NaNs losses")

        # Before starting, select the best replicas according to the proportion set in the __init__
        num_best = int(np.ceil(self._proportion * len(validation_loss)))
        best_indexes = np.argsort(validation_loss, axis=0)[:num_best]
        best_validation_losses = validation_loss[best_indexes]

        # Select the `N3PDF` models to be used to compute the hyperopt metrics. The models
        # are selected based on the validation losses using `self._proportion`.
        pdf_object_reduced = pdf_object.select_models(best_indexes)

        # Compute the different hyperopt metrics
        hypermetics: HyperoptMetrics = compute_hyperopt_metrics(
            n3pdf=pdf_object_reduced, experimental_data=experimental_data
        )

        # Extract & save the values of the hyperopt metrics
        hyper_chi2_per_fold = hypermetics.chi2  # computed with PDF covmat
        hyper_phi2_per_fold = hypermetics.phi2  # computed without PDF covmat
        hyper_logp_per_fold = hypermetics.logp  # computed with PDF covmat

        # Update hyperopt metrics history
        self._save_hyperopt_metrics(
            hyper_chi2_per_fold,
            hyper_phi2_per_fold,
            hyper_logp_per_fold,
            experimental_loss,
            penalties,
            fold_idx,
        )

        # Prepare the output loss, including penalties if necessary
        if self._penalties_in_loss:
            # include penalties to experimental loss
            experimental_loss += sum(penalties.values())
            # add penalties to `phi2` and `logp` in the form of a sum of per-replicas averages
            sum_penalties = sum(np.mean(penalty) for penalty in penalties.values())
            hyper_phi2_per_fold += sum_penalties
            hyper_logp_per_fold += sum_penalties

        # define loss for hyperopt according to the chosen loss_type
        if self.loss_type == "chi2":
            # calculate statistics of chi2 over replicas for a given k-fold_statistic

            # Construct the final loss as a sum of:
            # 1. The validation chi2
            # 2. The distance to 2 for the experimental chi2
            # In the hyperopt paper we used 80% and 10% respectively, as a proxy for:
            # "80% of the replicas should be good, but only a small % has to cover the folds"
            # Currently take reduce_proportion for a) and 1.0 - reduce_proportion for b)
            validation_loss_average = self.reduce_over_replicas(best_validation_losses)

            nselect = int(np.ceil((1.0 - self._proportion) * len(experimental_loss)))
            best_exp_losses = np.sort(experimental_loss, axis=0)[:nselect]
            exp_loss_average = self.reduce_over_replicas(best_exp_losses)

            loss = validation_loss_average + (max(exp_loss_average, 2.0) - 2.0)
        elif self.loss_type == "phi2":
            loss = hyper_phi2_per_fold
        elif self.loss_type == "logp":
            loss = hyper_logp_per_fold
        elif self.loss_type == "chi2p":
            loss = hyper_chi2_per_fold

        return loss

    def _save_hyperopt_metrics(
        self,
        hyper_chi2_per_fold: float,
        hyper_phi2_per_fold: float,
        hyper_logp_per_fold: float,
        exp_chi2_per_fold: np.ndarray,
        penalties: dict[str, np.ndarray],
        fold_idx: int = 0,
    ) -> None:
        """
        Save all the calculated metrics per replica and per fold, including penalties.

        Parameters
        ----------
            hyper_chi2_per_fold: float
                Computed chi2 for a given k-fold
            hyper_phi2_per_fold: float
                Computed phi2 for a given k-fold
            hyper_logp_per_fold: float
                Computed logp for a given k-fold
            exp_chi2_per_fold: np.ndarray
                Computed experimental chi2 for all the replica for a given k-fold
            penalties: Dict[str, np.ndarray]
                dictionary of all penalties with their names
            fold_idx: int
                k-fold index. Defaults to 0.
        """
        # reset chi2 and phi arrays for every trial
        if fold_idx == 0:
            self.hyper_chi2_vector = []
            self.hyper_phi2_vector = []
            self.hyper_logp_vector = []
            self.exp_chi2_matrix = []
            self.penalties = {}

        # populate chi2 matrix and phi vector calculated for a given k-fold
        self.hyper_chi2_vector.append(hyper_chi2_per_fold)
        self.hyper_phi2_vector.append(hyper_phi2_per_fold)
        self.hyper_logp_vector.append(hyper_logp_per_fold)
        self.exp_chi2_matrix.append(exp_chi2_per_fold)

        # save penalties per replica for a given k-fold
        for name, values in penalties.items():
            temp = self.penalties.get(name, [])
            temp.append(values)
            self.penalties[name] = temp

    def _parse_loss(self, loss_type: str) -> str:
        """
        Parse the type of loss and return the default if None.

        Parameters
        ----------
            loss_type: str
                The loss type to parse.

        Returns
        -------
            loss_type: str
                The parsed loss type.

        Raises
        ------
            ValueError: If an invalid loss type is provided.
        """
        if loss_type is None:
            loss_type = self._default_loss
            log.warning(f"No loss_type selected in HyperLoss, defaulting to {loss_type}")

        if loss_type not in IMPLEMENTED_LOSSES:
            valid_options = ", ".join(IMPLEMENTED_LOSSES)
            raise ValueError(f"Invalid loss type '{loss_type}'. Valid options are: {valid_options}")

        log.info(f"Setting '{loss_type}' as the loss type for hyperoptimization")

        return loss_type

    def _parse_statistic(self, statistic: str, target: str, default: str = "average") -> Callable:
        """
        Parse the statistic and return the default if None.


        Parameters
        ----------
            statistic: str
                The statistic to parse.
            target: str
                The target of the statistic (either replica or fold)

        Returns
        -------
            Callable: The parsed statistic method.

        Raises
        ------
            ValueError: If an invalid statistic is provided.

        Notes
        -----
            For loss type equal to phi2, the applied fold statistics is always the reciprocal of the selected stats.
        """
        if statistic is None:
            statistic = default
            log.warning(f"No {target} selected in HyperLoss, defaulting to {statistic}")

        if statistic not in IMPLEMENTED_STATS:
            valid_options = ", ".join(IMPLEMENTED_STATS.keys())
            raise ValueError(f"Invalid {target} '{statistic}'. Valid options are: {valid_options}")

        log.info(f"Using '{statistic}' as the {target} for hyperoptimization")

        selected_statistic = IMPLEMENTED_STATS[statistic]

        if self.loss_type == "chi2" or self.loss_type == "logp":
            return selected_statistic
        elif self.loss_type == "phi2":
            # In case of phi2, calculate the inverse of the applied statistics
            # This is only used when calculating statistics over folds
            return lambda x: np.reciprocal(selected_statistic(x))

        raise ValueError(f"{self.loss_type} is not a valid hyperopt loss.")


def fit_distance(pdfs_per_fold=None, **_kwargs):
    """Loss function for hyperoptimization based on the distance of
    the fits of all folds to the first fold
    """
    n3pdfs = _pdfs_to_n3pdfs(pdfs_per_fold)
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
        if key != "xgrid_integration":
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
