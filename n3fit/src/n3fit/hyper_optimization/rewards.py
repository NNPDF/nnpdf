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
    via :attr:`~HyperLoss.implemented_stats`; their name in the runcard must
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
from typing import Callable, Dict, List

import numpy as np

from n3fit.backends import MetaModel
from n3fit.vpinterface import N3PDF, compute_phi2
from validphys.core import DataGroupSpec
from validphys.pdfgrids import distance_grids, xplotting_grid

log = logging.getLogger(__name__)


def _pdfs_to_n3pdfs(pdfs_per_fold):
    """Convert a list of multi-replica PDFs to a list of N3PDFs"""
    return [N3PDF(pdf.split_replicas(), name=f"fold_{k}") for k, pdf in enumerate(pdfs_per_fold)]


class HyperLoss:
    """
    Class to compute the hyper_loss based on the individual replica losses.

    Computes the statistic over the replicas and then over the folds, both
    statistics default to the average.

    Parameters
    ----------
        loss_type: str
            the type of loss over the replicas to use.
            Options are "chi2" and "phi2".
        replica_statistic: str
            the statistic over the replicas to use, for per replica losses.
            Options are "average", "best_worst", and "std".
        fold_statistic: str
            the statistic over the folds to use.
            Options are "average", "best_worst", and "std".
    """

    def __init__(
        self, loss_type: str = None, replica_statistic: str = None, fold_statistic: str = None
    ):
        self.implemented_stats = {
            "average": self._average,
            "best_worst": self._best_worst,
            "std": self._std,
        }
        self.implemented_losses = ["chi2", "phi2"]

        self._default_statistic = "average"
        self._default_loss = "chi2"

        self.loss_type = self._parse_loss(loss_type)
        self.reduce_over_replicas = self._parse_statistic(replica_statistic, "replica_statistic")
        self.reduce_over_folds = self._parse_statistic(fold_statistic, "fold_statistic")

        self.phi2_vector = []
        self.chi2_matrix = []

        self.penalties = {}

    def compute_loss(
        self,
        penalties: Dict[str, np.ndarray],
        experimental_loss: np.ndarray,
        pdf_model: MetaModel,
        experimental_data: List[DataGroupSpec],
        fold_idx: int = 0,
    ) -> float:
        """
        Compute the loss, including added penalties, for a single fold.

        Parameters
        ----------
            penalties: Dict[str, NDArray(replicas)]
                Dict of penalties for each replica.
                Possible keys are 'saturation', 'patience' and 'integrability'
                as defined in 'penalties.py' and instantiated within :class:`~n3fit.model_trainer.ModelTrainer`.
            experimental_loss: NDArray(replicas)
                Experimental loss for each replica.
            pdf_model: :class:`n3fit.backends.MetaModel`
                N3fitted meta-model.
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
        >>> from validphys.loader import Loader
        >>> hyper = HyperLoss(loss_type="chi2", replica_statistic="average", fold_statistic="average")
        >>> penalties = {'saturation': np.array([1.0, 2.0]), 'patience': np.array([3.0, 4.0]), 'integrability': np.array([5.0, 6.0]),}
        >>> experimental_loss = np.array([0.1, 0.2])
        >>> ds = Loader().check_dataset("NMC", theoryid=399, cuts="internal")
        >>> experimental_data = [Loader().check_experiment("My DataGroupSpec", [ds])]
        >>> fake_fl = [{'fl' : i, 'largex' : [0,1], 'smallx': [1,2]} for i in ['u', 'ubar', 'd', 'dbar', 'c', 'g', 's', 'sbar']]
        >>> pdf_model = generate_pdf_model(nodes=[8], activations=['linear'], seed=0, num_replicas=2, flav_info=fake_fl, fitbasis="FLAVOUR")
        >>> loss = hyper.compute_loss(penalties, experimental_loss, pdf_model, experimental_data)
        """
        # calculate phi2 for a given k-fold using vpinterface and validphys
        phi2_per_fold = compute_phi2(N3PDF(pdf_model.split_replicas()), experimental_data)

        # update hyperopt metrics
        # these are saved in the phi2_vector and chi2_matrix attributes, excluding penalties
        self._save_hyperopt_metrics(phi2_per_fold, experimental_loss, penalties, fold_idx)

        # include penalties to experimental loss
        # this allows introduction of statistics also to penalties
        experimental_loss_w_penalties = experimental_loss + sum(penalties.values())

        # add penalties to phi2 in the form of a sum of per-replicas averages
        phi2_per_fold += sum(np.mean(penalty) for penalty in penalties.values())

        # define loss for hyperopt according to the chosen loss_type
        if self.loss_type == "chi2":
            # calculate statistics of chi2 over replicas for a given k-fold
            loss = self.reduce_over_replicas(experimental_loss_w_penalties)
        elif self.loss_type == "phi2":
            loss = phi2_per_fold

        return loss

    def _save_hyperopt_metrics(
        self,
        phi2_per_fold: float,
        chi2_per_fold: np.ndarray,
        penalties: Dict[str, np.ndarray],
        fold_idx: int = 0,
    ) -> None:
        """
        Save all chi2 and phi2 calculated metrics per replica and per fold, including penalties.

        Parameters
        ----------
            phi2_per_fold: float
                Computed phi2 for a given k-fold
            chi2_per_fold: np.ndarray
                Computed chi2 for each replica for a given k-fold
            penalties: Dict[str, np.ndarray]
                dictionary of all penalties with their names
            fold_idx: int
                k-fold index. Defaults to 0.
        """
        # reset chi2 and phi2 arrays for every trial
        if fold_idx == 0:
            self.phi2_vector = []
            self.chi2_matrix = []
            self.penalties = {}

        # populate chi2 matrix and phi2 vector calculated for a given k-fold
        self.chi2_matrix.append(chi2_per_fold)
        self.phi2_vector.append(phi2_per_fold)

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
        else:
            if loss_type not in self.implemented_losses:
                valid_options = ", ".join(self.implemented_losses)
                raise ValueError(
                    f"Invalid loss type '{loss_type}'. Valid options are: {valid_options}"
                )

        log.info(f"Setting '{loss_type}' as the loss type for hyperoptimization")

        return loss_type

    def _parse_statistic(self, statistic: str, name: str) -> Callable:
        """
        Parse the statistic and return the default if None.

        Parameters
        ----------
            statistic: str
                The statistic to parse.
            name: str
                The name of the statistic.

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
            statistic = self._default_statistic
            log.warning(f"No {name} selected in HyperLoss, defaulting to {statistic}")
        else:
            if statistic not in self.implemented_stats:
                valid_options = ", ".join(self.implemented_stats.keys())
                raise ValueError(
                    f"Invalid {name} '{statistic}'. Valid options are: {valid_options}"
                )

        log.info(f"Using '{statistic}' as the {name} for hyperoptimization")

        selected_statistic = self.implemented_stats[statistic]

        if self.loss_type == "chi2":
            return selected_statistic

        elif self.loss_type == "phi2":
            # In case of phi2, calculate the inverse of the applied statistics
            # This is only used when calculating statistics over folds
            return lambda x: np.reciprocal(selected_statistic(x))

    @staticmethod
    def _average(fold_losses: np.ndarray, axis: int = 0) -> float:
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

    @staticmethod
    def _best_worst(fold_losses: np.ndarray, axis: int = 0) -> float:
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

    @staticmethod
    def _std(fold_losses: np.ndarray, axis: int = 0) -> float:
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
