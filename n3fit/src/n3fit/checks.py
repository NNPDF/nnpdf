"""
This module contains checks to be perform by n3fit on the input
"""

import logging
import numbers
import os

from n3fit.hyper_optimization import penalties as penalties_module
from n3fit.hyper_optimization.rewards import IMPLEMENTED_LOSSES, IMPLEMENTED_STATS
from reportengine.checks import CheckError, make_argcheck
from validphys.loader import FallbackLoader
from validphys.pdfbases import check_basis

log = logging.getLogger(__name__)

NN_PARAMETERS = ["nodes_per_layer", "optimizer", "activation_per_layer"]


def _is_floatable(num):
    """Check that num is a number or, worst case scenario, a number that can
    be casted to a float (such as a tf scalar)"""
    if isinstance(num, numbers.Number):
        return True
    try:
        float(num)
        return True
    except (ValueError, TypeError):
        return False


# Checks on the NN parameters
def check_existing_parameters(parameters):
    """Check that non-optional parameters are defined and are not empty"""
    for param_name in NN_PARAMETERS:
        if param_name in parameters:
            val = parameters[param_name]
            if len(val) == 0:
                raise CheckError(f"The parameter {param_name} cannot be empty")
        else:
            raise CheckError(f"Missing {param_name} parameter in the runcard")
    # Check that positivity is not defined wrong
    if "pos_initial" in parameters or "pos_multiplier" in parameters:
        raise CheckError(
            "The definition of the positivity parameters is deprecated, please "
            "use instead:\nparameters:\n  positivity\n"
            "    multiplier: x\n    initial: y\n"
            "as can be seen in the example runcard: n3fit/runcards/Basic_runcard.yml"
        )


def check_consistent_layers(parameters):
    """Checks that all layers have an activation function defined
    and that a final-activation function is not being used half-way through.
    """
    final_activations = ["square_singlet"]
    npl = len(parameters["nodes_per_layer"])
    act_per_layer = parameters["activation_per_layer"]
    apl = len(act_per_layer)
    if npl != apl:
        raise CheckError(f"Number of layers ({npl}) does not match activation functions: {apl}")

    for fin_act in final_activations:
        if fin_act in act_per_layer[:-1]:
            raise CheckError(f"The activation {fin_act} can only be used as last layer")


def check_stopping(parameters):
    """Checks whether the stopping-related options are sane:
    stopping patience as a ratio between 0 and 1
    and positive number of epochs
    """
    spt = parameters.get("stopping_patience")
    if spt is not None and not 0.0 <= spt <= 1.0:
        raise CheckError(f"The stopping_patience must be between 0 and 1, got: {spt}")
    epochs = parameters["epochs"]
    if epochs < 1:
        raise CheckError(f"Needs to run at least 1 epoch, got: {epochs}")


def check_basis_with_layers(basis, validphys_basis, parameters):
    """Check that the last layer matches the number of flavours defined in the runcard.
    And that the activation functions are compatible with the basis.
    """
    number_of_flavours = len(basis)
    last_layer = parameters["nodes_per_layer"][-1]
    if number_of_flavours != last_layer:
        raise CheckError(
            f"The number of nodes in the last layer ({last_layer}) does not"
            f" match the number of flavours: ({number_of_flavours})"
        )

    flavours = [i["fl"] for i in basis]
    if parameters["activation_per_layer"][-1] == "square_singlet":
        if not (("sng" in flavours) and ("g" in flavours)):
            raise CheckError(
                "square_singlet can only be used when `gluon` (g) and `singlet` (sng) are being fitted"
            )
        if (val := validphys_basis.indexes.get("sng")) > 1:
            raise CheckError(
                f"When using square_singlet, \\Sigma must be either element 0 or 1, found {val}"
            )
        if (val := validphys_basis.indexes.get("g")) > 1:
            raise CheckError(
                f"When using square_singlet, gluon must be either element 0 or 1, found {val}"
            )


def check_optimizer(optimizer_dict):
    """Checks whether the optimizer setup is valid"""
    name_key = "optimizer_name"
    name = optimizer_dict[name_key]
    from n3fit.backends import MetaModel

    accepted_optimizers = MetaModel.accepted_optimizers
    optimizer_data = accepted_optimizers.get(name)
    if optimizer_data is None:
        raise CheckError(f"Optimizer {name} not accepted by MetaModel")
    # Get the dictionary of accepted parameters
    data = optimizer_data[1]
    for key in optimizer_dict.keys():
        if key not in data and key != name_key:
            raise CheckError(f"Optimizer {name} does not accept the option: {key}")


def check_initializer(initializer):
    """Checks whether the initializer is implemented"""
    from n3fit.backends import MetaLayer

    accepted_init = MetaLayer.initializers
    if initializer not in accepted_init:
        raise CheckError(f"Initializer {initializer} not accepted by {MetaLayer}")


def check_layer_type_implemented(parameters):
    """Checks whether the layer_type is implemented"""
    layer_type = parameters.get("layer_type")
    implemented_types = ["dense", "dense_per_flavour"]
    if layer_type not in implemented_types:
        raise CheckError(
            f"Layer type {layer_type} not implemented, must be one of {implemented_types}"
        )


def check_dropout(parameters):
    """Checks the dropout setup (positive and smaller than 1.0)"""
    dropout = parameters.get("dropout")
    if dropout is not None and not 0.0 <= dropout <= 1.0:
        raise CheckError(f"Dropout must be between 0 and 1, got: {dropout}")

    layer_type = parameters.get("layer_type")
    if dropout is not None and dropout > 0.0 and layer_type == "dense_per_flavour":
        raise CheckError(
            "Dropout is not compatible with the dense_per_flavour layer type, "
            "please use instead `parameters::layer_type: dense`"
        )


def check_tensorboard(tensorboard):
    """Check that the tensorbard callback can be enabled correctly"""
    if tensorboard is not None:
        # Check that Tensorflow is installed
        try:
            import tensorflow
        except ModuleNotFoundError as e:
            raise ModuleNotFoundError(
                "The tensorboard callback requires `tensorflow` to be installed"
            ) from e

        weight_freq = tensorboard.get("weight_freq", 0)
        if weight_freq < 0:
            raise CheckError(
                f"The frequency at which weights are saved must be greater than 0, received {weight_freq}"
            )


def check_lagrange_multipliers(parameters, key):
    """Checks the parameters in a lagrange multiplier dictionary
    are correct, e.g. for positivity and integrability"""
    lagrange_dict = parameters.get(key)
    if lagrange_dict is None:
        return
    multiplier = lagrange_dict.get("multiplier")
    if multiplier is not None and multiplier < 0:
        log.warning("The %s multiplier is below 0, it will produce a negative loss", key)
    elif multiplier is not None and multiplier == 0:
        log.warning("The %s multiplier is 0 so it won't contribute to the loss", key)
    threshold = lagrange_dict.get("threshold")
    if threshold is not None and not _is_floatable(threshold):
        raise CheckError(f"The {key}::threshold must be a number, received: {threshold}")


def check_model_file(save, load):
    """Checks whether the model_files given in the runcard are acceptable"""
    if save:
        if not isinstance(save, str):
            raise CheckError(f"Model file to save to: {save} not understood")
        # Since the file to save to will be found inside the replica folder, it should writable as all the others

    if load:
        if not isinstance(load, str):
            raise CheckError(f"Model file to load: {load} not understood, str expected")
        if not os.path.isfile(load):
            raise CheckError(f"Model file to load: {load} can not be opened, does it exist?")
        if not os.access(load, os.R_OK):
            raise CheckError(f"Model file to load: {load} cannot be read, permission denied")
        if os.stat(load).st_size == 0:
            raise CheckError(f"Model file {load} seems to be empty")


@make_argcheck
def wrapper_check_NN(tensorboard, save, load, parameters):
    """Wrapper function for all NN-related checks"""
    check_tensorboard(tensorboard)
    check_model_file(save, load)
    check_existing_parameters(parameters)
    check_consistent_layers(parameters)
    check_stopping(parameters)
    check_layer_type_implemented(parameters)
    check_dropout(parameters)
    check_lagrange_multipliers(parameters, "integrability")
    check_lagrange_multipliers(parameters, "positivity")
    # Checks that need to import the backend (and thus take longer) should be done last
    check_optimizer(parameters["optimizer"])
    check_initializer(parameters["initializer"])


def check_hyperopt_architecture(architecture):
    """Checks whether the scanning setup for the NN architecture works
    - Initializers are valid
    - Dropout setup is valid
    - No 'min' is greater than its corresponding 'max'
    """
    if architecture is None:
        return
    initializers = architecture.get("initializers")
    if initializers is not None:
        for init in initializers:
            check_initializer(init)
    # Check that max-min are correct
    dropout = architecture.get("max_drop")
    if dropout is not None and not 0.0 <= dropout <= 1.0:
        raise CheckError(f"max_drop must be between 0 and 1, got: {dropout}")
    min_u = architecture.get("min_units", 1)
    # Set a minimum number of units in case none is defined to check later if the maximum is sane
    if min_u <= 0:
        raise CheckError(f"All layers must have at least 1 unit, got min_units: {min_u}")
    max_u = architecture.get("max_units")
    if max_u is not None and max_u < min_u:
        raise CheckError(
            "The maximum number of units must be bigger than the minimum "
            f" but got min: {min_u}, max: {max_u}"
        )


def check_hyperopt_positivity(positivity_dict):
    """Checks that the positivity multiplier and initial values are sensible and valid"""
    if positivity_dict is None:
        return
    min_mul = positivity_dict.get("min_multiplier")
    max_mul = positivity_dict.get("max_multiplier")
    if max_mul is not None or min_mul is not None:
        if max_mul is None:
            raise CheckError("Need to set a maximum positivity multiplier if the minimum is set")
        if min_mul is not None and max_mul <= min_mul:
            raise CheckError("The minimum multiplier cannot be greater than the maximum")
    min_ini = positivity_dict.get("min_initial")
    max_ini = positivity_dict.get("max_initial")
    if max_ini is not None or min_ini is not None:
        if max_ini is None or min_ini is None:
            raise CheckError(
                "Need to set both the max_initial and the min_initial positivity values"
            )
        if max_ini <= min_ini:
            raise CheckError("The minimum initial value cannot be greater than the maximum")


def check_kfold_options(kfold):
    """Warns the user about potential bugs on the kfold setup"""
    threshold = kfold.get("threshold")
    if threshold is not None and threshold < 2.0:
        log.warning("The kfolding loss threshold might be too low: %f", threshold)
    penalty_selection = kfold.get("penalties", [])
    for penalty in penalty_selection:
        if not hasattr(penalties_module, penalty):
            raise CheckError(
                f"The penalty '{penalty}' is not recognized, ensure it is implemented in hyper_optimization/penalties.py"
            )

    loss_type = kfold.get("loss_type")
    if loss_type is not None:
        if loss_type not in IMPLEMENTED_LOSSES:
            raise CheckError(
                f"Loss type '{loss_type}' is not recognized, "
                "ensure it is implemented in the HyperLoss class in hyper_optimization/rewards.py."
                "Options so far are 'chi2' or 'phi2'."
            )
    replica_statistic = kfold.get("replica_statistic")
    if replica_statistic is not None:
        if replica_statistic not in IMPLEMENTED_STATS:
            raise CheckError(
                f"The replica statistic '{replica_statistic}' is not recognized, "
                "ensure it is implemented in the HyperLoss class in hyper_optimization/rewards.py"
            )
    fold_statistic = kfold.get("fold_statistic")
    if fold_statistic is not None:
        if fold_statistic not in IMPLEMENTED_STATS:
            raise CheckError(
                f"The fold statistic '{fold_statistic}' is not recognized, "
                "ensure it is implemented in the HyperLoss class in hyper_optimization/rewards.py"
            )

    partitions = kfold["partitions"]
    # Check specific errors for specific targets
    loss_target = kfold.get("fold_statistic")  # TODO: haven't updated this
    if loss_target == "fit_future_tests":
        if len(partitions) == 1:
            raise CheckError("Cannot use target 'fit_future_tests' with just one partition")
        if partitions[-1]["datasets"]:
            log.warning("Last partition in future test is not empty, some datasets will be ignored")


def check_correct_partitions(kfold, data):
    """Ensures that all experimennts in all partitions
    are included in  the fit definition"""
    # Get all datasets
    datasets = list(map(str, data))
    for partition in kfold["partitions"]:
        fold_sets = partition["datasets"]
        for dset in fold_sets:
            if dset not in datasets:
                raise CheckError(f"The k-fold defined dataset {dset} is not part of the fit")


def check_hyperopt_stopping(stopping_dict):
    """Checks that the options selected for the stopping are consistent"""
    if stopping_dict is None:
        return
    min_ep = stopping_dict.get("min_epochs")
    max_ep = stopping_dict.get("max_epochs")
    if max_ep is not None or min_ep is not None:
        if min_ep is None or max_ep is None:
            raise CheckError("Need to set both the max_epochs and the min_epochs")
        if min_ep < 1:
            raise CheckError(f"Can't run for less than 1 epoch: selected min_ep = {min_ep}")
        if max_ep <= min_ep:
            raise CheckError(f"min_epochs cannot be greater than max_epochs: ({min_ep} > {max_ep})")
    min_pat = stopping_dict.get("min_patience")
    max_pat = stopping_dict.get("max_patience")
    if min_pat is not None or max_pat is not None:
        if min_pat is not None and min_pat < 0.0:
            raise CheckError(f"min_patience cannot be less than 0.0: selected min_pat = {min_pat}")
        if max_pat is not None:
            if max_pat > 1.0:
                raise CheckError(
                    f"max_patience cannot be greater than 1.0: selected max_pat = {max_pat}"
                )
            if min_pat is not None and max_pat < min_pat:
                raise CheckError(
                    f"min_patience cannot be greater than max_patience: ({min_pat} > {max_pat})"
                )


@make_argcheck
def wrapper_hyperopt(hyperopt, hyperscan_config, kfold, data):
    """Wrapper function for all hyperopt-related checks
    No check is performed if hyperopt is not active
    """
    if not hyperopt:
        return
    if hyperscan_config is None:
        raise CheckError("Can't perform hyperoptimization without the hyperscan_config key")
    if kfold is None:
        raise CheckError("Can't perform hyperoptimization without folds")
    check_hyperopt_stopping(hyperscan_config.get("stopping"))
    check_hyperopt_architecture(hyperscan_config.get("architecture"))
    check_hyperopt_positivity(hyperscan_config.get("positivity"))
    check_kfold_options(kfold)
    check_correct_partitions(kfold, data)


def check_sumrules(sum_rules):
    """Checks that the chosen option for the sum rules are sensible"""
    if isinstance(sum_rules, bool):
        return
    accepted_options = ["ALL", "MSR", "VSR", "TSR", "ALLBUTCSR"]
    if sum_rules.upper() in accepted_options:
        return
    raise CheckError(f"The only accepted options for the sum rules are: {accepted_options}")


# Checks on the physics
@make_argcheck
def check_consistent_basis(sum_rules, fitbasis, basis, theoryid, parameters):
    """Checks the fitbasis setup for inconsistencies
    - Checks the sum rules can be imposed
    - Correct flavours for the selected basis
    - Correct ranges (min < max) for the small and large-x exponents
    - When feature scaling is active, the large_x interpolation is not set
    """
    check_sumrules(sum_rules)
    # Check that there are no duplicate flavours and that parameters are sane
    flavs = []
    for flavour_dict in basis:
        name = flavour_dict["fl"]
        smallx = flavour_dict["smallx"]
        if smallx[0] > smallx[1]:
            raise CheckError(f"Wrong smallx range for flavour {name}: {smallx}")
        if name in flavs:
            raise CheckError(f"Repeated flavour name: {name}. Check basis dictionary")
        flavs.append(name)

        # Large-x is allowed to not exist if feature scaling is enabled
        if parameters.get("feature_scaling_points") is not None:
            if "largex" in flavour_dict and not flavour_dict["largex"] == [0.0, 0.0]:
                raise CheckError("No largex exponent allowed when feature_scaling_points is set")
        else:
            largex = flavour_dict["largex"]
            if largex[0] > largex[1]:
                raise CheckError(f"Wrong largex range for flavour {name}: {largex}")

    # Finally check whether the basis considers or not charm
    # Check that the basis given in the runcard is one of those defined in validphys.pdfbases
    vp_basis = check_basis(fitbasis, flavs)["basis"]
    # Now check that basis and theory id are consistent
    has_c = vp_basis.has_element("c") or vp_basis.has_element("T15") or vp_basis.has_element("cp")
    if theoryid.get_description()["IC"] and not has_c:
        raise CheckError(f"{theoryid} (intrinsic charm) is incompatible with basis {fitbasis}")
    if not theoryid.get_description()["IC"] and has_c:
        raise CheckError(f"{theoryid} (perturbative charm) is incompatible with basis {fitbasis}")

    check_basis_with_layers(basis, vp_basis, parameters)


@make_argcheck
def check_consistent_parallel(parameters, parallel_models):
    """Checks whether the multiple-replica fit options are consistent among them
    i.e., that the trvl seed is fixed and the layer type is correct
    """
    if not parallel_models:
        return
    if parameters.get("layer_type") not in ("dense"):
        raise CheckError(
            "Parallelization has only been tested with layer_type=='dense', set `parallel_models: false`"
        )


@make_argcheck
def check_deprecated_options(fitting, parameters):
    """Checks whether the runcard is using deprecated options"""
    options_outside = ["trvlseed", "nnseed", "mcseed", "save", "load", "genrep", "parameters"]
    for option in options_outside:
        if option in fitting:
            raise CheckError(
                f"The key '{option}' should be top-level key and not part of the 'fitting' namespace"
            )
    if "epochs" in fitting:
        raise CheckError("The key 'epoch' should only appear as part of the 'parameters' namespace")
    nnfit_options = ["seed", "rnalgo", "fitmethod", "nmutants", "paramtype", "nnodes"]
    for option in nnfit_options:
        if option in fitting:
            log.warning("'fitting::%s' is an nnfit-only key, it will be ignored", option)
    if "interpolation_points" in parameters:
        raise CheckError(
            "`interpolation_points` no longer accepted, please change to `feature_scaling_points`"
        )


@make_argcheck
def check_multireplica_qed(replicas, fiatlux):
    if fiatlux is not None:
        if len(replicas) > 1:
            raise CheckError("At the moment, running a multireplica QED fits is not allowed.")
        
@make_argcheck
def check_photonQED_exists(theoryid, fiatlux):
    """Check that the Photon QED set for this theoryid and luxset exists"""
    if fiatlux is not None:
      luxset = fiatlux['luxset']
      try:
          _ = FallbackLoader().check_photonQED(theoryid.id, luxset)
          log.info(f"Photon QED set found for {theoryid.id} with luxset {luxset}.")
      except FileNotFoundError:
          log.warning(f"No Photon QED set found for {theoryid} with luxset {luxset} and "\
                      "will be compute using FiatLux. This may impact performance.")


@make_argcheck
def check_polarized_configs(fitting, fitbasis, positivity_bound):
    if fitbasis.startswith("POLARIZED_"):
        if positivity_bound is None:
            raise CheckError(
                "For polarized fits, the 'positivity_bound' key has to be defined in the runcard."
            )
        if positivity_bound.get("unpolarized_bc") is None:
            raise CheckError(
                "For polarized fits, the unpolarized PDF has to be defined in positivity_bound::unpolarized_bc."
            )
        if fitting.get("sum_rules", True) and fitting.get("sum_rules") != "TSR":
            raise CheckError("The 'sum_rules' key needs to be 'TSR' for polarised PDF fits.")


@make_argcheck
def check_eko_exists(theoryid):
    """Check that an eko for this theory exists.
    Since there might still be theories without an associated eko,
    this function raises a logger' error instead of an Exception."""
    try:
        _ = FallbackLoader().check_eko(theoryid.id)
    except FileNotFoundError:
        log.error(f"No eko found for {theoryid}")
