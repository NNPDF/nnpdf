"""
This module contains checks to be perform by n3fit on the input
"""
import logging
from reportengine.checks import make_argcheck, CheckError
from validphys.pdfbases import check_basis

log = logging.getLogger(__name__)

NN_PARAMETERS = ["nodes_per_layer", "optimizer", "activation_per_layer"]

# Checks on the NN parameters
def check_existing_parameters(parameters):
    """ Check that non-optional parameters are defined and are not empty """
    for param_name in NN_PARAMETERS:
        if param_name in parameters:
            val = parameters[param_name]
            if len(val) == 0:
                raise CheckError(f"The parameter {param_name} cannot be empty")
        else:
            raise CheckError(f"Missing {param_name} parameter in the runcard")


def check_consistent_layers(parameters):
    """ Checks that all layers have an activation function defined """
    npl = len(parameters["nodes_per_layer"])
    apl = len(parameters["activation_per_layer"])
    if npl != apl:
        raise CheckError(f"Number of layers ({npl}) does not match activation functions: {apl}")


def check_stopping(parameters):
    """ Checks whether the stopping-related options are sane:
    stopping patience as a ratio between 0 and 1
    and positive number of epochs
    """
    spt = parameters.get("stopping_patience")
    if spt is not None and not 0.0 <= spt <= 1.0:
        raise CheckError(f"The stopping_patience must be between 0 and 1, got: {spt}")
    epochs = parameters["epochs"]
    if epochs < 1:
        raise CheckError(f"Needs to run at least 1 epoch, got: {epochs}")


def check_basis_with_layers(fitting, parameters):
    """ Check that the last layer matches the number of flavours defined in the runcard"""
    number_of_flavours = len(fitting["basis"])
    last_layer = parameters["nodes_per_layer"][-1]
    if number_of_flavours != last_layer:
        raise CheckError(
            f"The number of nodes in the last layer ({last_layer}) does not"
            " match the number of flavours: ({number_of_flavours})"
        )


def check_optimizer(optimizer_dict):
    """ Checks whether the optimizer setup is valid """
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
    """ Checks whether the initializer is implemented """
    from n3fit.backends import MetaLayer

    accepted_init = MetaLayer.initializers
    if initializer not in accepted_init:
        raise CheckError(f"Initializer {initializer} not accepted by {MetaLayer}")


def check_dropout(parameters):
    """ Checks the dropout setup (positive and smaller than 1.0) """
    dropout = parameters.get("dropout")
    if dropout is not None and not 0.0 <= dropout <= 1.0:
        raise CheckError(f"Dropout must be between 0 and 1, got: {dropout}")


@make_argcheck
def wrapper_check_NN(fitting):
    """ Wrapper function for all NN-related checks """
    parameters = fitting["parameters"]
    check_existing_parameters(parameters)
    check_consistent_layers(parameters)
    check_basis_with_layers(fitting, parameters)
    check_stopping(parameters)
    check_dropout(parameters)
    # Checks that need to import the backend (and thus take longer) should be done last
    #     check_optimizer(parameters["optimizer"]) # this check is waiting for PR 783
    check_initializer(parameters["initializer"])


def check_hyperopt_architecture(architecture):
    """ Checks whether the scanning setup for the NN architecture works
    - Initializers are valid
    - Droput setup is valid
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
            "The maximum number of units must be bigger than the minimum"
            f" but got min: {min_u}, max: {max_u}"
        )


def check_hyperopt_positivity(positivity_dict):
    """ Checks that the positivity multiplier and initial values are sensible and valid
    """
    if positivity_dict is None:
        return
    min_mul = positivity_dict.get("min_multiplier")
    max_mul = positivity_dict.get("max_multiplier")
    if max_mul is not None or min_mul is not None:
        if max_mul is None:
            raise CheckError("Need to set a maximum positivity multiplier is the minimum is set")
        if min_mul is not None and max_mul <= min_mul:
            raise CheckError("The minimum multiplier cannot be greater than the maximum")
    min_ini = positivity_dict.get("min_initial")
    max_ini = positivity_dict.get("max_initial")
    if max_ini is not None or min_ini is not None:
        if max_ini is None or min_ini is None:
            raise CheckError(
                "Need to set both the maximum and minimum positivitiy initial value if any of the two is set"
            )
        if min_ini is not None and max_ini <= min_mul:
            raise CheckError("The minimum initial value cannot be greater than the maximum")


def check_kfold_options(kfold):
    """ Warns the user about potential bugs on the kfold setup"""
    threshold = kfold.get("threshold")
    if threshold is not None and threshold < 2.0:
        log.warning("The kfolding loss threshold might be too low: %f", threshold)


def check_correct_partitions(kfold, experiments):
    """ Ensures that all experimennts in all partitions
    are included in  the fit definition """
    # Get all datasets
    datasets = []
    for exp in experiments:
        datasets += [i.name for i in exp.datasets]
    for partition in kfold["partitions"]:
        fold_sets = partition["datasets"]
        for dset in fold_sets:
            if dset not in datasets:
                raise CheckError(f"The k-fold defined dataset {dset} is not part of the fit")


@make_argcheck
def wrapper_hyperopt(hyperopt, hyperscan, fitting, experiments):
    """ Wrapper function for all hyperopt-related checks
    No check is performed if hyperopt is not active
    """
    if not hyperopt:
        return None
    if fitting["genrep"]:
        raise CheckError("Generation of replicas is not accepted during hyperoptimization")
    if hyperscan is None:
        raise CheckError("Can't perform hyperoptimization without the hyperscan key")
    if "kfold" not in hyperscan:
        raise CheckError("The hyperscan::kfold dictionary is not defined")
    check_hyperopt_architecture(hyperscan.get("architecture"))
    check_hyperopt_positivity(hyperscan.get("positivity"))
    check_kfold_options(hyperscan["kfold"])
    check_correct_partitions(hyperscan["kfold"], experiments)


# Checks on the physics
@make_argcheck
def check_consistent_basis(fitting):
    """ Checks the fitbasis setup for inconsistencies
    - Correct flavours for the selected basis
    - Correct ranges (min < max) for the small and large-x exponents
    """
    fitbasis = fitting["fitbasis"]
    # Check that there are no duplicate flavours and that parameters are sane
    flavs = []
    for flavour_dict in fitting["basis"]:
        name = flavour_dict["fl"]
        smallx = flavour_dict["smallx"]
        if smallx[0] > smallx[1]:
            raise CheckError(f"Wrong smallx range for flavour {name}: {smallx}")
        largex = flavour_dict["largex"]
        if largex[0] > largex[1]:
            raise CheckError(f"Wrong largex range for flavour {name}: {largex}")
        if name in flavs:
            raise CheckError(f"Repeated flavour name: {name}. Check basis dictionary")
        flavs.append(name)
    # Check that the basis given in the runcard is one of those defined in validphys.pdfbases
    check_basis(fitbasis, flavs)
