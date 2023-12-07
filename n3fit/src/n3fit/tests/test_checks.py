"""
    Test for the n3fit checks
"""
import pytest

from n3fit import checks
from n3fit.backends import MetaModel
from reportengine.checks import CheckError


def test_existing_parameters():
    """Test the parameters check"""
    needed_params = checks.NN_PARAMETERS
    keys_in = dict.fromkeys(needed_params, "Something")
    # Check that with all parameters pass
    checks.check_existing_parameters(keys_in)
    # Check that if something is missing, it fails
    with pytest.raises(CheckError):
        keys_in.pop(needed_params[0])
        checks.check_existing_parameters(keys_in)


def test_consistent_layers():
    """Test the check fails for inconsistent layers-activation groups"""
    params = {"nodes_per_layer": [3, 4], "activation_per_layer": ["sigmoid"]}
    with pytest.raises(CheckError):
        checks.check_consistent_layers(params)


def test_check_stopping():
    """Test the stopping related options correctly fail"""
    params = {"stopping_patience": 1.5}
    with pytest.raises(CheckError):
        checks.check_stopping(params)
    params["stopping_patience"] = 0.5
    params["epochs"] = -4
    with pytest.raises(CheckError):
        checks.check_stopping(params)
    params["epochs"] = 4
    checks.check_stopping(params)


def test_check_basis_with_layer():
    """Test that it fails with layers that do not match flavours"""
    flavs = ["g", "ubar"]
    layers = [4, 5, 9]
    with pytest.raises(CheckError):
        checks.check_basis_with_layers({"basis": flavs}, {"nodes_per_layer": layers})


def test_check_optimizer():
    """Test that the optimizer check is correct"""
    params = {"optimizer_name": "fake_non_existing"}
    with pytest.raises(CheckError):
        checks.check_optimizer(params)
    params["optimizer_name"] = "RMSprop"
    params["wrong_parameters"] = 7
    with pytest.raises(CheckError):
        checks.check_optimizer(params)


def test_check_initializer():
    """Test that it fails with a wrong initializer"""
    with pytest.raises(CheckError):
        checks.check_initializer("Wrong_one")


def test_check_dropout():
    """Test the dropout checks"""
    with pytest.raises(CheckError):
        checks.check_dropout({"dropout": 1.5})
    with pytest.raises(CheckError):
        checks.check_dropout({"dropout": -0.5})
    checks.check_dropout({"dropout": 0.5})


def test_check_hyperopt_architecture():
    """Test the checks for the hyperopt architecture"""
    params = {"initializers": ["Fake_bad_non"]}

    def autocheck():
        with pytest.raises(CheckError):
            checks.check_hyperopt_architecture(params)

    # Check that it fails on wrong initializer
    autocheck()
    params["initializers"] = None
    # Check it fails on wrong dropouts
    params["max_drop"] = 1.5
    autocheck()  # dropout > 1.0
    params["max_drop"] = -0.5
    autocheck()  # dropout < 0.0
    params["max_drop"] = 0.5
    # Check that it fails in wrong number of units
    params["min_units"] = 0
    autocheck()
    params["min_units"] = -4
    autocheck()  # units <= 0
    params["min_units"] = 5
    params["max_units"] = 4
    autocheck()  # min > max


def test_check_hyperopt_positivity():
    """Test the positivity settings in hyperopt"""
    params = {"min_multiplier": 1.0}

    def autocheck():
        with pytest.raises(CheckError):
            checks.check_hyperopt_positivity(params)

    autocheck()  # missing max
    params["max_multiplier"] = 0.4
    autocheck()  # max < min
    params["max_multiplier"] = 2.0
    params["min_initial"] = 0.4
    autocheck()  # missing max initial
    params["min_initial"] = None
    params["max_initial"] = 0.4
    autocheck()  # missing min
    params["min_initial"] = 1.4
    autocheck()  # min > max


def test_check_kfold_options():
    """Test the kfold option (not including datasets)"""
    params = {"penalties": ["Fake_penalty_doesnt_exists"]}
    with pytest.raises(CheckError):
        checks.check_kfold_options(params)


def test_check_hyperopt_stopping():
    """Test stopping hyperscan options"""
    params = {"max_epochs": 4}

    def autocheck():
        with pytest.raises(CheckError):
            checks.check_hyperopt_stopping(params)

    autocheck()  # min missing
    params["max_epochs"] = None
    params["min_epochs"] = -4
    autocheck()  # max missing
    params["max_epochs"] = 10
    autocheck()  # negative min
    params["min_epochs"] = 15
    autocheck()  # min > max
    params["max_epochs"] = 1000
    params["min_patience"] = -10
    autocheck()  # negative min patience
    params["min_patience"] = 0.5
    params["max_patience"] = 0.3
    autocheck()  # max > min
    params["max_patience"] = 1.3
    autocheck()  # max > 1.0
