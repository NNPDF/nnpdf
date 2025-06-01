"""
Test for the model generation

These tests check that the generated NN are as expected
It checks that both the number of layers and the shape
of the weights of the layers are what is expected
"""

from dataclasses import asdict

import pytest

from n3fit.backends import Input
from n3fit.model_gen import _generate_nn, _ReplicaSettings
from nnpdf_data.utils import parse_input

INSIZE = 16
OUT_SIZES = (4, 7, 3)
BASIS_SIZE = OUT_SIZES[-1]


def _common_generation(architecture, dropout=0.0):
    """Generate a NN with shared configuration with only
    some free parameters"""
    config = {
        "architecture": architecture,
        "nodes": OUT_SIZES,
        "activations": ["sigmoid", "sigmoid", "tanh"],
        "initializer": "glorot_uniform",
        "seed": 27,
        "dropout_rate": dropout,
        "regularizer": None,
        "regularizer_args": {},
    }
    xin = Input(shape=(None, INSIZE), batch_size=1)
    return _generate_nn(xin, **config)


def test_generate_dense_network():
    nn_w_dropout = _common_generation("dense", dropout=0.4)
    nn = _common_generation("dense")

    # The number of layers should be input layer + len(OUT_SIZES)
    assert len(nn.layers) == len(OUT_SIZES) + 1
    # And one more with dropout
    assert len(nn_w_dropout.layers) == len(OUT_SIZES) + 2

    # Check that the number of parameters is as expected
    expected_sizes = [(INSIZE, OUT_SIZES[0])]
    for i, oz in enumerate(OUT_SIZES[:-1]):
        expected_sizes.append((oz,))
        expected_sizes.append((oz, OUT_SIZES[i + 1]))
    expected_sizes.append((OUT_SIZES[-1],))
    for weight, esize in zip(nn.weights, expected_sizes):
        assert weight.shape == esize


def test_generate_dense_per_flavour_network():
    nn = _common_generation("dense_per_flavour")

    # The number of layers should be input + BASIS_SIZE*len(OUT_SIZES) + concatenate
    assert len(nn.layers) == BASIS_SIZE * len(OUT_SIZES) + 2
    # The shape for this network of denses for flavours will depend on the basis_size
    expected_sizes = []
    expected_sizes += BASIS_SIZE * [(INSIZE, OUT_SIZES[0]), (OUT_SIZES[0],)]
    for i, oz in enumerate(OUT_SIZES[1:-1]):
        expected_sizes += BASIS_SIZE * [(OUT_SIZES[i], oz), (oz,)]
    expected_sizes += BASIS_SIZE * [(oz, 1), (1,)]

    for weight, esize in zip(nn.weights, expected_sizes):
        assert weight.shape == esize


def test_replica_settings():
    """Checks that the _ReplicaSettings object works as expected and
    that it matches the input of _generate_nn"""
    config = {
        "seed": 8,
        "nodes": [4, 10],
        "activations": ["linear"] * 2,
        "architecture": "dense",
        "initializer": "glorot_uniform",
        "dropout_rate": 0.4,
    }

    rsettings = parse_input(config, _ReplicaSettings)

    with pytest.raises(ValueError):
        ctmp = {**config, "regularizer_args": {"some": 4}}
        _ReplicaSettings(**ctmp)

    with pytest.raises(ValueError):
        ctmp = {**config, "nodes": [2]}
        _ReplicaSettings(**ctmp)

    x = Input(shape=(None, 2))
    nn = _generate_nn(x, 0, **asdict(rsettings))

    nn.layers == (1 + len(rsettings.nodes) + 1 * rsettings.dropout_rate)
