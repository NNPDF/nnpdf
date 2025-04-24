"""
Test for the model generation

These tests check that the generated NN are as expected
It checks that both the number of layers and the shape
of the weights of the layers are what is expected
"""

from n3fit.backends import NN_PREFIX, Input
from n3fit.model_gen import _generate_nn

INSIZE = 16
OUT_SIZES = (4, 3)
BASIS_SIZE = OUT_SIZES[-1]


def _common_generation(architecture, dropout=0.0):
    """Generate a NN with shared configuration with only
    some free parameters"""
    config = {
        "architecture": architecture,
        "nodes": OUT_SIZES,
        "activations": ["sigmoid", "tanh"],
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
    expected_sizes = [(INSIZE, OUT_SIZES[0]), (OUT_SIZES[0]), (*OUT_SIZES,), (OUT_SIZES[1])]
    for weight, esize in zip(nn.weights, expected_sizes):
        assert weight.shape == esize


def test_generate_dense_per_flavour_network():
    nn = _common_generation("dense_per_flavour")

    # The number of layers should be input + BASIS_SIZE*len(OUT_SIZES) + concatenate
    assert len(nn.layers) == BASIS_SIZE * len(OUT_SIZES) + 2
    # The shape for this network of denses for flavours will depend on the basis_size
    expected_sizes = []
    expected_sizes += BASIS_SIZE * [(INSIZE, OUT_SIZES[0]), (OUT_SIZES[0],)]
    expected_sizes += BASIS_SIZE * [(OUT_SIZES[0], 1), (1,)]
    for weight, esize in zip(nn.weights, expected_sizes):
        assert weight.shape == esize
