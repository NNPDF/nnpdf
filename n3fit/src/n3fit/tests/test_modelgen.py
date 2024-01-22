""" 
    Test for the model generation

    These tests check that the generated NN are as expected
    It checks that both the number of layers and the shape
    of the weights of the layers are what is expected
"""
from n3fit.backends.keras.metamodel import NN_PREFIX
from n3fit.model_gen import generate_nn

INSIZE = 16
OUT_SIZES = (4, 3)
BASIS_SIZE = 3

COMMON_ARGS = {
    "nodes_in": INSIZE,
    "nodes": OUT_SIZES,
    "activations": ["sigmoid", "tanh"],
    "initializer_name": "glorot_uniform",
    "replica_seeds": [0],
    "dropout": 0.0,
    "regularizer": None,
    "regularizer_args": {},
    "last_layer_nodes": BASIS_SIZE,
}


def test_generate_dense_network():
    nn = generate_nn("dense", **COMMON_ARGS).get_layer(f"{NN_PREFIX}_0")

    # The number of layers should be input layer + len(OUT_SIZES)
    assert len(nn.layers) == len(OUT_SIZES) + 1
    # Check that the number of parameters is as expected
    # We expect 4 weights where the two first ones are
    # (INSIZE, OUT_SIZE[0]) (OUT_SIZE[0],)
    # and the second one
    # (OUT_SIZE[0], OUT_SIZE[1]) (OUT_SIZE[1],)
    expected_sizes = [(INSIZE, OUT_SIZES[0]), (OUT_SIZES[0],), OUT_SIZES, (OUT_SIZES[1],)]
    for weight, esize in zip(nn.weights, expected_sizes):
        assert weight.shape == esize


def test_generate_dense_per_flavour_network():
    nn = generate_nn("dense_per_flavour", **COMMON_ARGS).get_layer(f"{NN_PREFIX}_0")

    # The number of layers should be input + BASIS_SIZE*len(OUT_SIZES) + concatenate
    assert len(nn.layers) == BASIS_SIZE * len(OUT_SIZES) + 2
    # The shape for this network of denses for flavours will depend on the basis_size
    expected_sizes = []
    expected_sizes += BASIS_SIZE * [(INSIZE, OUT_SIZES[0]), (OUT_SIZES[0],)]
    expected_sizes += BASIS_SIZE * [(OUT_SIZES[0], 1), (1,)]
    for weight, esize in zip(nn.weights, expected_sizes):
        assert weight.shape == esize
