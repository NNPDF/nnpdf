""" 
    Test for the model generation

    These tests check that the generated NN are as expected
    It checks that both the number of layers and the shape
    of the weights of the layers are what is expected
"""
import numpy as np
import n3fit.model_gen
from n3fit.backends import MetaModel
from n3fit.backends import operations as op

INSIZE = 16
OUT_SIZES = (4, 3)
BASIS_SIZE = 3


def test_generate_dense_network():
    nodes_in = INSIZE
    nodes_out = OUT_SIZES
    activations = ["sigmoid", "tanh"]
    layers = n3fit.model_gen.generate_dense_network(nodes_in, nodes_out, activations)
    arr = np.random.rand(1, INSIZE)
    input_layer = op.numpy_to_input(arr)
    curr_layer = input_layer
    for layer in layers:
        curr_layer = layer(curr_layer)
    modelito = MetaModel({"input": input_layer}, curr_layer)
    # The number of layers should be input layer + len(OUT_SIZES)
    assert len(modelito.layers) == len(OUT_SIZES) + 1
    # Check that the number of parameters is as expected
    # We expect 4 weights where the two first ones are
    # (INSIZE, OUT_SIZE[0]) (OUT_SIZE[0],)
    # and the second one
    # (OUT_SIZE[0], OUT_SIZE[1]) (OUT_SIZE[1],)
    expected_sizes = [
        (INSIZE, OUT_SIZES[0]),
        (OUT_SIZES[0],),
        OUT_SIZES,
        (OUT_SIZES[1],),
    ]
    for weight, esize in zip(modelito.weights, expected_sizes):
        assert weight.shape == esize


def test_generate_dense_per_flavour_network():
    nodes_in = INSIZE
    nodes_out = OUT_SIZES
    activations = ["sigmoid", "tanh"]
    layers = n3fit.model_gen.generate_dense_per_flavour_network(
        nodes_in, nodes_out, activations, basis_size=BASIS_SIZE
    )
    arr = np.random.rand(1, INSIZE)
    input_layer = op.numpy_to_input(arr)
    curr_layer = input_layer
    for layer in layers:
        curr_layer = layer(curr_layer)
    modelito = MetaModel({"input": input_layer}, curr_layer)
    # The number of layers should be input + BASIS_SIZE*len(OUT_SIZES) + concatenate
    assert len(modelito.layers) == BASIS_SIZE * len(OUT_SIZES) + 2
    # The shape for this network of denses for flavours will depend on the basis_size
    expected_sizes = []
    expected_sizes += BASIS_SIZE * [(INSIZE, OUT_SIZES[0]), (OUT_SIZES[0],)]
    expected_sizes += BASIS_SIZE * [(OUT_SIZES[0], 1), (1,)]
    for weight, esize in zip(modelito.weights, expected_sizes):
        assert weight.shape == esize
