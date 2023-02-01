"""
    Test the n3fit-validphys interface
"""

import numpy as np
from hypothesis import given, settings, example
from hypothesis.strategies import integers
from validphys.pdfgrids import xplotting_grid, distance_grids
from n3fit.vpinterface import N3PDF, integrability_numbers, compute_arclength
from n3fit.model_gen import pdfNN_layer_generator


def generate_n3pdf(layers=1, members=1, name="n3fit"):
    """Generate a N3PDF model"""
    fake_fl = [
        {"fl": i, "largex": [0, 1], "smallx": [1, 2]}
        for i in ["u", "ubar", "d", "dbar", "c", "g", "s", "sbar"]
    ]
    nodes = list(np.random.randint(1, 10, size=layers)) + [8]
    activations = ["tanh"] * layers + ["linear"]
    pdf_model = pdfNN_layer_generator(
        nodes=nodes,
        activations=activations,
        seed=np.random.randint(100),
        flav_info=fake_fl,
        parallel_models=members,
        fitbasis="FLAVOUR"
    )
    return N3PDF(pdf_model, name=name)


@given(integers(1, 3), integers(0, 3))
@example(1, 3)
@example(3, 0)
@settings(deadline=12000, max_examples=5)
def test_N3PDF(members, layers):
    """Test the N3PDF class produces the right members and atributes
    for several different combinations of the number of layers and members
    """
    xsize = np.random.randint(2, 20)
    xx = np.random.rand(xsize)
    n3pdf = generate_n3pdf(layers=layers, members=members)
    assert len(n3pdf) == members
    w = n3pdf.get_nn_weights()
    assert len(w) == members
    assert len(w[0]) == 16 + (layers + 1) * 2  # 16=8*2 preprocessing
    ret = n3pdf(xx)
    assert ret.shape == (members, xsize, 14)
    int_numbers = integrability_numbers(n3pdf)
    if members == 1:
        assert int_numbers.shape == (5,)
    else:
        assert int_numbers.shape == (members, 5)
    assert compute_arclength(n3pdf).shape == (5,)
    # Try to get a plotting grid
    res = xplotting_grid(n3pdf, 1.6, xx)
    assert res.grid_values.data.shape == (members, 8, xsize)


def test_vpinterface():
    """Test several uses of the n3fit - VP interface"""
    fit_1 = generate_n3pdf(layers=2, members=5, name="fit_1")
    fit_2 = generate_n3pdf(layers=4, members=3, name="fit_2")
    xgrid = np.concatenate([np.logspace(-6, -1, 20), np.linspace(0.15, 0.9, 20)])
    res_1 = xplotting_grid(fit_1, 1.6, xgrid)
    res_2 = xplotting_grid(fit_2, 1.6, xgrid)
    distances = distance_grids([fit_1, fit_2], [res_1, res_2], 0)
    assert len(distances) == 2
    assert distances[0].grid_values.data.shape == (1, 8, 40)
    assert distances[1].grid_values.data.shape == (1, 8, 40)
    np.testing.assert_allclose(distances[0].grid_values.data, 0.0)
    assert not np.allclose(distances[1].grid_values.data, 0.0)
