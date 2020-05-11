"""
    Tests for the layers of n3fit
    This module checks that the layers do what they would do with numpy
"""

import numpy as np
from n3fit.backends import operations as op
import n3fit.layers as layers
from validphys.pdfbases import rotation


FLAVS = 3
XSIZE = 4
NDATA = 3
THRESHOLD = 1e-6

# Helper functions
def generate_input_had(flavs=3, xsize=2, ndata=4, n_combinations=None):
    """ Generates fake input (fktable and array of combinations) for the hadronic convolution

    Parameters
    ----------
    flavs: int
        number of flavours to consider
    xsize: int
        size of the grid on x
    ndata: int
        number of experimental datapoints
    n_combinations: int
        number of combinations of flavours to take into account
        default: flavs*flavs (all)
    """
    # If n_combinations = -1 combinations will include al possible combinations
    af = np.arange(flavs)
    all_combinations = np.array(np.meshgrid(af, af)).T.reshape(-1, 2)
    # See whether we want to only keep some of the combinations
    if n_combinations is None:
        combinations = all_combinations
        lc = len(combinations)
    else:
        bchoice = sorted(
            np.random.choice(
                np.arange(flavs * flavs), size=n_combinations, replace=False
            )
        )
        combinations = [all_combinations[i] for i in bchoice]
        lc = n_combinations
    # Generate random FK table and PDF
    fktable = np.random.rand(ndata, lc, xsize, xsize)
    return fktable, np.array(combinations)


# Generate an FK table, PDF and combinations list for DIS
def generate_input_DIS(flavs=3, xsize=2, ndata=5, n_combinations=-1):
    """ Generates fake input (fktable and array of combinations) for the DIS convolution

    Parameters
    ----------
    flavs: int
        number of flavours to consider
    xsize: int
        size of the grid on x
    ndata: int
        number of experimental datapoints
    n_combinations: int
        number of combinations of flavours to take into account
        default: flavs (all)
    """
    all_combinations = np.arange(flavs)
    # Now see whether we want to only keep some of the combinations
    if n_combinations is None:
        combinations = all_combinations
        lc = len(combinations)
    else:
        combinations = sorted(
            np.random.choice(all_combinations, size=n_combinations, replace=False)
        )
        lc = n_combinations
    # Generate random FK table and PDF
    fktable = np.random.rand(ndata, lc, xsize)
    return fktable, np.array(combinations)


# Tests
def test_DIS_basis():
    fk, comb = generate_input_DIS(
        flavs=FLAVS, xsize=XSIZE, ndata=NDATA, n_combinations=FLAVS - 1
    )
    obs_layer = layers.DIS(NDATA, fk, basis=comb, nfl=FLAVS)
    # Get the basis from the layer
    result = obs_layer.basis
    # Compute the basis with numpy
    reference = np.zeros(FLAVS, dtype=bool)
    for i in comb:
        reference[i] = True
    assert np.alltrue(result == reference)


def test_DY_basis():
    fk, comb = generate_input_had(
        flavs=FLAVS, xsize=XSIZE, ndata=NDATA, n_combinations=FLAVS
    )
    obs_layer = layers.DY(NDATA, fk, basis=comb, nfl=FLAVS)
    # Get the basis from the layer
    result = obs_layer.basis
    # Compute the basis with numpy
    reference = comb.reshape(-1, 2)
    assert np.alltrue(result == reference)


def test_DIS():
    # Input values
    fk, comb = generate_input_DIS(
        flavs=FLAVS, xsize=XSIZE, ndata=NDATA, n_combinations=FLAVS - 1
    )
    pdf = np.random.rand(XSIZE, FLAVS)
    kp = op.numpy_to_tensor(np.expand_dims(pdf, 0))
    # generate the n3fit results
    obs_layer = layers.DIS(NDATA, fk, basis=comb, nfl=FLAVS)
    result_tensor = obs_layer(kp)
    result = op.evaluate(result_tensor)
    # Compute the numpy version of this layer
    basis = obs_layer.basis
    pdf_masked = pdf.T[basis].T
    reference = np.tensordot(fk, pdf_masked, axes=[[2, 1], [0, 1]])
    assert np.allclose(result, reference, THRESHOLD)


def test_DY():
    # Input values
    fk, comb = generate_input_had(
        flavs=FLAVS, xsize=XSIZE, ndata=NDATA, n_combinations=FLAVS - 1
    )
    pdf = np.random.rand(XSIZE, FLAVS)
    kp = op.numpy_to_tensor(np.expand_dims(pdf, 0))
    # generate the n3fit results
    obs_layer = layers.DY(NDATA, fk, basis=comb, nfl=FLAVS)
    result_tensor = obs_layer(kp)
    result = op.evaluate(result_tensor)
    # Compute the numpy version of this layer
    mask = np.zeros((FLAVS, FLAVS), dtype=bool)
    for i, j in obs_layer.basis:
        mask[i, j] = True
    lumi = np.tensordot(pdf, pdf, axes=0)
    lumi_perm = np.moveaxis(lumi, [1, 3], [0, 1])
    lumi_masked = lumi_perm[mask]
    reference = np.tensordot(fk, lumi_masked, axes=3)
    assert np.allclose(result, reference, THRESHOLD)


def test_rotation():
    # Input dictionary to build the rotation matrix using vp2 functions
    flav_info = [
    {'fl':'u'},
    {'fl':'ubar'},
    {'fl':'d'},
    {'fl':'dbar'},
    {'fl':'s'},
    {'fl':'sbar'},
    {'fl':'c'},
    {'fl':'g'},
]
    # Apply the rotation using numpy tensordot
    x = np.ones(8) # Vector in the flavour basis v_i
    x = np.expand_dims(x,axis=[0,1]) # Give to the input the shape (1,1,8)
    mat = rotation(flav_info) # Rotation matrix R_ij, i=flavour, j=evolution
    res_np = np.tensordot(x,mat,(2,0)) # Vector in the evolution basis u_j=R_ij*vi
    
    # Apply the rotation through the rotation layer
    x = op.numpy_to_tensor(x)
    rotmat = layers.FlavourToEvolution(flav_info)
    res_layer = rotmat(x)
    assert np.alltrue(res_np==res_layer)
    