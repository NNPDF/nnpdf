"""
Tests for the layers of n3fit
This module checks that the layers do what they would do with numpy
"""

import dataclasses

import numpy as np

from n3fit.backends import operations as op
import n3fit.layers as layers
from validphys.loader import Loader
from validphys.pdfbases import fitbasis_to_NN31IC

FLAVS = 3
XSIZE = 4
NDATA = 3
THRESHOLD = 1e-6

PARAMS = {
    "dataset_name": "NULL",
    "operation_name": "NULL",
    "nfl": FLAVS,
    "boundary_condition": None,
}


@dataclasses.dataclass
class _fake_FKTableData:
    """Fake validphys.coredata.FKTableData to be used in the tests"""

    fktable: np.array
    luminosity_mapping: np.array
    xgrid: np.array
    convolution_types: tuple = ("UnpolPDF",)

    @property
    def hadronic(self):
        return len(self.convolution_types) == 2


# Helper functions
def generate_input_had(flavs=3, xsize=2, ndata=4, n_combinations=None):
    """Generates fake input (fktable and array of combinations) for the hadronic convolution

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
            np.random.choice(np.arange(flavs * flavs), size=n_combinations, replace=False)
        )
        combinations = [all_combinations[i] for i in bchoice]
        lc = n_combinations
    # Generate random FK table and PDF
    fktable = np.random.rand(ndata, lc, xsize, xsize)
    return fktable, np.array(combinations)


# Generate an FK table, PDF and combinations list for DIS
def generate_input_DIS(flavs=3, xsize=2, ndata=5, n_combinations=-1):
    """Generates fake input (fktable and array of combinations) for the DIS convolution

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


def generate_DIS(nfk=1):
    fktables = []
    for _ in range(nfk):
        fk, comb = generate_input_DIS(
            flavs=FLAVS, xsize=XSIZE, ndata=NDATA, n_combinations=FLAVS - 1
        )
        fktables.append(_fake_FKTableData(fk, comb, np.ones((1, XSIZE))))
    return fktables


def generate_had(nfk=1):
    fktables = []
    for _ in range(nfk):
        fk, comb = generate_input_had(flavs=FLAVS, xsize=XSIZE, ndata=NDATA, n_combinations=FLAVS)
        fktables.append(
            _fake_FKTableData(
                fk, comb, np.ones((1, XSIZE)), convolution_types=("UnpolPDF", "UnpolPDF")
            )
        )
    return fktables


# Tests
def test_DIS_basis():
    fktables = generate_DIS(2)
    fks = [i.fktable for i in fktables]
    obs_layer = layers.DIS(fktables, fks, **PARAMS)
    # Get the masks from the layer
    all_masks = obs_layer.all_masks
    for result, fk in zip(all_masks, fktables):
        comb = fk.luminosity_mapping
        # Compute the basis with numpy
        reference = np.zeros(FLAVS, dtype=bool)
        for i in comb:
            reference[i] = True
        np.testing.assert_allclose(result, reference)


def test_DY_basis():
    fktables = generate_had(2)
    fks = [i.fktable for i in fktables]
    obs_layer = layers.DY(fktables, fks, **PARAMS)
    # Get the mask from the layer
    all_masks = obs_layer.all_masks
    for result, fk in zip(all_masks, fktables):
        comb = fk.luminosity_mapping
        reference = np.zeros((FLAVS, FLAVS))
        for i, j in comb:
            reference[i, j] = True
        np.testing.assert_allclose(result, reference)


def test_DIS():
    tests = [(2, "ADD"), (1, "NULL")]
    for nfk, ope in tests:
        # Input values
        kwargs = dict(PARAMS)
        kwargs["operation_name"] = ope
        fktables = generate_DIS(nfk)
        fks = [i.fktable for i in fktables]
        obs_layer = layers.DIS(fktables, fks, **kwargs)
        pdf = np.random.rand(XSIZE, FLAVS)
        kp = op.numpy_to_tensor([[pdf]])  # add batch and replica dimension
        # generate the n3fit results
        result_tensor = obs_layer(kp)
        result = op.tensor_to_numpy_or_python(result_tensor)
        # Compute the numpy version of this layer
        all_masks = obs_layer.all_masks
        if len(all_masks) < nfk:
            all_masks *= nfk
        reference = 0
        for fktabledata, mask in zip(fktables, all_masks):
            fk = fktabledata.fktable
            pdf_masked = pdf.T[mask.numpy()].T
            reference += np.tensordot(fk, pdf_masked, axes=[[2, 1], [0, 1]])
        assert np.allclose(result, reference, THRESHOLD)


def test_DY():
    tests = [(2, "ADD"), (1, "NULL")]
    for nfk, ope in tests:
        # Input values
        kwargs = dict(PARAMS)
        kwargs["operation_name"] = ope
        fktables = generate_had(nfk)
        fks = [i.fktable for i in fktables]
        obs_layer = layers.DY(fktables, fks, **kwargs)
        pdf = np.random.rand(XSIZE, FLAVS)
        kp = op.numpy_to_tensor([[pdf]])  # add batch and replica dimension
        # generate the n3fit results
        result_tensor = obs_layer(kp)
        result = op.tensor_to_numpy_or_python(result_tensor)
        # Compute the numpy version of this layer
        all_masks = obs_layer.all_masks
        if len(all_masks) < nfk:
            all_masks *= nfk
        reference = 0
        for fktabledata, mask in zip(fktables, all_masks):
            fk = fktabledata.fktable
            lumi = np.tensordot(pdf, pdf, axes=0)
            lumi_perm = np.moveaxis(lumi, [1, 3], [0, 1])
            lumi_masked = lumi_perm[mask.numpy()]
            reference += np.tensordot(fk, lumi_masked, axes=3)
        assert np.allclose(result, reference, THRESHOLD)


def test_rotation_flavour():
    # Input dictionary to build the rotation matrix using vp2 functions
    flav_info = [
        {"fl": "u"},
        {"fl": "ubar"},
        {"fl": "d"},
        {"fl": "dbar"},
        {"fl": "s"},
        {"fl": "sbar"},
        {"fl": "c"},
        {"fl": "g"},
    ]
    # Apply the rotation using numpy tensordot
    pdf = np.ones(8)  # Vector in the flavour basis v_i
    pdf = np.expand_dims(pdf, axis=[0, 1, 2])  # Add batch, replica, x dimensions
    mat = fitbasis_to_NN31IC(flav_info, "FLAVOUR")  # Rotation matrix R_ij, i=flavour, j=evolution
    res_np = np.tensordot(pdf, mat, (3, 0))  # Vector in the evolution basis u_j=R_ij*vi

    # Apply the rotation through the rotation layer
    pdf = op.numpy_to_tensor(pdf)
    rotmat = layers.FlavourToEvolution(flav_info, "FLAVOUR")
    res_layer = rotmat(pdf)
    np.testing.assert_allclose(res_np, res_layer)


def test_rotation_evol():
    # Input dictionary to build the rotation matrix using vp2 functions
    flav_info = [
        {"fl": "sng"},
        {"fl": "v"},
        {"fl": "v3"},
        {"fl": "v8"},
        {"fl": "t3"},
        {"fl": "t8"},
        {"fl": "t15"},
        {"fl": "g"},
    ]
    # Apply the rotation using numpy tensordot
    pdf = np.ones(8)  # Vector in the flavour basis v_i
    pdf = np.expand_dims(pdf, axis=[0, 1, 2])  # Add batch, replica, x dimensions
    mat = fitbasis_to_NN31IC(flav_info, "EVOL")  # Rotation matrix R_ij, i=flavour, j=evolution
    res_np = np.tensordot(pdf, mat, (3, 0))  # Vector in the evolution basis u_j=R_ij*vi

    # Apply the rotation through the rotation layer
    pdf = op.numpy_to_tensor(pdf)
    rotmat = layers.FlavourToEvolution(flav_info, "EVOL")
    res_layer = rotmat(pdf)
    np.testing.assert_allclose(res_np, res_layer)


def test_mask():
    """Test the mask layer"""
    batch_size, replicas, points = 1, 1, 100
    shape = (batch_size, replicas, points)
    fi = np.random.rand(*shape)
    # Check that the multiplier works
    vals = [0.0, 2.0, np.random.rand()]
    for val in vals:
        masker = layers.Mask(c=val)
        ret = masker(fi)
        np.testing.assert_allclose(ret, val * fi, rtol=1e-5)
    # Check that the boolean works
    np_mask = np.random.randint(0, 2, size=shape[1:], dtype=bool)
    masker = layers.Mask(bool_mask=np_mask)
    ret = masker(fi)
    masked_fi = fi[np.newaxis, :, np_mask]
    np.testing.assert_allclose(ret, masked_fi, rtol=1e-5)
    # Check that the combination works!
    rn_val = vals[-1]
    masker = layers.Mask(bool_mask=np_mask, c=rn_val)
    ret = masker(fi)
    np.testing.assert_allclose(ret, masked_fi * rn_val, rtol=1e-5)


def test_addphoton_init():
    """Test AddPhoton class."""
    addphoton = layers.AddPhoton(photons=None)
    np.testing.assert_equal(addphoton._photons_generator, None)
    addphoton = layers.AddPhoton(photons=1234)
    np.testing.assert_equal(addphoton._photons_generator, 1234)
    np.testing.assert_equal(addphoton._pdf_ph, None)


class FakePhoton:
    def __call__(self, xgrid):
        return [np.exp(-xgrid)]


def test_compute_photon():
    photon = FakePhoton()
    addphoton = layers.AddPhoton(photons=photon)
    xgrid = np.geomspace(1e-4, 1.0, 10)
    addphoton.register_photon(xgrid)
    np.testing.assert_allclose(addphoton._pdf_ph, [np.exp(-xgrid)])


def test_computation_bc():
    """Test the computation of the boundary conditions."""
    n_replicas = 25
    xgrid = np.geomspace(1e-4, 1.0, num=100)
    pdf = Loader().check_pdf("NNPDF40_nnlo_as_01180")
    respdf_bc = layers.observable.compute_pdf_boundary(
        pdf=pdf, q0_value=10.0, xgrid=xgrid, n_std=0.0, n_replicas=n_replicas
    )
    exp_shape = [1, n_replicas, xgrid.size, 14]  # (batch, replicas, x, flavours)
    np.testing.assert_allclose(respdf_bc.shape.as_list(), exp_shape)
