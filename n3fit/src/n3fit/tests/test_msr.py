import numpy as np

from n3fit.backends import operations as op
from n3fit.layers import MSR_Normalization


def apply_layer_to_fixed_input(layer):
    np.random.seed(422)
    pdf_integrated = op.numpy_to_tensor(np.random.normal(size=(1, 1, 14)))

    photon_integral = op.numpy_to_tensor([np.random.normal(size=(1, 1))])

    return layer(pdf_integrated, photon_integral)


def test_all():
    layer = MSR_Normalization(mode='ALL')
    output = apply_layer_to_fixed_input(layer)
    known_output = op.numpy_to_tensor(
        [
            [
                [
                    1.0,
                    1.0,
                    6.7740397,
                    -3.8735497,
                    15.276561,
                    5.9522753,
                    8.247783,
                    -3.8735497,
                    -3.8735497,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                ]
            ]
        ]
    )
    np.testing.assert_allclose(output, known_output, rtol=1e-5)


def test_msr():
    layer = MSR_Normalization(mode='MSR')
    output = apply_layer_to_fixed_input(layer)
    known_output = op.numpy_to_tensor(
        [[[1.0, 1.0, 6.7740397, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]]]
    )
    np.testing.assert_allclose(output, known_output, rtol=1e-5)


def test_vsr():
    layer = MSR_Normalization(mode='VSR')
    output = apply_layer_to_fixed_input(layer)
    known_output = op.numpy_to_tensor(
        [
            [
                [
                    1.0,
                    1.0,
                    1.0,
                    -3.8735497,
                    15.276561,
                    5.9522753,
                    8.247783,
                    -3.8735497,
                    -3.8735497,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                ]
            ]
        ]
    )
    np.testing.assert_allclose(output, known_output, rtol=1e-5)
