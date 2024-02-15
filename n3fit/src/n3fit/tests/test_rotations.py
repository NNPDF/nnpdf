import numpy as np

from n3fit.backends import operations as op
from n3fit.layers import FkRotation, FlavourToEvolution


def test_fk():
    rotation = FkRotation()
    gridpoints = 2
    np.random.seed(0)
    pdf = op.numpy_to_tensor(np.random.rand(1, 1, gridpoints, 9))
    pdf_rotated = rotation(pdf)
    # extract single replica
    pdf_rotated = pdf_rotated[:, 0]
    pdf_rotated_known = op.numpy_to_tensor(
        [
            [
                [
                    0.0,
                    0.5488135,
                    0.71518934,
                    0.60276335,
                    0.5448832,
                    0.4236548,
                    0.96366274,
                    0.60276335,
                    0.60276335,
                    0.6458941,
                    0.4375872,
                    -3.0182784,
                    0.5488135,
                    0.5488135,
                ],
                [
                    0.0,
                    0.3834415,
                    0.79172504,
                    0.5288949,
                    0.56804454,
                    0.92559665,
                    0.83261985,
                    0.5288949,
                    0.5288949,
                    0.07103606,
                    0.0871293,
                    0.30256793,
                    0.3834415,
                    0.3834415,
                ],
            ]
        ]
    )
    np.testing.assert_allclose(pdf_rotated.numpy(), pdf_rotated_known.numpy(), rtol=1e-5)
