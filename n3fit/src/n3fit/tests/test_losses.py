"""
    Test the losses layers
"""
import numpy as np
from n3fit.layers import losses
from .test_backend import are_equal, DIM

ARR1 = np.random.rand(DIM)
ARR2 = np.random.rand(DIM)
C = np.random.rand(DIM, DIM)
INVCOVMAT = np.linalg.inv(C @ C.T)

# Tests loss functions
def test_l_invcovmat():
    loss_f = losses.LossInvcovmat(INVCOVMAT, ARR1)
    # Add a replica and batch dimension to T2
    result = loss_f(np.expand_dims(ARR2, [0, 1]))
    y = ARR1 - ARR2
    tmp = np.dot(INVCOVMAT, y)
    reference = np.dot(y, tmp)
    are_equal(result, reference, threshold=1e-4)


def test_l_positivity():
    alpha = 1e-7
    loss_f = losses.LossPositivity(alpha=alpha)
    result = loss_f(np.expand_dims(ARR2, [0, 1]))

    def elu_sum(yarr_in):
        """Applies Exponential Linear Unit
        to an array and sums it up"""
        yarr = -yarr_in
        res = 0.0
        for y in yarr:
            if y > 0:
                res += y
            else:
                res += alpha * (np.exp(y) - 1)
        return res

    reference = elu_sum(ARR1)
    are_equal(result, reference)
