"""
    This module contains layers acting on the x-grid input of the NN

    The three operations included are:
        - ``xDivide``
        - ``xIntegrator``

    The names are self-describing. The only subtlety is that they do not act equally
    for all flavours. The choice of flavours on which to act in a different way is given
    as an input argument.
"""

from n3fit.backends import MetaLayer
from n3fit.backends import operations as op

BASIS_SIZE=14

class xDivide(MetaLayer):
    """
    Divide some PDFs by x

    By default it utilizes the 14-flavour FK basis and divides [v, v3, v8, v15]
    which corresponds to indices (3,4,5, 6) from
    (photon, sigma, g, v, v3, v8, v15, v24, v35, t3, t8, t15, t24, t35)

    Parameters:
    -----------
        output_dim: int
            dimension of the pdf
        div_list: list
            list of indices to be divided by X (by default [3,4,5, 6]; [v, v3, v8, v15]
    """

    def __init__(self, output_dim=BASIS_SIZE, div_list=None, **kwargs):
        if div_list is None:
            div_list = [3, 4, 5, 6]
        self.output_dim = output_dim
        self.div_list = div_list
        super().__init__(**kwargs)

    def call(self, x):
        out_array = []
        one = op.tensor_ones_like(x)
        for i in range(self.output_dim):
            if i in self.div_list:
                res = one / x
            else:
                res = one
            out_array.append(res)
        out_tensor = op.concatenate(out_array)
        return out_tensor


class xIntegrator(MetaLayer):
    """
    This layer performs a sum of the input layer/tensor on the first axis

    Receives as input a rank-n (n > 1) tensor `x` (batch_dims ..., xpoints, flavours)
    and returns a summation on the `xpoints` index (i.e., index -2)
    weighted by the weights of the grid

    Parameters
    ----------
        grid_weights: np.array
            weights of the grid
    """

    def __init__(self, grid_weights, output_dim=BASIS_SIZE, **kwargs):
        grid_weights_tensor = op.numpy_to_tensor(grid_weights)
        # Open up the grid weights
        self.grid_weights = op.many_replication(grid_weights_tensor, output_dim, axis=1)
        super().__init__(**kwargs)

    def call(self, x):
        xx = x * self.grid_weights
        return op.sum(xx, axis=-2)
