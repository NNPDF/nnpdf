"""
    This module contains layers acting on the x-grid input of the NN

    The three operations included are:
        - `xDivide`
        - `xMultiply`
        - `xIntegrator`

    The names are self-describing. The only subtlety is that they do not act equally
    for all flavours. The choice of flavours on which to act in a different way is given
    as an input argument.
"""

from n3fit.backends import MetaLayer


class xDivide(MetaLayer):
    """
        Divide the pdf by x
        By default: dimension = 8 dividing the entries corresponding to v, v3, v8

        Parameters:
        -----------
            output_dim: int
                dimension of the pdf
            div_list: list
                list of indices to be divided by X (by default [2,3,4]; [v, v3, v8]
    """

    def __init__(self, output_dim=8, div_list=None, **kwargs):
        if div_list is None:
            div_list = [2, 3, 4]
        self.output_dim = output_dim
        self.div_list = div_list
        super(MetaLayer, self).__init__(**kwargs)

    def call(self, x):
        out_array = []
        one = self.tensor_ones_like(x)
        for i in range(self.output_dim):
            if i in self.div_list:
                res = one / x
            else:
                res = one
            out_array.append(res)
        out_tensor = self.concatenate(out_array, axis=-1)
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

    def __init__(self, grid_weights, **kwargs):
        grid_weights_tensor = self.np_to_tensor(grid_weights)
        # Open up the grid weights
        self.grid_weights = self.many_replication(grid_weights_tensor, replications=8, axis=1)
        super(MetaLayer, self).__init__(**kwargs)

    def call(self, x):
        xx = x * self.grid_weights
        return self.sum(xx, axis=-2)
