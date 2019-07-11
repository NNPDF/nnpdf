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

        # Arguments:
            - `output_dim`: dimension of the pdf
            - `div_list`: list of indices to be divided by X (by default [2,3,4]; [v, v3, v8]
    """

    def __init__(self, output_dim=8, div_list=None, **kwargs):
        if div_list is None:
            div_list = [2, 3, 4]
        self.output_dim = output_dim
        self.div_list = div_list
        super(MetaLayer, self).__init__(**kwargs)

    def compute_output_shape(self, input_shape):
        return (input_shape[0], self.output_dim)

    def call(self, x):
        out_array = []
        one = self.tensor_ones_like(x)
        for i in range(self.output_dim):
            if i in self.div_list:
                res = one / x
            else:
                res = one
            out_array.append(res)
        out_tensor = self.concatenate(out_array, axis=1)
        return out_tensor


class xMultiply(MetaLayer):
    """
        Multiply the pdf by x
        By default: dimension = 8 multiply all entries but v, v3, v8

        # Arguments:
            - `output_dim`: dimension of the pdf
            - `not_mul_list`: list of indices *not* to multiply by X (by default [2,3,4]; [v, v3, v8]
    """

    def __init__(self, output_dim=8, not_mul_list=None, **kwargs):
        if not_mul_list is None:
            not_mul_list = [2, 3, 4]
        self.output_dim = output_dim
        self.not_mul_list = not_mul_list
        super(MetaLayer, self).__init__(**kwargs)

    def compute_output_shape(self, input_shape):
        return (input_shape[0], self.output_dim)

    def call(self, x):
        out_array = []
        one = self.tensor_ones_like(x)
        for i in range(self.output_dim):
            if i in self.not_mul_list:
                res = one
            else:
                res = one * x
            out_array.append(res)
        out_tensor = self.concatenate(out_array, axis=1)
        return out_tensor


class xIntegrator(MetaLayer):
    """
    This layer performs a sum of the input layer/tensor on the first axis
    """

    def __init__(self, grid_weights, output_dim=8, **kwargs):
        grid_weights_tensor = self.np_to_tensor(grid_weights)
        # Open up the grid weights
        self.grid_weights = self.many_replication(grid_weights_tensor, replications=8, axis=1)
        self.output_dim = output_dim
        super(MetaLayer, self).__init__(**kwargs)

    def compute_output_shape(self, input_shape):
        return (self.output_dim,)

    def call(self, x):
        """
            Receives as input a rank-2 tensor `x` (xpoints, flavours)
            and returns a summation on the first index (xpoints) of tensor `x`, weighted by the
            weights of the grid (in the most common case, 1/grid_points)
        """
        xx = x * self.grid_weights
        return self.sum(xx, axis=0)
