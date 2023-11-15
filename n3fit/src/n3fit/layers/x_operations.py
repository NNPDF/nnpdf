"""
    This module contains layers acting on the x-grid input of the NN

    The two operations included are:
        - ``xDivide``
        - ``xIntegrator``

    The names are self-describing. The only subtlety is that they do not act equally
    for all flavours. The choice of flavours on which to act in a different way is given
    as an input argument.
"""
from typing import List, Optional

from n3fit.backends import MetaLayer
from n3fit.backends import operations as op

BASIS_SIZE = 14


class xDivide(MetaLayer):
    """
    Create tensor of either 1/x or ones depending on the flavour,
    to be used to divide some PDFs by x by multiplying with the result.

    By default it utilizes the 14-flavour FK basis and divides [v, v3, v8, v15]
    which corresponds to indices (3, 4, 5, 6) from
    (photon, sigma, g, v, v3, v8, v15, v24, v35, t3, t8, t15, t24, t35)

    Parameters:
    -----------
        output_dim: int
            dimension of the pdf
        div_list: list
            list of indices to be divided by X (by default [3, 4, 5, 6]; [v, v3, v8, v15]
    """

    def __init__(
        self, output_dim: int = BASIS_SIZE, div_list: Optional[List[int]] = None, **kwargs
    ):
        if div_list is None:
            div_list = [3, 4, 5, 6]
        self.output_dim = output_dim
        self.div_list = div_list
        super().__init__(**kwargs)

        self.powers = [-1 if i in div_list else 0 for i in range(output_dim)]

    def call(self, x):
        return op.pow(x, self.powers)

    def get_config(self):
        config = super().get_config()
        config.update({"output_dim": self.output_dim, "div_list": self.div_list})
        return config


class xIntegrator(MetaLayer):
    """
    This layer performs a sum of the input layer/tensor on the axis corresponding to the x-grid
    weighted by the weights of the grid.

    The output shape is the input shape with the x-axis removed.

    Parameters
    ----------
        grid_weights: np.array
            weights of the grid
        x_axis: int (default=1)
            axis of the input tensor that corresponds to the x-grid
    """

    def __init__(self, grid_weights, x_axis=1, **kwargs):
        self.x_axis = x_axis
        self.grid_weights = op.flatten(op.numpy_to_tensor(grid_weights))
        super().__init__(**kwargs)

    def call(self, pdf):
        return op.tensor_product(pdf, self.grid_weights, axes=[self.x_axis, 0])
