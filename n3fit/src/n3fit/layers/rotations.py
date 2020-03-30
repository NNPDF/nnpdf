"""
    This module includes rotation layers
"""

from n3fit.backends import MetaLayer

class FlavourToEvolution(MetaLayer):
    """ 
        Rotates from the evolution basis to
        the evolution basis
    """
    # TODO: add a custom __init__ that defines the rotation

    def call(self, x_raw):
        # Let's decide that the input is
        # u, ubar, d, dbar, s, sbar, c, g
        # TODO: it needs to match
        #       
        x_flav = self.transpose(x_raw)
        u    = x_flav[0]
        ubar = x_flav[1]
        d    = x_flav[2]
        dbar = x_flav[3]
        s    = x_flav[4]
        sbar = x_flav[5]
        c    = x_flav[6]
        g    = x_flav[7]
        cbar = c
        sigma = u + ubar + d + dbar + s + sbar + c + cbar
        v = u - ubar + d - dbar + s - sbar + c - cbar
        v3 = u - ubar - d + dbar
        v8 = u - ubar + d - dbar - 2*s + 2*sbar
        t3 = u + ubar - d - dbar
        t8 = u + ubar + d + dbar - 2*s - 2*sbar
        pdf_evol = [sigma, g, v, v3, v8, t3, t8, c+cbar]
        ret = self.concatenate(pdf_evol, target_shape=x_raw.shape)
        return ret


class Rotation(MetaLayer):
    """
        Applies a transformation from the dimension-8 fit basis
        to the dimension-14 evolution basis
    """

    def __init__(self, output_dim=14, **kwargs):
        self.output_dim = output_dim
        super(MetaLayer, self).__init__(**kwargs, name="evolution")

    def call(self, x_raw):
        x = self.transpose(x_raw)
        pdf_raw_list = [
            0 * x[0],  # photon
            x[0],  # sigma
            x[1],  # g
            x[2],  # v
            x[3],  # v3
            x[4],  # v8
            x[2],  # v15
            x[2],  # v24
            x[2],  # v35
            x[5],  # t3
            x[6],  # t8
            x[0] - 4 * x[7],  # t15 (c-)
            x[0],  # t24
            x[0],  # t35
        ]
        return self.concatenate(pdf_raw_list)
