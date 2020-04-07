"""
    This module includes rotation layers
"""
from n3fit.backends import MetaLayer
from validphys.pdfbases import rotation

class FlavourToEvolution(MetaLayer):
    """ 
        Rotates from the flavour basis to
        the evolution basis. 
    """
    def __init__(
        self,
        flav_info,
        **kwargs,
    ):
        rotation_matrix = rotation(flav_info)
        self.rotation_matrix = self.np_to_tensor(rotation_matrix)   
        super().__init__(**kwargs)
    
    def call(self, x_raw):    
        x_flav = self.transpose(x_raw)
        pdf_evol = self.tensor_product(self.rotation_matrix, x_flav, 1)
        return self.reshape(pdf_evol, x_raw.shape)      
        


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
