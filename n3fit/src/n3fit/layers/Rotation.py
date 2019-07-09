from n3fit.backends import MetaLayer

class Rotation(MetaLayer):
    """
        Applies a transformation from the dimension-8 fit basis
        to the dimension-14 evolution basis
    """
    def __init__(self, output_dim=14, **kwargs):
        self.output_dim = output_dim
        super(MetaLayer, self).__init__(**kwargs, name="evolution")

    def compute_output_shape(self, input_shape):
        return (input_shape[0], self.output_dim)

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
        pdf = self.concatenate(pdf_raw_list, target_shape=(self.output_dim, -1))
        return self.transpose(pdf)
