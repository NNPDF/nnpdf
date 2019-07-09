from n3fit.backends import MetaLayer

class MSR_Normalization(MetaLayer):
    """
        Applies the normalisation so that the PDF output fullfills the sum rules
    """
    def __init__(self, output_dim=14, **kwargs):
        self.output_dim = output_dim
        self.one = self.tensor_ones((1, 1))
        self.three = 3 * self.tensor_ones((1, 1))
        super(MSR_Normalization, self).__init__(**kwargs, name="normalizer")

    def compute_output_shape(self, input_shape):
        return (self.output_dim,)

    def call(self, x):
        """
            Receives as input a tensor with the value of the MSR for each PDF
            and returns a rank-1 tensor with the normalization factor A_i of each flavour
        """
        pdf_sr = self.concatenate(
            [
                self.one,  # photon
                self.one,  # sigma
                (self.one - x[0]) / x[1],  # g
                self.three / x[2],  # v
                self.one / x[3],  # v3
                self.three / x[4],  # v8
                self.three / x[2],  # v15
                self.three / x[2],  # v24
                self.three / x[2],  # v35
                self.one,  # t3
                self.one,  # t8
                self.one,  # t15 (c-)
                self.one,  # t24
                self.one,  # t35
            ]
        )
        return pdf_sr
