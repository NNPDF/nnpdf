from n3fit.backends import MetaLayer
from n3fit.backends import operations as op


class MSR_Normalization(MetaLayer):
    """
    Applies the normalisation so that the PDF output fullfills the sum rules
    """

    def __init__(self, output_dim=14, mode="ALL", **kwargs):
        self.output_dim = output_dim
        self.one = op.numpy_to_tensor([[1.0]])
        self.three = op.numpy_to_tensor([[3.0]])
        if mode == True or mode.upper() == "ALL":
            self._call_f = self._impose_all
        elif mode.upper() == "MSR":
            self._call_f = self._impose_msr
        elif mode.upper() == "VSR":
            self._call_f = self._impose_vsr

        super(MSR_Normalization, self).__init__(**kwargs, name="normalizer")

    def _impose_all(self, x):
        """
        Receives as input a tensor with the value of the MSR for each PDF
        and returns a rank-1 tensor with the normalization factor A_i of each flavour
        """
        pdf_sr = op.concatenate(
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

    def _impose_msr(self, x):
        """
        Imposes only the sum rules for the gluon
        """
        pdf_sr = op.concatenate(
            [
                self.one,  # photon
                self.one,  # sigma
                (self.one - x[0]) / x[1],  # g
                self.one,  # v
                self.one,  # v3
                self.one,  # v8
                self.one,  # v15
                self.one,  # v24
                self.one,  # v35
                self.one,  # t3
                self.one,  # t8
                self.one,  # t15 (c-)
                self.one,  # t24
                self.one,  # t35
            ]
        )
        return pdf_sr

    def _impose_vsr(self, x):
        """
        Imposes only the valence sum rules
        """
        pdf_sr = op.concatenate(
            [
                self.one,  # photon
                self.one,  # sigma
                self.one,  # g
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

    def call(self, xgrid):
        x = op.flatten(xgrid)
        return self._call_f(x)
