from n3fit.backends import MetaLayer
from n3fit.backends import operations as op

IDX = {
    'photon': 0,
    'sigma': 1,
    'g': 2,
    'v': 3,
    'v3': 4,
    'v8': 5,
    'v15': 6,
    'v24': 7,
    'v35': 8,
}
MSR_INDICES = [IDX['g']]
VSR_INDICES = [IDX[f] for f in ['v', 'v35', 'v24', 'v3', 'v8', 'v15']]


class MSR_Normalization(MetaLayer):
    """
    Applies the normalisation so that the PDF output fullfills the sum rules
    """

    _msr_enabled = False
    _vsr_enabled = False

    def __init__(self, mode="ALL", **kwargs):
        if mode == True or mode.upper() == "ALL":
            self._msr_enabled = True
            self._vsr_enabled = True
        elif mode.upper() == "MSR":
            self._msr_enabled = True
        elif mode.upper() == "VSR":
            self._vsr_enabled = True
        else:
            raise ValueError(f"Mode {mode} not accepted for sum rules")

        self.indices = []
        if self._msr_enabled:
            self.indices += MSR_INDICES
        if self._vsr_enabled:
            self.indices += VSR_INDICES
        self.indices = [[i] for i in self.indices]

        super().__init__(**kwargs)

    def build(self, input_shape):
        out_shape = input_shape[1:]
        self._out_scatter = lambda pdf_integrated: op.scatter_to_one(
            pdf_integrated, indices=self.indices, output_shape=out_shape
        )
        super().build(input_shape)

    def call(self, pdf_integrated):
        """Imposes the valence and momentum sum rules:
        A_g = (1-sigma-photon)/g
        A_v = A_v24 = A_v35 = 3/V
        A_v3 = 1/V_3
        A_v8 = 3/V_8
        A_v15 = 3/V_15

        Note that both the input and the output are in the 14-flavours fk-basis

        Parameters
        ----------
        pdf_integrated: (Tensor(1, 14))
            the integrated PDF

        Returns
        -------
        normalization_factor: Tensor(14)
            The normalization factors per flavour.
        """
        y = pdf_integrated[0]  # get rid of the batch dimension
        norm_constants = []

        if self._msr_enabled:
            n_ag = [(1.0 - y[IDX['sigma']] - y[IDX['photon']]) / y[IDX['g']]]
            norm_constants += n_ag

        if self._vsr_enabled:
            n_av = [3.0 / y[IDX['v']]] * 3
            n_av3 = [1.0 / y[IDX['v3']]]
            n_av8 = [3.0 / y[IDX['v8']]]
            n_av15 = [3.0 / y[IDX['v15']]]
            norm_constants += n_av + n_av3 + n_av8 + n_av15

        return self._out_scatter(norm_constants)
