from n3fit.backends import MetaLayer
from n3fit.backends import operations as op

V_IDX = [[3], [7], [8]]

indices = {
    'photon': 0,
    'sigma': 1,
    'g': 2,
    'v': 3,
    'v3': 4,
    'v8': 5,
    'v15': 6,
    'v35': 7,
    'v24': 8,
}



class MSR_Normalization(MetaLayer):
    """
    Applies the normalisation so that the PDF output fullfills the sum rules
    """

    _msr_enabled = False
    _vsr_enabled = False

    def __init__(self, output_dim=14, mode="ALL", **kwargs):
        if mode == True or mode.upper() == "ALL":
            self._msr_enabled = True
            self._vsr_enabled = True
        elif mode.upper() == "MSR":
            self._msr_enabled = True
        elif mode.upper() == "VSR":
            self._vsr_enabled = True
        else:
            raise ValueError(f"Mode {mode} not accepted for sum rules")

        idx = []
        if self._msr_enabled:
            idx += [indices['g']]
        if self._vsr_enabled:
            idx += [indices[f] for f in ['v', 'v35', 'v24', 'v3', 'v8', 'v15']]
        idx = [[i] for i in idx]

        self._out_scatter = op.as_layer(
            op.scatter_to_one, op_kwargs={"indices": idx, "output_dim": output_dim}
        )

        super().__init__(**kwargs)

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
        pdf_integrated: (Tensor(1,None,14))
            the integrated PDF

        Returns
        -------
        normalization_factor: Tensor(14)
            The normalization factors per flavour.
        """
        y = op.flatten(pdf_integrated)
        norm_constants = []

        if self._msr_enabled:
            n_ag = [(1.0 - y[indices['sigma']] - y[indices['photon']]) / y[indices['g']]]
            norm_constants += n_ag

        if self._vsr_enabled:
            n_av = [3.0 / y[indices['v']]] * 3
            n_av3 = [1.0 / y[indices['v3']]]
            n_av8 = [3.0 / y[indices['v8']]]
            n_av15 = [3.0 / y[indices['v15']]]
            norm_constants += n_av + n_av3 + n_av8 + n_av15

        return self._out_scatter(norm_constants)
