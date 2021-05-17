from n3fit.backends import MetaLayer
from n3fit.backends import operations as op

GLUON_IDX = [[2]]
V_IDX = [[3], [6], [7], [8]]
V3_IDX = [[4]]
V8_IDX = [[5]]


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
            idx += GLUON_IDX
        if self._vsr_enabled:
            idx += V_IDX + V3_IDX + V8_IDX

        self._out_scatter = op.as_layer(
            op.scatter_to_one, op_kwargs={"indices": idx, "output_dim": output_dim}
        )

        super().__init__(**kwargs, name="normalizer")

    def call(self, pdf_integrated):
        """Imposes the valence and momentum sum rules:
        A_g = (1-sigma)/g
        A_v  = A_v15 = A_v24 = A_v35 = 3/V
        A_v3 = 1/V_3
        A_v8 = 3/V_8

        Note that both the input and the output are in the 14-flavours fk-basis
        """
        y = op.flatten(pdf_integrated)
        norm_constants = []

        if self._msr_enabled:
            n_ag = [(1.0 - y[GLUON_IDX[0][0]-1]) / y[GLUON_IDX[0][0]]] * len(GLUON_IDX)
            norm_constants += n_ag

        if self._vsr_enabled:
            n_av = [3.0 / y[V_IDX[0][0]]] * len(V_IDX)
            n_av3 = [1.0 / y[V3_IDX[0][0]]] * len(V3_IDX)
            n_av8 = [3.0 / y[V8_IDX[0][0]]] * len(V8_IDX)
            norm_constants += n_av + n_av3 + n_av8

        return self._out_scatter(norm_constants)
