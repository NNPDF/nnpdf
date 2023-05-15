import numpy as np
from n3fit.layers import xDivide, xIntegrator
from n3fit.backends import MetaLayer, Lambda
from n3fit.backends import operations as op

GLUON_IDX = [[2]]
V_IDX = [[3], [7], [8]]
V3_IDX = [[4]]
V8_IDX = [[5]]
V15_IDX = [[6]]


class MSR_Normalization(MetaLayer):
    """
    Applies the normalisation so that the PDF output fullfills the sum rules
    """

    _msr_enabled = False
    _vsr_enabled = False

    def __init__(self, output_dim=14, mode="ALL", nx=int(2e3), scaler=None, **kwargs):
        self.nx = nx
        self.scaler = scaler
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
            idx += V_IDX + V3_IDX + V8_IDX + V15_IDX

        self._out_scatter = op.as_layer(
            op.scatter_to_one, op_kwargs={"indices": idx, "output_dim": output_dim}
        )

        self._gen_integration_input()
        # 2. Prepare the pdf for integration
        #    for that we need to multiply several flavours with 1/x
        self.get_original = Lambda(lambda x: op.op_gather_keep_dims(x, -1, axis=-1), name="x_original_integ")
        self.divide_by_x = xDivide()
        # 3. Now create the integration layer (the layer that will simply integrate, given some weight
        self.integrator = xIntegrator(self.weights_array, input_shape=(self.nx,))

        self.compute_integrand = Lambda(op.op_multiply, name="pdf_integrand")
        self.compute_normalized_pdf = Lambda(lambda pdf_norm: pdf_norm[0] * pdf_norm[1], name="pdf_normalized")

        super().__init__(**kwargs, name="normalizer")

    def call(self, pdf_integrated):
        """Imposes the valence and momentum sum rules:
        A_g = (1-sigma)/g
        A_v = A_v24 = A_v35 = 3/V
        A_v3 = 1/V_3
        A_v8 = 3/V_8
        A_v15 = 3/V_15

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
            n_av15 = [3.0 / y[V15_IDX[0][0]]] * len(V15_IDX)
            norm_constants += n_av + n_av3 + n_av8 + n_av15

        return self._out_scatter(norm_constants)

    def msr_impose(self):
        """
            This function receives:
            Generates a function that applies a normalization layer to the fit.
                - fit_layer: the 8-basis layer of PDF which we fit
            The normalization is computed from the direct output of the NN (so the 7,8-flavours basis)
                - final_layer: the 14-basis which is fed to the fktable
            and it is applied to the input of the fktable (i.e., to the 14-flavours fk-basis).
            It uses pdf_fit to compute the sum rule and returns a modified version of
            the final_pdf layer with a normalisation by which the sum rule is imposed

            Parameters
            ----------
                scaler: scaler
                    Function to apply to the input. If given the input to the model
                    will be a (1, None, 2) tensor where dim [:,:,0] is scaled 
        """

        # 4. Now create the normalization by selecting the right integrations

        # Finally prepare a function which will take as input the output of the PDF model
        # and will return it appropiately normalized.
        def apply_normalization(layer_pdf):
            """
                layer_pdf: output of the PDF, unnormalized, ready for the fktable
            """
            pdf_xgrid_integration = layer_pdf(self.xgrid_integration)

            def ultimate_pdf(x):
                pdf_xgrid = layer_pdf(x)
                return self.tempcall([pdf_xgrid, pdf_xgrid_integration])

            return ultimate_pdf

        return apply_normalization

    def tempcall(self, pdfx_pdfinteg):
        pdf_xgrid, pdf_xgrid_integration = pdfx_pdfinteg

        x_original = self.get_original(self.xgrid_integration)
        x_divided = self.divide_by_x(x_original)
        pdf_integrand = self.compute_integrand([x_divided, pdf_xgrid_integration])
        pdf_integrated = self.integrator(pdf_integrand)
        normalization_factor = self(pdf_integrated)

        pdf_normalized = self.compute_normalized_pdf([pdf_xgrid, normalization_factor])
        return pdf_normalized

    def _gen_integration_input(self):
        """
        Generates a np.array (shaped (nx,1)) of nx elements where the
        nx/2 first elements are a logspace between 0 and 0.1
        and the rest a linspace from 0.1 to 0
        """
        lognx = int(self.nx / 2)
        linnx = int(self.nx - lognx)
        xgrid_log = np.logspace(-9, -1, lognx + 1)
        xgrid_lin = np.linspace(0.1, 1, linnx)
        xgrid = np.concatenate([xgrid_log[:-1], xgrid_lin]).reshape(self.nx, 1)

        spacing = [0.0]
        for i in range(1, self.nx):
            spacing.append(np.abs(xgrid[i - 1] - xgrid[i]))
        spacing.append(0.0)

        weights = []
        for i in range(self.nx):
            weights.append((spacing[i] + spacing[i + 1]) / 2.0)
        self.weights_array = np.array(weights).reshape(self.nx, 1)

        # 1. Generate the fake input which will be used to integrate
        # 1b If a scaler is provided, scale the input xgrid
        if self.scaler:
            xgrid = self.scaler(xgrid)

        # 5. Make the xgrid array into a backend input layer so it can be given to the normalization
        self.xgrid_integration = op.numpy_to_input(xgrid, name="integration_grid")

