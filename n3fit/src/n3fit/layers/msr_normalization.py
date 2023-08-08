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
MSR_COMPONENTS = ['g']
MSR_DENOMINATORS = {
    'g': 'g',
}
# The VSR normalization factor of component f is given by
# VSR_CONSTANTS[f] / VSR_DENOMINATORS[f]
VSR_COMPONENTS = ['v', 'v35', 'v24', 'v3', 'v8', 'v15']
VSR_CONSTANTS = {
    'v': 3.0,
    'v35': 3.0,
    'v24': 3.0,
    'v3': 1.0,
    'v8': 3.0,
    'v15': 3.0,
}
VSR_DENOMINATORS = {
    'v': 'v',
    'v35': 'v',
    'v24': 'v',
    'v3': 'v3',
    'v8': 'v8',
    'v15': 'v15',
}


class MSR_Normalization(MetaLayer):
    """
    Computes the normalisation factors for the sum rules of the PDFs.
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
        self.divisor_indices = []
        if self._msr_enabled:
            self.indices += [IDX[c] for c in MSR_COMPONENTS]
            self.divisor_indices += [IDX[MSR_DENOMINATORS[c]] for c in MSR_COMPONENTS]
        if self._vsr_enabled:
            self.divisor_indices += [IDX[VSR_DENOMINATORS[c]] for c in VSR_COMPONENTS]
            self.indices += [IDX[c] for c in VSR_COMPONENTS]
            self.vsr_factors = [VSR_CONSTANTS[c] for c in VSR_COMPONENTS]
        # Need this extra dimension for the scatter_to_one operation
        self.indices = [[i] for i in self.indices]

        super().__init__(**kwargs)

    def call(self, pdf_integrated, photon_integral):
        """
        Computes the normalization factors for the PDFs:
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
        photon_integral: (Tensor(1, 1))
            the integrated photon PDF

        Returns
        -------
        normalization_factor: Tensor(14)
            The normalization factors per flavour.
        """
        y = pdf_integrated[0]  # get rid of the batch dimension
        photon_integral = photon_integral[0]  # get rid of the batch dimension
        norm_constants = []

        if self._msr_enabled:
            norm_constants += [(1.0 - y[IDX['sigma']] - photon_integral[0])]

        if self._vsr_enabled:
            norm_constants += self.vsr_factors

        divisors = op.gather(y, self.divisor_indices, axis=0)
        norm_constants = norm_constants / divisors

        # Fill in the rest of the flavours with 1
        norm_constants = op.scatter_to_one(
            norm_constants, indices=self.indices, output_shape=y.shape
        )

        return norm_constants
