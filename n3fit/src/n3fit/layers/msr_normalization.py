"""
    Definition of the imposition of the Momentum Sum Rule and Valence Sum Rules to in the PDF fit.

    In the module level constants ``{MSR/VSR}_COMPONENTS`` the flavours affected by the MSR and VSR are defined.
    For the Valence Sum Rule instead `VSR_DENOMINATOR` defines the integral of which flavour are used
    to compute the normalization. Note that for a Nf=4 fit  `v35=v24=v`.
    If the number of flavours were to be changed in the future, this would need to be updated accordingly.
"""
import numpy as np

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
    't3': 9,
    't8': 10,
    't15': 11,
    't24': 12,
    't35': 13,
}
MSR_COMPONENTS = ['g']
MSR_DENOMINATORS = {'g': 'g'}
# The VSR normalization factor of component f is given by
# VSR_CONSTANTS[f] / VSR_DENOMINATORS[f]
VSR_COMPONENTS = ['v', 'v35', 'v24', 'v3', 'v8', 'v15']
VSR_CONSTANTS = {'v': 3.0, 'v35': 3.0, 'v24': 3.0, 'v3': 1.0, 'v8': 3.0, 'v15': 3.0}
VSR_DENOMINATORS = {'v': 'v', 'v35': 'v', 'v24': 'v', 'v3': 'v3', 'v8': 'v8', 'v15': 'v15'}

CSR_COMPONENTS = ['v', 'v35', 'v24']
CSR_DENOMINATORS = {'v': 'v', 'v35': 'v', 'v24': 'v'}
NOV15_COMPONENTS = ['v3', 'v8']
NOV15_CONSTANTS = {'v3': 1.0, 'v8': 3.0}
NOV15_DENOMINATORS = {'v3': 'v3', 'v8': 'v8'}

# The following lays out the SR for Polarised PDFs
TSR_COMPONENTS = ['t3', 't8']
TSR_DENOMINATORS = {'t3': 't3', 't8': 't8'}
# Sum Rules defined as in PDG 2023
TSR_CONSTANTS = {'t3': 1.2756, 't8': 0.5850}
TSR_CONSTANTS_UNC = {'t3': 0.0013, 't8': 0.025}


def sample_tsr(v: dict, e: dict, t: list, nr: int) -> list:
    """
    Sample the Triplets sum rules according to the PDG uncertainties.

    Parameters
    ----------
    v: dict
        dictionary that maps the triplet component to its PDG value
    e: dict
        dictionary that maps the triplet component to its denominator
    t: list
        list of triplet component for which SR should be applied
    nr: int
        number of replicas that are fitter simultaneously

    Returns
    -------
    list:
        list of sum rule values sampled according to a normal distribution
    """
    return [[np.random.normal(v[c], e[c]) for _ in range(nr)] for c in t]


class MSR_Normalization(MetaLayer):
    """
    Computes the normalisation factors for the sum rules of the PDFs.
    """

    _msr_enabled = False
    _vsr_enabled = False
    _tsr_enabled = False
    _csr_enabled = False

    def __init__(self, mode: str = "ALL", replicas: int = 1, **kwargs):
        if mode == True or mode.upper() == "ALL":
            self._msr_enabled = True
            self._vsr_enabled = True
        elif mode.upper() == "MSR":
            self._msr_enabled = True
        elif mode.upper() == "VSR":
            self._vsr_enabled = True
        elif mode.upper() == "TSR":
            self._tsr_enabled = True
        elif mode.upper() == "ALLBUTCSR":
            self._msr_enabled = True
            self._csr_enabled = True
        else:
            raise ValueError(f"Mode {mode} not accepted for sum rules")

        self.replicas = replicas
        indices = []
        self.divisor_indices = []
        if self._msr_enabled:
            indices += [IDX[c] for c in MSR_COMPONENTS]
            self.divisor_indices += [IDX[MSR_DENOMINATORS[c]] for c in MSR_COMPONENTS]
        if self._vsr_enabled:
            self.divisor_indices += [IDX[VSR_DENOMINATORS[c]] for c in VSR_COMPONENTS]
            indices += [IDX[c] for c in VSR_COMPONENTS]
            self.vsr_factors = op.numpy_to_tensor(
                [np.repeat(VSR_CONSTANTS[c], replicas) for c in VSR_COMPONENTS]
            )
        if self._tsr_enabled:
            self.divisor_indices += [IDX[TSR_DENOMINATORS[c]] for c in TSR_COMPONENTS]
            indices += [IDX[c] for c in TSR_COMPONENTS]
            self.tsr_factors = op.numpy_to_tensor(
                sample_tsr(TSR_CONSTANTS, TSR_CONSTANTS_UNC, TSR_COMPONENTS, replicas)
            )

        if self._csr_enabled:
            # modified vsr for V, V24, V35
            indices += [IDX[c] for c in CSR_COMPONENTS]
            self.divisor_indices += [IDX[CSR_DENOMINATORS[c]] for c in CSR_COMPONENTS]
            # no V15 vsr
            self.divisor_indices += [IDX[NOV15_DENOMINATORS[c]] for c in NOV15_COMPONENTS]
            indices += [IDX[c] for c in NOV15_COMPONENTS]
            self.vsr_factors = op.numpy_to_tensor(
                [np.repeat(NOV15_CONSTANTS[c], replicas) for c in NOV15_COMPONENTS]
            )
        # Need this extra dimension for the scatter_to_one operation
        self.indices = [[i] for i in indices]

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
        pdf_integrated: (Tensor(1, replicas, 14))
            the integrated PDF
        photon_integral: (Tensor(1, replicas, 1))
            the integrated photon PDF

        Returns
        -------
        normalization_factor: Tensor(replicas, 1, 14)
            The normalization factors per flavour.
        """
        # get rid of batch dimension and put replicas last
        reshape = lambda x: op.transpose(x[0])
        y = reshape(pdf_integrated)
        photon_integral = reshape(photon_integral)

        numerators = []

        if self._msr_enabled:
            numerators += [
                op.batchit(1.0 - y[IDX['sigma']] - photon_integral[0], batch_dimension=0)
            ]
        if self._vsr_enabled:
            numerators += [self.vsr_factors]
        if self._tsr_enabled:
            numerators += [self.tsr_factors]
        if self._csr_enabled:
            numerators += len(CSR_COMPONENTS) * [
                op.batchit(4.0 - 1.0 / 3.0 * y[IDX['v15']], batch_dimension=0)
            ]
            numerators += [self.vsr_factors]

        numerators = op.concatenate(numerators, axis=0)
        divisors = op.gather(y, self.divisor_indices, axis=0)

        # Fill in the rest of the flavours with 1
        num_flavours = y.shape[0]
        norm_constants = op.scatter_to_one(
            numerators / divisors, indices=self.indices, output_shape=(num_flavours, self.replicas)
        )

        return op.batchit(op.transpose(norm_constants), batch_dimension=1)
