"""
    n3fit interface to validphys

    Example
    -------

    >>> import numpy as np
    >>> from n3fit.vpinterface import N3PDF
    >>> from n3fit.model_gen import pdfNN_layer_generator
    >>> from validphys.pdfgrids import xplotting_grid
    >>> fake_fl = [{'fl' : i, 'largex' : [0,1], 'smallx': [1,2]} for i in ['u', 'ubar', 'd', 'dbar', 'c', 'cbar', 's', 'sbar']]
    >>> fake_x = np.linspace(1e-3,0.8,3)
    >>> pdf_model = pdfNN_layer_generator(nodes=[8], activations=['linear'], seed=0, flav_info=fake_fl)
    >>> n3pdf = N3PDF(pdf_model)
    >>> res = xplotting_grid(n3pdf, 1.6, fake_x)
    >>> res.grid_values.error_members().shape
    (1, 8, 3)


"""
import logging
from collections.abc import Iterable
import numpy as np
import numpy.linalg as la
from validphys.core import PDF, MCStats
from validphys.pdfbases import ALL_FLAVOURS, check_basis
from validphys.arclength import integrability_number, arc_lengths

log = logging.getLogger(__name__)
# Order of the evolution basis output from n3fit
EVOL_LIST = [
    "photon",
    "sigma",
    "gluon",
    "V",
    "V3",
    "V8",
    "V15",
    "V24",
    "V35",
    "T3",
    "T8",
    "T15",
    "T24",
    "T35",
]


class N3PDF(PDF):
    """
    Creates a N3PDF object, extension of the validphys PDF object to perform calculation
    with a n3fit generated model.

    Parameters
    ----------
        pdf_models: :py:class:`n3fit.backends.MetaModel` (or list thereof)
            PDF trained with n3fit, x -> f(x)_{i} where i are the flavours in the evol basis
        fit_basis: list(dict)
            basis of the training, used for reporting
        name: str
            name of the N3PDF object
    """

    def __init__(self, pdf_models, map_pdfs, fit_basis=None, name="n3fit"):
        self.map_pdfs = map_pdfs
        self.fit_basis = fit_basis
        self.basis = check_basis("evolution", EVOL_LIST)["basis"]

        if isinstance(pdf_models, Iterable):
            self._models = pdf_models
        else:
            self._models = [pdf_models]

        super().__init__(name)
        # The number of members will be minimum 2, needed for legacy compatibility
        # Note that the member 0 does not exist as an actual model as it corresponds
        # to the average of all others
        self._info = {"ErrorType": "replicas", "NumMembers": len(self._models) + 1}

    @property
    def stats_class(self):
        """The stats class of N3PDF is always MCStats"""
        return MCStats

    def load(self):
        """Many vp functions ask for a LHAPDF pdf
        from nnpdflib, this class fakes it until a time in which vp is free from C++
        """
        return self

    def get_member_model(self, replica):
        """Return the member corresponding to the given replica"""
        return self._models[replica - 1]

    def get_nn_weights(self):
        """Outputs all weights of the NN as numpy.ndarrays"""
        return [model.get_weights() for model in self._models]

    def get_preprocessing_factors(self, replica=None):
        """Loads the preprocessing alpha and beta arrays from the PDF trained model.
        If a ``fit_basis`` given in the format of ``n3fit`` runcards is given it will be used
        to generate a new dictionary with the names, the exponent and whether they are trainable
        otherwise outputs a Nx2 array where [:,0] are alphas and [:,1] betas
        """
        # If no replica is explicitly requested, get the preprocessing layer for the first model
        if replica is None:
            replica = 1
        preprocessing_layers = self.get_member_model(replica).get_layer_re(r"pdf_prepro_\d")
        if len(preprocessing_layers) != 1:
            # We really don't want to fail at this point, but print a warning at least...
            log.warning("More than one preprocessing layer found within the model!")
        preprocessing_layer = preprocessing_layers[0]

        alphas_and_betas = None
        if self.fit_basis is not None:
            output_dictionaries = []
            for idx in range(len(self.map_pdfs)):
                for d in self.fit_basis:
                    flavour = d["fl"] + "_idx_" + str(idx)
                    alpha = preprocessing_layer.get_weight_by_name(f"alpha_{flavour}")
                    beta = preprocessing_layer.get_weight_by_name(f"beta_{flavour}")
                    if alpha is not None:
                        alpha = float(alpha.numpy())
                    if beta is not None:
                        beta = float(beta.numpy())
                    output_dictionaries.append(
                        {
                            "fl": flavour,
                            "smallx": alpha,
                            "largex": beta,
                            "trainable": d.get("trainable", True),
                        }
                    )
            alphas_and_betas = output_dictionaries
        return alphas_and_betas

    def __call__(self, xarr, flavours=None, replica=None):
        """Uses the internal model to produce pdf values.
        The output is on the evolution basis.

        Parameters
        ----------
            xarr: numpy.ndarray
                x-points with shape (xgrid_size,) (size-1 dimensions are removed)
            flavours: list
                list of flavours to output
            replica: int
                replica whose value must be returned (by default return all members)
                replica 0 corresponds to the central value

        Returns
        -------
            numpy.ndarray
                (xgrid_size, flavours) pdf result
        """
        if flavours is None:
            flavours = EVOL_LIST
        # Ensures that the input has the shape the model expect, no matter the input
        # as the scaling is done by the model itself
        mod_xgrid = xarr.reshape(1, -1, 1)

        if replica is None or replica == 0:
            # We need generate output values for all replicas
            result = np.concatenate([m.predict([mod_xgrid]) for m in self._models], axis=0)
            if replica == 0:
                # We want _only_ the central value
                result = np.mean(result, axis=0, keepdims=True)
        else:
            result = self.get_member_model(replica).predict([mod_xgrid])

        # Split PDF according to the number of nuclei involved
        # TODO: Make the following split more consistent. The below might break!!
        split_size = int(result.shape[-1] / len(EVOL_LIST))
        pdf_per_amn = np.split(result, split_size, axis=-1)

        if flavours != "n3fit":
            # Ensure that the result has its flavour in the basis-defined order
            ii = self.basis._to_indexes(flavours)
            copy_pdf_per_amn = []
            for reslt in pdf_per_amn:
                reslt[:, :, ii] = reslt
                copy_pdf_per_amn.append(reslt)
            # Overwrite previous result
            pdf_per_amn = copy_pdf_per_amn
        return pdf_per_amn

    def grid_values(self, flavours, xarr, qmat=None):
        """
        Parameters
        ----------
            flavours: numpy.ndarray
                flavours to compute
            xarr: numpy.ndarray
                x-points to compute, dim: (xgrid_size,)
            qmat: numpy.ndarray
                q-points to compute (not used by n3fit, used only for shaping purposes)

        Returns
        ------
            numpy.ndarray
            array of shape (replicas, flavours, xgrid_size, qmat) with the values of
                the ``pdf_model``(s) evaluated in ``xarr``
        """
        n3fit_result = self(xarr.reshape(1, -1, 1))

        # The results of n3fit are always in the 14-evolution basis used in fktables
        # the calls to grid_values always assume the result will be LHAPDF flavours
        # we need then to rotate them to the LHAPDF-flavour basis,
        # we don't care that much for precision here
        to_flav = la.inv(self.basis.from_flavour_mat)
        to_flav[np.abs(to_flav) < 1e-12] = 0.0
        flav_result = np.tensordot(n3fit_result, to_flav, axes=(-1, 1))
        # Now drop the indices that are not requested
        requested_flavours = [ALL_FLAVOURS.index(i) for i in flavours]
        flav_result = flav_result[..., requested_flavours]

        # Swap the flavour and xgrid axis for vp-compatibility
        ret = flav_result.swapaxes(-2, -1)
        # If given, insert as many copies of the grid as q values
        ret = np.expand_dims(ret, -1)
        if qmat is not None:
            lq = len(qmat)
            ret = ret.repeat(lq, -1)
        return ret

    # # Utilities
    # def integrability_numbers(self, q0=1.65, flavours=None):
    #     """Compute the integrability numbers for the current PDF
    #     using the corresponding validphys action

    #     Parameters
    #     ----------
    #         q0: float
    #             energy at which the integrability is computed
    #         flavours: list
    #             flavours for which the integrability is computed

    #     Returns
    #     -------
    #         np.array(float)
    #             Value for the integrability for each of the flavours

    #     Example
    #     -------
    #     >>> from n3fit.vpinterface import N3PDF
    #     >>> from n3fit.model_gen import pdfNN_layer_generator
    #     >>> fake_fl = [
    #         {'fl' : i, 'largex' : [0,1], 'smallx': [1,2]} for i in ['u', 'ubar', 'd', 'dbar', 'c', 'cbar', 's', 'sbar']
    #         ]
    #     >>> pdf_model = pdfNN_layer_generator(nodes=[8], activations=['linear'], seed=0, flav_info=fake_fl)
    #     >>> n3pdf = N3PDF(pdf_model)
    #     >>> res = n3pdf.integrability_numbers()
    #     """
    #     if flavours is None:
    #         flavours = ["V", "T3", "V3", "T8", "V8"]
    #     return integrability_number(self, [q0], flavours=flavours)

    # def compute_arclength(self, q0=1.65, basis="evolution", flavours=None):
    #     """
    #     Given the layer with the fit basis computes the arc length
    #     using the corresponding validphys action

    #     Parameters
    #     ----------
    #         pdf_function: function
    #             pdf function has received by the writer or ``pdf_model``
    #         q0: float
    #             energy at which the arc length is computed
    #         basis: str
    #             basis in which to compute the arc length
    #         flavours: list
    #             output flavours

    #     Example
    #     -------
    #     >>> from n3fit.vpinterface import N3PDF
    #     >>> from n3fit.model_gen import pdfNN_layer_generator
    #     >>> fake_fl = [
    #         {'fl' : i, 'largex' : [0,1], 'smallx': [1,2]} for i in ['u', 'ubar', 'd', 'dbar', 'c', 'cbar', 's', 'sbar']
    #     ]
    #     >>> pdf_model = pdfNN_layer_generator(nodes=[8], activations=['linear'], seed=0, flav_info=fake_fl)
    #     >>> n3pdf = N3PDF(pdf_model)
    #     >>> res = n3pdf.compute_arclength()
    #     """
    #     if flavours is None:
    #         flavours = ["sigma", "gluon", "V", "V3", "V8"]
    #     ret = arc_lengths(self, [q0], basis, flavours)
    #     return ret.stats.central_value()
