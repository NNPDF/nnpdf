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
    >>> res.grid_values.shape
    (1, 8, 3)


"""
import numpy as np
import numpy.linalg as la
from validphys.core import PDF, MCStats
from validphys.pdfbases import ALL_FLAVOURS, evolution


class N3PDF(PDF):
    def __init__(self, pdf, name="n3fit"):
        self.model = pdf
        # Set the number of members to two for legacy compatibility
        # in this case replica 0 and replica 1 are the same
        self.NumMembers = 2
        super().__init__(name)

    @property
    def stats_class(self):
        """ The stats class of N3PDF is always be MCStats """
        return MCStats

    def load(self):
        """ Many vp functions ask for a LHAPDF pdf
        from nnpdflib, this class fakes it until a time in which vp is free from C++
        """
        return self

    def predict(self, xarr):
        """ Uses the internal model to produce pdf values
        Parameters
        ----------
            xarr: numpy.array
                x-points with shape (1, xgrid_size, 1)
        """
        if hasattr(self.model, 'predict'):
            return self.model.predict([xarr])
        else:
            ret = self.model(xarr.squeeze(0))
            # add a batch dimension for n3fit-compatibility
            return np.expand_dims(ret, 0)

    def grid_values(self, flavs, xarr, qmat=None):
        """
        Parameters
        ----------
            flavs: numpy.array
                flavours to compute
            xarr: numpy.array
                x-points to compute, dim: (xgrid_size,)
            qmat: numpy.array
                q-points to compute (not used by n3fit, used only for shaping purposes)

        Returns
        ------
            numpy.array
            array of shape (1, flavs, xgrid_size, qmat) with the values of the ``pdf_model``
            evaluated in ``xarr``
        """
        evol_result = self.predict(xarr.reshape(1, -1, 1))
        # The results of n3fit are always in the 14-evolution basis used in fktables
        # we need to rotate them to the LHAPDF-flavour basis,
        # we don't care that much for precision here
        evol2flav = la.inv(evolution.from_flavour_mat)
        evol2flav[np.abs(evol2flav) < 1e-12] = 0.0
        flav_result = np.tensordot(evol_result, evol2flav, axes=[-1, 1])
        # Now drop the indices that are not requested
        requested_flavours = [ALL_FLAVOURS.index(i) for i in flavs]
        flav_result = flav_result[:, :, requested_flavours]
        # Swap the flavour and xgrid axis for vp-compatibility
        ret = flav_result.swapaxes(1, 2)
        # If given, insert as many copies of the grid as q values
        ret = np.expand_dims(ret, -1)
        if qmat is not None:
            lq = len(qmat)
            ret = ret.repeat(lq, -1)
        return ret
