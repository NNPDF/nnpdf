"""
Data containers backed by Python managed memory (Numpy arrays and Pandas
dataframes).
"""

import dataclasses
import logging

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


@dataclasses.dataclass(eq=False)
class FKTableData:
    """
    Data contained in an FKTable

    Parameters
    ----------
    hadronic : bool
        Whether a hadronic (two PDFs) or a DIS (one PDF) convolution is needed.

    Q0 : float
        The scale at which the PDFs should be evaluated (in GeV).

    ndata : int
        The number of data points in the grid.

    xgrid : array, shape (nx)
        The points in x at which the PDFs should be evaluated.

    sigma : pd.DataFrame
        For hadronic data, the columns are the indexes in the ``NfxNf`` list of
        possible flavour combinations of two PDFs.  The MultiIndex contains
        three keys, the data index, an index into ``xgrid`` for the first PDF
        and an idex into ``xgrid`` for the second PDF, indicating if the points in
        ``x`` where the PDF should be evaluated.

        For DIS data, the columns are indexes in the ``Nf`` list of flavours.
        The MultiIndex contains two keys, the data index and an index into
        ``xgrid`` indicating the points in ``x`` where the PDF should be
        evaluated.

    convolution_types: list(pineappl.convolutions.Conv)
        The type of convolution that the FkTable is expecting for each of the
        functions to be convolved with (usually the two types of PDF from the two
        incoming hadrons).

    metadata : dict
        Other information contained in the FKTable.

    protected: bool
        When a fktable is protected cuts will not be applied.
        The most common use-case is when a total cross section is used
        as a normalization table for a differential cross section,
        in legacy code (<= NNPDF4.0) both fktables would be cut using the differential index.
    """

    hadronic: bool
    Q0: float
    ndata: int
    xgrid: np.ndarray
    sigma: pd.DataFrame
    convolution_types: list | None = None
    metadata: dict = dataclasses.field(default_factory=dict, repr=False)
    protected: bool = False

    def with_cfactor(self, cfactor):
        """Returns a copy of the FKTableData object with cfactors applied to the fktable"""
        if all(c == 1.0 for c in cfactor):
            return self
        if len(cfactor) != self.ndata:
            if self.protected:
                cfactor = cfactor[0]
            else:
                name = self.metadata.get("target_dataset")
                raise ValueError(
                    f"The length of cfactor for {name} differs from the number of datapoints in the grid"
                )
        new_sigma = self.sigma.multiply(pd.Series(cfactor), axis=0, level=0)
        return dataclasses.replace(self, sigma=new_sigma)

    def with_cuts(self, cuts):
        """Return a copy of the FKTable with the cuts applied.  The data index
        of the sigma operator (the outermost level), contains the data point
        that have been kept. The ndata property is updated to reflect the new
        number of datapoints. If cuts is None, return the object unmodified.

        Parameters
        ----------
        cuts : array_like or validphys.core.Cuts or None.
            The cuts to be applied.

        Returns
        -------
        res : FKTableData
            A copy of the FKtable with the cuts applies.

        Notes
        -----
        The original number of points can be accessed with
        ``table.metadata['GridInfo'].ndata``.

        Examples
        --------

        >>> from validphys.fkparser import load_fktable
        ... from validphys.loader import Loader
        ... l = Loader()
        ... ds = l.check_dataset('ATLASTTBARTOT', theoryid=53, cfac=('QCD',))
        ... table = load_fktable(ds.fkspecs[0])
        ... newtable = table.with_cuts([0,1])
        >>> assert set(newtable.sigma.index.get_level_values(0)) == {0,1}
        >>> assert newtable.ndata == 2
        >>> assert newtable.metadata['GridInfo'].ndata == 3
        """
        if hasattr(cuts, "load"):
            cuts = cuts.load()
        if cuts is None or self.protected:
            return self
        newndata = len(cuts)
        try:
            newsigma = self.sigma.loc[cuts]
        except KeyError as e:
            # This will be an ugly erorr msg, but it should be scary anyway
            log.error(f"Problem applying cuts to {self.metadata}")
            raise e
        return dataclasses.replace(self, ndata=newndata, sigma=newsigma)

    @property
    def luminosity_mapping(self):
        """Return the flavour combinations that contribute to the fktable
        in the form of a single array

        The return shape is:
            (nbasis,) for DIS
            (nbasis*2,) for hadronic
        """
        basis = self.sigma.columns.to_numpy()
        if self.hadronic:
            ret = np.zeros(14 * 14, dtype=bool)
            ret[basis] = True
            basis = np.array(np.where(ret.reshape(14, 14))).T.reshape(-1)
        return basis

    def get_np_fktable(self):
        """Returns the fktable as a dense numpy array that can be directly
        manipulated with numpy

        The return shape is:
            (ndata, nx, nbasis) for DIS
            (ndata, nx, nx, nbasis) for hadronic
        where nx is the length of the xgrid
        and nbasis the number of flavour contributions that contribute
        """
        # Read up the shape of the output table
        ndata = self.ndata
        nx = len(self.xgrid)
        nbasis = self.sigma.shape[1]

        if ndata == 0:
            if self.hadronic:
                return np.zeros((ndata, nbasis, nx, nx))
            return np.zeros((ndata, nbasis, nx))

        # Make the dataframe into a dense numpy array

        # First get the data index out of the way
        # this is necessary because cuts/shifts and for performance reasons
        # otherwise we will be putting things in a numpy array in very awkward orders
        ns = self.sigma.unstack(level=("data",), fill_value=0)
        x1 = ns.index.get_level_values(0)

        if self.hadronic:
            x2 = ns.index.get_level_values(1)
            fk_raw = np.zeros((nx, nx, ns.shape[1]))
            fk_raw[x2, x1, :] = ns.values

            # The output is (ndata, basis, x1, x2)
            fktable = fk_raw.reshape((nx, nx, nbasis, ndata)).T
        else:
            fk_raw = np.zeros((nx, ns.shape[1]))
            fk_raw[x1, :] = ns.values

            # The output is (ndata, basis, x1)
            fktable = fk_raw.reshape((nx, nbasis, ndata)).T

        return fktable

    def determine_pdfs(self, pdf):
        """Determine the PDF (or PDFs) that should be used to be convoluted with this fktable.
        Uses the `convolution_types` key to decide the PDFs.
        If `convolution_types` is not defined, it returns the pdf object.
        """
        if self.convolution_types is None:
            if self.hadronic:
                return [pdf, pdf]
            return [pdf]

        conv_pdfs = []
        for convolution_type in self.convolution_types:
            # Check the polarization of the current convolution
            polarized = convolution_type.conv_type.polarized
            # Check the type of convolutions that the fktable is asking for and match it to the PDF
            if not polarized:
                if pdf.is_polarized:
                    if pdf.unpolarized_bc is None:
                        raise ValueError(
                            "The FKTable asked for an unpolarized PDF but received only polarized PDFs"
                        )
                    conv_pdfs.append(pdf.unpolarized_bc.make_only_cv())
                else:
                    conv_pdfs.append(pdf)
            elif polarized:
                if not pdf.is_polarized:
                    raise ValueError(
                        """The FKTable asked for a polarized PDF, but the PDF received cannot be understood
                        as polarized. When using a polarized PDF make sure to include a boundary condition
                        `unpolarized_bc: <pdf name>` whenever needed (`t0`, `dataspecs`...)."""
                    )
                conv_pdfs.append(pdf)
            else:  # Other scenarios (such as `time_like`) should be implemented as another `elif` statement
                raise ValueError("The convolution type is not recognized!")

        return conv_pdfs


@dataclasses.dataclass(eq=False)
class CFactorData:
    """
    Data contained in a CFactor

    Parameters
    ----------

    description : str
        Information on how the data was obtained.

    central_value : array, shape(ndata)
        The value of the cfactor for each data point.

    uncertainty : array, shape(ndata)
        The absolute uncertainty on the cfactor if available.
    """

    description: str
    central_value: np.array
    uncertainty: np.array
