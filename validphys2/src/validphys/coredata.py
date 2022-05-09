"""
Data containers backed by Python managed memory (Numpy arrays and Pandas
dataframes).  This module is intended to substitute large parts of the C++
wrappers.

"""
import dataclasses
import numpy as np
import pandas as pd


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
                raise ValueError(f"The length of cfactor for {name} differs from the number of datapoints in the grid")
        new_sigma = self.sigma.multiply(pd.Series(cfactor), axis=0, level=0)
        return dataclasses.replace(self, sigma=new_sigma)

    def with_cuts(self, cuts):
        """Return a copy of the FKTabe with the cuts applied.  The data index
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
        newsigma = self.sigma.loc[cuts]
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


@dataclasses.dataclass(eq=False)
class CommonData:
    """
    Data contained in Commondata files, relevant cuts applied.

    Parameters
    ----------

    setname : str
        Name of the dataset

    ndata : int
        Number of data points

    commondataproc : str
        Process type, one of 21 options

    nkin : int
        Number of kinematics specified

    kinematics : list of str with length nkin
        Kinematic variables kin1, kin2, kin3 ...

    nsys : int
        Number of systematics

    sysid : list of str with length nsys
        ID for systematic

    commondata_table : pd.DataFrame
        Pandas dataframe containing the commondata

    systype_table : pd.DataFrame
        Pandas dataframe containing the systype index
        for each systematic alongside the uncertainty
        type (ADD/MULT/RAND) and name
        (CORR/UNCORR/THEORYCORR/SKIP)
    """

    setname: str
    ndata: int
    commondataproc: str
    nkin: int
    nsys: int
    commondata_table: pd.DataFrame = dataclasses.field(repr=False)
    systype_table: pd.DataFrame = dataclasses.field(repr=False)
    systematics_table: pd.DataFrame = dataclasses.field(init=None, repr=False)

    def __post_init__(self):
        self.systematics_table = self.commondata_table.drop(
            columns=["process", "kin1", "kin2", "kin3", "data", "stat"]
        )

    def with_cuts(self, cuts):
        """A method to return a CommonData object where
        an integer mask has been applied, keeping only data
        points which pass cuts.

        Note if the first data point passes cuts, the first entry
        of ``cuts`` should be ``0``.

        Paramters
        ---------
        cuts: list or validphys.core.Cuts or None
        """
        # Ensure that the cuts we're applying applies to this dataset
        # only check, however, if the cuts is of type :py:class:`validphys.core.Cuts`
        if hasattr(cuts, "name") and self.setname != cuts.name:
            raise ValueError(
                f"The cuts provided are for {cuts.name} which does not apply "
                f"to this commondata file: {self.setname}"
            )

        if hasattr(cuts, "load"):
            cuts = cuts.load()
        if cuts is None:
            return self

        # We must shift the cuts up by 1 since a cut of 0 implies the first data point
        # while commondata indexing starts at 1.
        cuts = list(map(lambda x: x + 1, cuts))

        newndata = len(cuts)
        new_commondata_table = self.commondata_table.loc[cuts]
        return dataclasses.replace(self, ndata=newndata, commondata_table=new_commondata_table)

    @property
    def central_values(self):
        return self.commondata_table["data"]

    @property
    def stat_errors(self):
        return self.commondata_table["stat"]

    @property
    def multiplicative_errors(self):
        """Returns the systematics which are multiplicative (systype is MULT)
        in a percentage format, with SKIP uncertainties removed.

        """
        mult_systype = self.systype_table[self.systype_table["type"] == "MULT"]
        # NOTE: Index with list here so that return is always a DataFrame, even
        # if N_sys = 1 (else a Series could be returned)
        mult_table = self.systematics_table.loc[:, ["MULT"]]
        # Minus 1 because iloc starts from 0, while the systype counting starts
        # from 1.
        mult_table = mult_table.iloc[:, mult_systype.index - 1]
        mult_table.columns = mult_systype["name"].to_numpy()
        return mult_table.loc[:, mult_table.columns != "SKIP"]

    @property
    def additive_errors(self):
        """Returns the systematics which are additive (systype is ADD) as
        absolute uncertainties (same units as data), with SKIP uncertainties
        removed.

        """
        add_systype = self.systype_table[self.systype_table["type"] == "ADD"]
        # NOTE: Index with list here so that return is always a DataFrame, even
        # if N_sys = 1 (else a Series could be returned)
        add_table = self.systematics_table.loc[:, ["ADD"]]
        # Minus 1 because iloc starts from 0, while the systype counting starts
        # from 1.
        add_table = add_table.iloc[:, add_systype.index - 1]
        add_table.columns = add_systype["name"].to_numpy()
        return add_table.loc[:, add_table.columns != "SKIP"]

    def systematic_errors(self, central_values=None):
        """Returns all systematic errors as absolute uncertainties, with a
        single column for each uncertainty. Converts
        :py:attr:`multiplicative_errors` to units of data and then appends
        onto :py:attr:`additive_errors`. By default uses the experimental
        central values to perform conversion, but the user can supply a
        1-D array of central values, with length :py:attr:`self.ndata`, to use
        instead of the experimental central values to calculate the absolute
        contribution of the multiplicative systematics.

        Parameters
        ----------
        central_values: None, np.array
            1-D array containing alternative central values to combine with
            multiplicative uncertainties. This array must have length equal
            to :py:attr:`self.ndata`. By default ``central_values`` is None, and
            the central values of the commondata are used.

        Returns
        -------
        systematic_errors: pd.DataFrame
            Dataframe containing systematic errors.

        """
        if central_values is None:
            central_values = self.central_values.to_numpy()
        converted_mult_errors = self.multiplicative_errors * central_values[:, np.newaxis] / 100
        return pd.concat((self.additive_errors, converted_mult_errors), axis=1)
