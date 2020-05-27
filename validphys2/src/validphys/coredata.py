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
    """

    hadronic: bool
    Q0: float
    ndata: int
    xgrid: np.array
    sigma: pd.DataFrame
    metadata: dict = dataclasses.field(default_factory=dict, repr=False)

    # TODO: When we move to something other than the current fktable format,
    # we should apply the cuts directly before loading the file.
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
        if hasattr(cuts, 'load'):
            cuts = cuts.load()
        if cuts is None:
            return self
        newndata = len(cuts)
        newsigma = self.sigma.loc[cuts]
        return dataclasses.replace(self, ndata=newndata, sigma=newsigma)


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
class SystematicError:
    add: float
    mult: float
    sys_type: str #e.g ADD
    name: str #e.g UNCORR

    def __repr__(self):
        return (f"{self.__class__.__name__}(add={self.add}, mult={self.mult},"
                "sys_type={self.sys_type}, name={self.name})")


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
    commondata_table: pd.DataFrame
    systype_table: pd.DataFrame

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
        if hasattr(cuts, 'name') and self.setname != cuts.name:
            raise ValueError(f"The cuts provided are for {cuts.name} which does not apply "
                    f"to this commondata file: {self.setname}")

        if hasattr(cuts, 'load'):
            cuts = cuts.load()
        if cuts is None:
            return self

        # We must shift the cuts up by 1 since a cut of 0 implies the first data point
        # while commondata indexing starts at 1.
        cuts = list(map(lambda x: x + 1, cuts))

        newndata = len(cuts)
        new_commondata_table = self.commondata_table.loc[cuts]
        return dataclasses.replace(
            self, ndata=newndata, commondata_table=new_commondata_table
        )

    @property
    def central_values(self):
        return self.commondata_table["data"]

    @property
    def stat_errors(self):
        return self.commondata_table["stat"]

    @property
    def sys_errors(self):
        sys_table = self.commondata_table.drop(
            columns=["process", "kin1", "kin2", "kin3", "data", "stat"]
        )
        table = [
            [
                SystematicError(
                    add=sys_table[f"sys.add.{j}"][i],
                    mult=sys_table[f"sys.mult.{j}"][i],
                    sys_type=self.systype_table["type"][j],
                    name=self.systype_table["name"][j],
                )
                for j in self.systype_table.index
            ]
            for i in self.commondata_table.index
        ]
        return pd.DataFrame(
            table,
            columns=[f"sys.{i}" for i in self.systype_table.index],
            index=self.commondata_table.index,
        )

class DataSet:
    """
    Add a docstring
    """
    # XXX: Note the C++ class has 2 constructors which accepts
    # a boolean mask and thus handles cuts appropriately

    def __init__(self, cd, fkset, weight):
        self.cd = cd
        self.fkset = fkset
        self.weight = weight

        self.ndata = cd.ndata
        self.nsys = cd.nsys
        self.central_values = cd.central_values
        self.stat_errors = cd.stat_errors
        self.sys_errors = cd.sys_errors


    def __getitem__(self, mask):
        """ Use the __getitem__ dunder method to apply
        a cuts mask to the dataset object.

        Example
        -------
        """
        # TODO: add a check such that the mask is same length as
        # commondata
        pass

    # XXX: Maybe change this to a normal method since some computation occurs
    @property
    def covmat(self):
        systype = self.cd.systype_table["name"]
        skip_systematics = systype.where(systype=="SKIP").dropna().index
        sys_errors = self.sys_errors.drop(columns=[f"sys.{i}" for i in skip_systematics])

        def to_additive(sys_errors):
            result = []
            for idx, sys_error in enumerate(sys_errors):
                if sys_error.sys_type == "ADD":
                    result.append(sys_error.add)
                elif sys_error.sys_type == "MULT":
                    result.append(self.central_values.iloc[idx] * sys_error.mult / 100)
            return result

        additive_mat = sys_errors.apply(to_additive)

        # Matrix where all non-zero values are due to correlated systematics
        correlated_mat = additive_mat.copy()
        for i, j in enumerate(self.cd.systype_table["name"]):
            if j not in ("UNCORR", "THEORYUNCORR"):
                continue
            correlated_mat.iloc[:, i] = 0

        correlated_mat = correlated_mat.to_numpy()

        # Systematic uncertainties converted to additive uncertainties
        additive_mat = additive_mat.to_numpy()

        # Diagonal matrix containing the statistical uncertainty for each
        # data point
        stat_mat = np.diag(self.stat_errors.to_numpy())

        if not correlated_mat.size > 0:
            # Some datasets e.g CMSWCHARMRAT have zero systematic
            # uncertainties. We treat this special edge case
            # by creating a matrix of zeros so it can be broadcast
            # with the statistical error matrix at the end
            correlated_mat = additive_mat = np.zeros_like(stat_mat)

        # Avoid double counting in case the i = j element is a correlated
        # systematic we take care of this case in the einsum below
        off_diagonals = correlated_mat @ additive_mat.T
        np.fill_diagonal(off_diagonals, 0)

        cov_mat = (
            stat_mat ** 2
            + off_diagonals
            + np.diag(np.einsum("ij, ij -> i", additive_mat, additive_mat))
        )
        return cov_mat
