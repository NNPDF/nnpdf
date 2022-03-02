"""
This module implements parsers for FKtable  and CFactor files into useful
datastructures, contained in the :py:mod:`validphys.coredata` module, which are
not backed by C++ managed memory, and so they can be easily pickled and
interfaces with common Python libraries.  The integration of these objects into
the codebase is currently work in progress, and at the moment this module
serves as a proof of concept.

Most users will be interested in using the high level interface
:py:func:`load_fktable`.  Given a :py:class:`validphys.core.FKTableSpec`
object, it returns an instance of :py:class:`validphys.coredata.FKTableData`,
an object with the required information to compute a convolution, with the
CFactors applied.

.. code-block:: python

    from validphys.fkparser import load_fktable
    from validphys.loader import Loader
    l = Loader()
    fk = l.check_fktable(setname="ATLASTTBARTOT", theoryID=53, cfac=('QCD',))
    res = load_fktable(fk)
"""
import io
import functools
import tarfile
import dataclasses
import logging

import numpy as np
import pandas as pd

from validphys.coredata import FKTableData, CFactorData
from reportengine.compat import yaml

log = logging.getLogger(__name__)


class BadCFactorError(Exception):
    """Exception raised when an CFactor cannot be parsed correctly"""


class BadFKTableError(Exception):
    """Exception raised when an FKTable cannot be parsed correctly"""


@dataclasses.dataclass(frozen=True)
class GridInfo:
    """Class containing the basic properties of an FKTable grid."""

    setname: str
    hadronic: bool
    ndata: int
    nx: int


def load_fktable(spec):
    """Load the data corresponding to a FKSpec object. The cfactors
    will be applied to the grid.
    If we have a new-type fktable, call directly `load()`, otherwise
    fallback to the old parser
    """
    if spec.legacy:
        with open_fkpath(spec.fkpath) as handle:
            tabledata = parse_fktable(handle)
    else:
        tabledata = spec.load()

    if not spec.cfactors:
        return tabledata

    ndata = tabledata.ndata
    cfprod = np.ones(ndata)
    for cf in spec.cfactors:
        with open(cf, "rb") as f:
            cfdata = parse_cfactor(f)
        if len(cfdata.central_value) != ndata:
            if tabledata.metadata.get("repetition_flag") or spec.norm:
                cfprod *= cfdata.central_value[0]
                continue
            raise BadCFactorError(
                "Length of cfactor data does not match the length of the fktable."
            )
        cfprod *= cfdata.central_value
    # TODO: Find a way to do this in place
    tabledata.sigma = tabledata.sigma.multiply(pd.Series(cfprod), axis=0, level=0)
    return tabledata


def _get_compressed_buffer(path):
    archive = tarfile.open(path)
    members = archive.getmembers()
    l = len(members)
    if l != 1:
        raise BadFKTableError(f"Archive {path} should contain one file, but it contains {l}.")
    return archive.extractfile(members[0])


def open_fkpath(path):
    """Return a file-like object from the fktable path, regardless of whether
    it is compressed

    Parameters
    ..........
    path: Path or str
        Path like file containing a valid FKTable. It can be either inside a
        tarball or in plain text.

    Returns
    -------
    f: file
        A file like object for further processing.
    """
    if tarfile.is_tarfile(path):
        return _get_compressed_buffer(path)
    return open(path, "rb")


def _is_header_line(line):
    return line.startswith((b"_", b"{"))


def _bytes_to_bool(x):
    return bool(int(x))


def _parse_fk_options(line_and_stream, value_parsers=None):
    """Parse a sequence of lines of the form
    *OPTION: VALUE
    into a dictionary.
    """
    res = {}
    if value_parsers is None:
        value_parsers = {}
    for lineno, next_line in line_and_stream:
        if _is_header_line(next_line):
            return res, lineno, next_line
        if not next_line.startswith(b"*"):
            raise BadFKTableError(f"Error on line {lineno}: Expecting an option starting with '*'")
        try:
            keybytes, valuebytes = next_line.split(b":", maxsplit=1)
        except ValueError:
            raise BadFKTableError(f"Error on line {lineno}: Expecting an option containing ':'")
        key = keybytes[1:].strip().decode()
        if key in value_parsers:
            try:
                value = value_parsers[key](valuebytes)
            except Exception as e:
                raise BadFKTableError(f"Could not parse key {key} on line {lineno}") from e
        else:
            value = valuebytes.strip().decode()
        res[key] = value

    raise BadFKTableError("FKTable should end with FastKernel spec, not with a set of options")


def _segment_parser(f):
    @functools.wraps(f)
    def f_(line_and_stream):
        buf = io.BytesIO()
        for lineno, next_line in line_and_stream:
            if _is_header_line(next_line):
                processed = f(buf)
                return processed, lineno, next_line
            buf.write(next_line)
        raise BadFKTableError("FKTable should end with FastKernel spec, not with a segment string")

    return f_


@_segment_parser
def _parse_string(buf):
    return buf.getvalue().decode()


@_segment_parser
def _parse_flavour_map(buf):
    buf.seek(0)
    return np.loadtxt(buf, dtype=bool)


@_segment_parser
def _parse_xgrid(buf):
    return np.fromstring(buf.getvalue(), sep="\n")


# This used a different interface from segment parser because we want it to
# be fast.
# We assume it is going to be the last section.
def _parse_hadronic_fast_kernel(f):
    """Parse the FastKernel secrion of an hadronic FKTable into a DataFrame.
    ``f`` should be a stream containing only the section"""
    # Note that we need the slower whitespace here because it turns out
    # that there are fktables where space and tab are used as separators
    # within the same table.
    df = pd.read_csv(f, sep=r"\s+", header=None, index_col=(0, 1, 2))
    df.columns = list(range(14 * 14))
    df.index.names = ["data", "x1", "x2"]
    return df


def _parse_dis_fast_kernel(f):
    """Parse the FastKernel section of a DIS FKTable into a DataFrame.
    ``f`` should be a stream containing only the section"""
    df = pd.read_csv(f, sep=r"\s+", header=None, index_col=(0, 1))
    df.columns = list(range(14))
    df.index.names = ["data", "x"]
    return df


def _parse_gridinfo(line_and_stream):
    dict_result, line_number, next_line = _parse_fk_options(
        line_and_stream, value_parsers={"HADRONIC": _bytes_to_bool, "NDATA": int, "NX": int}
    )
    gi = GridInfo(**{k.lower(): v for k, v in dict_result.items()})
    return gi, line_number, next_line


def _parse_header(lineno, header):
    if not _is_header_line(header):
        raise BadFKTableError(
            f"Bad header at line {lineno}: First " "character should be either '_' or '{'"
        )
    try:
        endname = header.index(b"_", 1)
    except ValueError:
        raise BadFKTableError(f"Bad header at line {lineno}: Expected '_' after name") from None
    header_name = header[1:endname]
    # Note: This is not the same as header[0]. Bytes iterate as ints.
    return header[0:1], header_name.decode()


def _build_sigma(f, res):
    gi = res["GridInfo"]
    fm = res["FlavourMap"]
    table = _parse_hadronic_fast_kernel(f) if gi.hadronic else _parse_dis_fast_kernel(f)
    # Filter out empty flavour indices
    table = table.loc[:, fm.ravel()]
    return table


_KNOWN_SEGMENTS = {
    "GridDesc": _parse_string,
    "VersionInfo": _parse_fk_options,
    "GridInfo": _parse_gridinfo,
    "FlavourMap": _parse_flavour_map,
    "xGrid": _parse_xgrid,
    "TheoryInfo": functools.partial(
        _parse_fk_options,
        value_parsers={
            "ID": int,
            "PTO": int,
            "DAMP": _bytes_to_bool,
            "IC": _bytes_to_bool,
            "XIR": float,
            "XIF": float,
            "NfFF": int,
            "MaxNfAs": int,
            "MaxNfPdf": int,
            "Q0": float,
            "alphas": float,
            "Qref": float,
            "QED": _bytes_to_bool,
            "alphaqed": float,
            "Qedref": float,
            "SxRes": _bytes_to_bool,
            "mc": float,
            "Qmc": float,
            "kcThr": float,
            "mb": float,
            "Qmb": float,
            "kbThr": float,
            "mt": float,
            "Qmt": float,
            "ktThr": float,
            "MZ": float,
            "MW": float,
            "GF": float,
            "SIN2TW": float,
            "TMC": _bytes_to_bool,
            "MP": float,
            "global_nx": int,
            "EScaleVar": _bytes_to_bool,
        },
    ),
}


def _check_required_sections(res, lineno):
    """Check that we have found all the required sections by the time we
    reach 'FastKernel'"""
    for section in _KNOWN_SEGMENTS:
        if section not in res:
            raise BadFKTableError(f"{section} must come before 'FastKernel' section at {lineno}")


def parse_fktable(f):
    """Parse an open byte stream into an FKTableData. Raise a BadFKTableError
    if problems are encountered.

    Parameters
    ----------
    f : file
        Open file-like object. See :func:`open_fkpath`to obtain it.

    Returns
    -------
    fktable : FKTableData
        An object containing the FKTable data and information.

    Notes
    -----
    This function operates at the level of a single file, and therefore it does
    not apply CFactors (see :py:func:`load_fktable` for that) or handle operations
    within COMPOUND ensembles.
    """
    line_and_stream = enumerate(f, start=1)
    res = {}
    lineno, header = next(line_and_stream)
    while True:
        marker, header_name = _parse_header(lineno, header)
        if header_name == "FastKernel":
            _check_required_sections(res, lineno)
            Q0 = res["TheoryInfo"]["Q0"]
            sigma = _build_sigma(f, res)
            hadronic = res["GridInfo"].hadronic
            ndata = res["GridInfo"].ndata
            xgrid = res.pop("xGrid")
            return FKTableData(
                sigma=sigma,
                ndata=ndata,
                Q0=Q0,
                metadata=res,
                hadronic=hadronic,
                xgrid=xgrid,
            )
        elif header_name in _KNOWN_SEGMENTS:
            parser = _KNOWN_SEGMENTS[header_name]
        elif marker == b"{":
            parser = _parse_string
        elif marker == b"_":
            parser = _parse_fk_options
        else:
            raise RuntimeError("Should not be here")
        try:
            out, lineno, header = parser(line_and_stream)
        except Exception as e:
            # Note that the old lineno is the one we want
            raise BadFKTableError(f"Failed processing header {header_name} on line {lineno}") from e
        res[header_name] = out


def parse_cfactor(f):
    """Parse an open byte stream into a :py:class`CFactorData`. Raise a
    BadCFactorError if problems are encountered.

    Parameters
    ----------
    f : file
        Binary file-like object

    Returns
    -------
    cfac : CFactorData
        An object containing the data on the cfactor for each point.
    """
    stars = f.readline()
    if not stars.startswith(b"*"):
        raise BadCFactorError("First line should start with '*'.")
    descbytes = io.BytesIO()
    for line in f:
        if line.startswith(b"*"):
            break
        descbytes.write(line)
    description = descbytes.getvalue().decode()
    try:
        data = np.loadtxt(f)
    except Exception as e:
        raise BadCFactorError(e) from e
    data = data.reshape(-1, 2)
    central_value = data[:, 0]
    uncertainty = data[:, 1]
    return CFactorData(
        description=description, central_value=central_value, uncertainty=uncertainty
    )


##### New fktable loader
EXT = "pineappl.lz4"


class PineAPPLEquivalentNotKnown(Exception):
    pass


class YamlFileNotFound(FileNotFoundError):
    pass


@dataclasses.dataclass
class FKTableComposer:
    """In the most general case, an FKTable for a given observable
    is formed by a number of pineappl grids that then need to be concatenated
    """

    grid_paths: list
    cfactor: str


def _load_yaml(yaml_file):
    """Load a dataset.yaml file"""
    if not yaml_file.exists():
        raise YamlFileNotFound(yaml_file)
    ret = yaml.safe_load(yaml_file.read_text())
    # Make sure the operations are upper-cased for compound-compatibility
    ret["operation"] = "NULL" if ret["operation"] is None else ret["operation"].upper()
    return ret


def get_yaml_information(yaml_file, theorypath, check_pineappl=False):
    """
    Given a yaml_file, return the corresponding dictionary
    with all information and an extra field "paths"
    with all the grids to be loaded for the given dataset.
    Checks whether the grid is apfelcomb or pineappl
    if check_pineappl is True this function will raise PineAPPLEquivalentNotKnown
    if a pineappl grid is not found.
    Parameters
    ----------
        yaml_file: Path
            path of the yaml file for the given dataset
        theorypath: Path
            path of the theory folder where to find the grids
    Returns
    -------
        yaml_content: dict
            Metadata prepared for the FKTables
        paths: list(list(path))
            List (of lists) with all the grids that will need to be loaded
    """
    yaml_content = _load_yaml(yaml_file)
    grids_folder = theorypath / "pineappls"

    if yaml_content.get("appl") and check_pineappl:
        # This might be useful to use the "legacy loader" when there is no actual pineappl available
        raise PineAPPLEquivalentNotKnown(yaml_content["target_dataset"])

    # Turn the operans and the members into paths (and check all of them exist)
    ret = []
    for operand in yaml_content["operands"]:
        tmp = []
        for member in operand:
            p = grids_folder / f"{member}.{EXT}"
            if not p.exists():
                raise FileNotFoundError(f"Failed to find {p}")
            tmp.append(p)
        ret.append(tmp)

    return yaml_content, ret


def pineappl_reader(fkspec):
    """
    Receives a fkspec, which contains the appropiate references to
    the pineappl grid to load and concatenate

    For more information: https://pineappl.readthedocs.io/en/latest/modules/pineappl/pineappl.html#pineappl.pineappl.PyFkTable
    """
    from pineappl.fk_table import FkTable

    pines = [FkTable.read(i) for i in fkspec.fkpath]

    # Use the first fktable to get some metadata from pineappl
    pp = pines[0]
    Q0 = np.sqrt(pp.muf2())
    xgrid = pp.x_grid()
    # Hadronic means in practice that not all luminosity combinations are just electron X proton
    non_hadronic = {11, 12, -11, -12}
    hadronic = all(non_hadronic.isdisjoint(i) for i in pp.lumi())
    # Now prepare the concatenation of grids
    if fkspec.norm:
        fktables = [np.sum(p.table(), axis=0, keepdims=True) for p in pines]
    else:
        fktables = [(p.table().T/p.bin_normalizations()).T for p in pines]
    fktable = np.concatenate(fktables, axis=0)

    # To create a dataframe compatible with that validphys creates from the old fktables we need to:
    # Step 1), make the luminosity into a 14x14 mask for the evolution basis
    # Step 2) prepare the indices for the dataframe
    eko_numbering_scheme = (22, 100, 21, 200, 203, 208, 215, 224, 235, 103, 108, 115, 124, 135)
    xi = np.arange(len(xgrid))
    # note that this is the same ordering that was used in fktables
    # the difference is that the fktables were doing this silently and now
    # we have the information made explicit
    if hadronic:
        flavour_map = np.zeros((14, 14), dtype=bool)
        for i, j in pp.lumi():
            idx = eko_numbering_scheme.index(i)
            jdx = eko_numbering_scheme.index(j)
            flavour_map[idx, jdx] = True

        co = np.where(flavour_map.ravel())[0]
        xdivision = (xgrid[:, None] * xgrid[None, :]).flatten()
    else:
        try:
            co = [eko_numbering_scheme.index(i) for _, i in pp.lumi()]
        except ValueError:
            co = [eko_numbering_scheme.index(i) for i, _ in pp.lumi()]
        xdivision = xgrid
    # The fktables for pineappl have an extra factor of x that we need to remove
    # hence xdivision

    # Step 3) Now put the flavours at the end and flatten
    # The output of pineappl is (ndata, flavours, x, x)
    lf = len(co)
    ndata = fktable.shape[0]
    xfktable = fktable.reshape(ndata, lf, -1) / xdivision
    fkmod = np.moveaxis(xfktable, 1, -1)

    # TODO: due to apfelcomb-pineappl incompatibilities
    # note: if it can be simplified it will be simplified but it should hopefully not get out
    #       of this function, hopefully fixing this should be done in the theory-pipeline
    #       before they get here

    # This needs to be always be applied first
    if fkspec.metadata.get("apfelcomb_norm"):
        total_norm = []
        for norms, operands in zip(
            fkspec.metadata.get("apfelcomb_norm"), fkspec.metadata["operands"]
        ):
            # Now check, in case of an operation, what is our index in this operation
            if len(operands) == len(fkspec.fkpath):
                for operand, fkpath, norm, fk in zip(operands, fkspec.fkpath, norms, fktables):
                    if fkpath.name == f"{operand}.{EXT}":
                        total_norm += [norm] * fk.shape[0]
        fkmod *= np.array(total_norm)[:, None, None]

    # we need to play some games with the dataframe
    # Look at the metadata to see whether we should apply a repetition cut here
    # Note: repetition always happen to the denominator when there's only 1
    if fkspec.metadata.get("repetition_flag"):
        valid_targets = []
        for operand, flag_state in zip(
            fkspec.metadata["operands"], fkspec.metadata.get("repetition_flag")
        ):
            if flag_state:
                valid_targets.append(f"{operand[0]}.{EXT}")
        # Now check whether the current fktable is part of the valid targets
        if fkspec.fkpath[0].name in valid_targets:
            fkmod = fkmod[0:1]
            ndata = 1

    ni = np.arange(ndata)
    if fkspec.metadata.get("shifts"):
        # Again, once we are here, anything different from this better fail
        shifts = fkspec.metadata.get("shifts")[0]
        ndata = 0
        ni = []
        for fk, shift in zip(fktables, shifts):
            if shift is not None:
                ndata += shift
            new_ndata = ndata + fk.shape[0]
            ni += list(range(ndata, new_ndata))
            ndata = new_ndata

    if hadronic:
        mi = pd.MultiIndex.from_product([ni, xi, xi], names=["data", "x1", "x2"])
    else:
        mi = pd.MultiIndex.from_product([ni, xi], names=["data", "x"])

    fkframe = fkmod.reshape(-1, lf)
    df = pd.DataFrame(fkframe, index=mi, columns=co)
    # Now prepare the FKTableData object
    # For the metadata use the one we have kept in the fkspec
    fkdata = FKTableData(
        sigma=df, ndata=ndata, Q0=Q0, metadata=fkspec.metadata, hadronic=hadronic, xgrid=xgrid
    )

    return fkdata
