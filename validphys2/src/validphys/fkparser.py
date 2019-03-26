"""
fkparser.py

Parse FKtables and CFactors into useful datastructures.
"""
import io
import functools
import tarfile
import dataclasses

import numpy as np
import pandas as pd


@dataclasses.dataclass(frozen=True)
class GridInfo:
    setname: str
    hadronic: bool
    ndata: int
    nx: int

@dataclasses.dataclass(eq=False)
class FKTableData:
    hadronic: bool
    Q0: float
    ndata: int
    xgrid: np.array
    sigma: pd.DataFrame
    metadata: dict = dataclasses.field(default_factory=dict, repr=False)

@dataclasses.dataclass(eq=False)
class CFactorData:
    description: str
    central_value: np.array
    uncertainty: np.array

class BadCFactorError(Exception):
    """Exception raised when an CFactor cannot be parsed correctly"""


class BadFKTableError(Exception):
    """Exception raised when an FKTable cannot be parsed correctly"""

def load_fktable(spec):
    """Load the data corresponding to a FKSpec object. The cfactors
    will be applied to the grid."""
    with open_fkpath(spec.fkpath) as handle:
        tabledata = parse_fktable(handle)
    if not spec.cfactors:
        return tabledata

    ndata = tabledata.ndata
    cfprod = np.ones(ndata)
    for cf in spec.cfactors:
        with open(cf, "rb") as f:
            cfdata = parse_cfactor(f)
        if len(cfdata.central_value) != ndata:
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
        raise BadFKTableError(
            f"Archive {path} should contain one file, but it contains {l}.")
    return archive.extractfile(members[0])


def open_fkpath(path):
    """Return a file-like object from the fktable path, regardless of whether
    it is compressed"""
    if tarfile.is_tarfile(path):
        return _get_compressed_buffer(path)
    return open(path, 'rb')


def _is_header_line(line):
    return line.startswith((b'_', b'{'))

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
        if not next_line.startswith(b'*'):
            raise BadFKTableError(f"Error on line {lineno}: Expecting an option starting with '*'")
        try:
            keybytes, valuebytes = next_line.split(b':', maxsplit=1)
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
    return np.fromstring(buf.getvalue(), sep='\n')

# This used a differen interface from segment parser because we want it to
# be fast.
# We assume it is going to be the last section.
def _parse_hadronic_fast_kernel(f):
    """Parse the FastKernel secrion of an hadronic FKTable into a DataFrame.
    ``f`` should be a stream containing only the section"""
    # Note that we need the slower whitespace here because it turns out
    # that there are fktables where space and tab are used as separators
    # within the same table.
    df = pd.read_csv(f, sep=r'\s+', header=None, index_col=(0,1,2))
    df.columns = list(range(14*14))
    df.index.names = ['data', 'x1', 'x2']
    return df

def _parse_dis_fast_kernel(f):
    """Parse the FastKernel section of a DIS FKTable into a DataFrame.
    ``f`` should be a stream containing only the section"""
    df = pd.read_csv(f, sep=r'\s+', header=None, index_col=(0,1))
    df.columns = list(range(14))
    df.index.names = ['data', 'x']
    return df


def _parse_gridinfo(line_and_stream):
    d, l, h = _parse_fk_options(
        line_and_stream,
        value_parsers={
            "HADRONIC": _bytes_to_bool,
            "NDATA": int,
            "NX": int
        })
    gi = GridInfo(**{k.lower(): v for k, v in d.items()})
    return gi, l, h



def _parse_header(lineno, header):
    if not _is_header_line(header):
        raise BadFKTableError(f"Bad header at line {lineno}: First "
                "character should be either '_' or '{'")
    try:
        endname = header.index(b'_', 1)
    except ValueError:
        raise BadFKTableError(f"Bad header at line {lineno}: Expected '_' after name") from None
    header_name = header[1:endname]
    #Note: This is not the same as header[0]. Bytes iterate as ints.
    return header[0:1], header_name.decode()


def _build_sigma(lineno, f, res):
    gi = res["GridInfo"]
    fm = res["FlavourMap"]
    table = (
        _parse_hadronic_fast_kernel(f) if gi.hadronic else _parse_dis_fast_kernel(f)
    )
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
            raise BadFKTableError(
                f"{section} must come before 'FastKernel' section at {lineno}"
            )

def parse_fktable(f):
    line_and_stream = enumerate(f, start=1)
    res = {}
    lineno, header = next(line_and_stream)
    while True:
        marker, header_name = _parse_header(lineno, header)
        if header_name == "FastKernel":
            _check_required_sections(res, lineno)
            Q0 = res['TheoryInfo']['Q0']
            sigma = _build_sigma(lineno, f, res)
            hadronic = res['GridInfo'].hadronic
            ndata = res['GridInfo'].ndata
            xgrid = res.pop('xGrid')
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
        elif marker == b'{':
            parser = _parse_string
        elif marker == b'_':
            parser = _parse_fk_options
        else:
            raise RuntimeError("Should not be here")
        try:
            out, lineno, header = parser(line_and_stream)
        except Exception as e:
            # Note that the old lineno is the one we want
            raise BadFKTableError(
                f"Failed processing header {header_name} on line {lineno}"
            ) from e
        res[header_name] = out


def parse_cfactor(f):
    """Parse an open byte stream into a ``CFactorData"""
    stars = f.readline()
    if not stars.startswith(b'*'):
        raise BadCFactorError("First line should start with '*'.")
    descbytes = io.BytesIO()
    for line in f:
        if line.startswith(b'*'):
            break
        descbytes.write(line)
    description = descbytes.getvalue().decode()
    try:
        data = np.loadtxt(f)
    except Exception as e:
        raise BadCFactorError(e) from e
    central_value = data[:, 0]
    uncertainty = data[:, 1]
    return CFactorData(
        description=description, central_value=central_value, uncertainty=uncertainty
    )
