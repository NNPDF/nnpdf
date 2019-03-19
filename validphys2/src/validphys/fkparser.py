"""
fkparser.py

Parse FKtables into useful datastructures
"""
import io
import functools
import tarfile

import numpy as np


class BadFKTableError(Exception):
    pass

def load_fktable(spec):
    with open_fkpath(spec.fkpath) as handle:
        return parse_fktable(handle)


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
def _parse_fast_kernel(f):
    return np.loadtxt(f)


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


_KNOWN_SEGMENTS = {
    "GridDesc": _parse_string,
    "VersionInfo": _parse_fk_options,
    "GridInfo": functools.partial(
        _parse_fk_options, value_parsers={"HADRONIC": bool, "NDATA": int, "NX": int}
    ),
    "FlavourMap": _parse_flavour_map,
    "xGrid": _parse_xgrid,
}

def parse_fktable(f):
    line_and_stream = enumerate(f, start=1)
    res = {}
    lineno, header = next(line_and_stream)
    while True:
        marker, header_name = _parse_header(lineno, header)
        if header_name in _KNOWN_SEGMENTS:
            parser = _KNOWN_SEGMENTS[header_name]
        elif marker == b'{':
            parser = _parse_string
        elif marker == b'_':
            parser = _parse_fk_options
        else:
            raise RuntimeError("Should not be here")
        if header_name == 'FastKernel':
            res['FastKernel'] = _parse_fast_kernel(f)
            return res
        try:
            out, lineno, header = parser(line_and_stream)
        except Exception as e:
            #Note that the old lineno is the one we want
            raise BadFKTableError(f"Failed processing header {header_name} on line {lineno}") from e
        res[header_name] = out
        if not header:
            break
