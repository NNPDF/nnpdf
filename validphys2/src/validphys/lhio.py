"""
A module that reads and writes LHAPDF grids.
"""

import logging
import os
import os.path as osp
import pathlib
import shutil

import lhapdf
import numpy as np
import pandas as pd

from reportengine.compat import yaml
from validphys import lhaindex
from validphys.core import PDF

log = logging.getLogger(__name__)

def split_sep(f):
    for line in f:
        if line.startswith(b'---'):
            break
        yield line

def read_xqf_from_file(f):

    lines = split_sep(f)
    try:
        (xtext, qtext, ftext) = [next(lines) for _ in range(3)]
    except StopIteration:
        return None
    xvals = np.fromstring(xtext, sep = " ")
    qvals = np.fromstring(qtext, sep = " ")
    fvals = np.fromstring(ftext, sep = " ", dtype=int)
    vals = np.fromstring(b''.join(lines), sep= " ")
    return pd.Series(vals, index = pd.MultiIndex.from_product((xvals, qvals, fvals)))


def read_xqf_from_lhapdf(pdf, replica, kin_grids):
    indexes = tuple(kin_grids.index)
    #Use LHAPDF directly to avoid the insanely deranged replica 0 convention
    #of libnnpdf.
    #TODO: Find a way around this
    import lhapdf

    xfxQ = lhapdf.mkPDF(pdf.name, int(replica)).xfxQ

    vals = []
    for x in indexes:
        #TODO: Change this for a faster grid_values call
        vals += [xfxQ(x[3],x[1],x[2])]
    return pd.Series(vals, index = kin_grids.index)

def read_all_xqf(f):
    while True:
        result = read_xqf_from_file(f)
        if result is None:
            return
        yield result

def load_replica(pdf, rep, kin_grids=None):

    suffix = str(rep).zfill(4)

    pdf_name = str(pdf)

    path = osp.join(lhaindex.finddir(pdf_name),
                    pdf_name + "_" + suffix + ".dat")

    log.debug("Loading replica {rep} at {path}".format(rep=rep,
                                                       path=path))

    with open(path, 'rb') as inn:
        header = b"".join(split_sep(inn))

        if kin_grids is not None:
            xfqs = read_xqf_from_lhapdf(pdf, rep, kin_grids)
        else:
            xfqs = list(read_all_xqf(inn))
            xfqs = pd.concat(xfqs, keys=range(len(xfqs)))
    return header, xfqs

#Split this to debug easily
def _rep_to_buffer(out, header, subgrids):
    sep = b'---'
    out.write(header)
    out.write(sep)
    for _,g in subgrids.groupby(level=0):
        out.write(b'\n')
        ind = g.index.get_level_values(1).unique()
        np.savetxt(out, ind, fmt='%.7E',delimiter=' ', newline=' ')
        out.write(b'\n')
        ind = g.index.get_level_values(2).unique()
        np.savetxt(out, ind, fmt='%.7E',delimiter=' ', newline=' ')
        out.write(b'\n')
        #Integer format
        ind = g.index.get_level_values(3).unique()
        np.savetxt(out, ind, delimiter=' ', fmt="%d",
                      newline=' ')
        out.write(b'\n ')
        #Reshape so printing is easy
        reshaped = g.values.reshape((len(g.groupby(level=1))*len(g.groupby(level=2)),
                              len(g.groupby(level=3))))
        np.savetxt(out, reshaped, delimiter=" ", newline="\n", fmt='%14.7E')
        out.write(sep)

def write_replica(rep, set_root, header, subgrids):
    suffix = str(rep).zfill(4)
    target_file = set_root / f'{set_root.name}_{suffix}.dat'
    if target_file.is_file():
        log.warning(f"Overwriting replica file {target_file}")
    with open(target_file, 'wb') as out:
        _rep_to_buffer(out, header, subgrids)

def load_all_replicas(pdf, db=None):
    if db is not None:
        #removing str() will crash as it casts to unicode due to pdf name
        key = str("(load_all_replicas, %s)" % pdf.get_key())
        if key in db:
            return db[key]
    rep0headers, rep0grids = load_replica(pdf, 0)

    headers, grids = zip(*[load_replica(pdf, rep, rep0grids)
                         for rep in range(1, len(pdf))])
    result = [rep0headers] + list(headers), [rep0grids] + list(grids)
    if db is not None:
        db[key] = result
    return result

def big_matrix(gridlist):
    """Return a properly indexes matrix of the differences between each member
    and the central value"""
    central_value = gridlist[0]
    X = pd.concat(gridlist[1:], axis=1,
                 keys=range(1,len(gridlist)+1), #avoid confusion with rep0
                 ).subtract(central_value, axis=0)
    if np.any(X.isnull()) or X.shape[0] != len(central_value):
        raise ValueError("Incompatible grid specifications")
    return X

def rep_matrix(gridlist):
    """Return a properly indexes matrix of all the members"""
    X = pd.concat(gridlist, axis=1,
                 keys=range(1,len(gridlist)+1), #avoid confusion with rep0
                 )
    if np.ravel(pd.isnull(X)).any():
        raise ValueError("Found null values in grid")
    return X

def _index_to_path(set_folder, set_name,  index):
    return set_folder/('%s_%04d.dat' % (set_name, index))

def generate_replica0(pdf, kin_grids=None, extra_fields=None):
    """ Generates a replica 0 as an average over an existing set of LHAPDF
        replicas and outputs it to the PDF's parent folder

    Parameters
    -----------
    pdf : validphys.core.PDF
        An existing validphys PDF object from which the average replica will be
        (re-)computed

    kin_grids: Grids in (x,Q) used to print replica0 upon. If None, the grids
        of the source replicas are used.
    """

    if extra_fields is not None:
        raise NotImplementedError()

    set_info = pdf.infopath
    set_root = set_info.parent
    if not set_root.exists():
        raise RuntimeError(f"Target directory {set_root} does not exist")

    loaded_grids = {}
    grids = []

    for irep in range(1, len(pdf)):
        if irep in loaded_grids:
            grid = loaded_grids[irep]
        else:
            _header, grid = load_replica(pdf, irep, kin_grids=kin_grids)
            loaded_grids[irep] = grid
        grids.append(grid)
    # This takes care of failing if headers don't match
    try:
        M = rep_matrix(grids)
    except ValueError as e:
        raise ValueError("Null values found in replica grid matrix. "
                         "This may indicate that the headers don't match"
                         "If this is intentional try using use_rep0grid=True") from e
    header = b'PdfType: central\nFormat: lhagrid1\n'
    write_replica(0, set_root, header, M.mean(axis=1))

def new_pdf_from_indexes(
        pdf, indexes, set_name=None, folder=None,
        extra_fields=None, installgrid=False, use_rep0grid=False):
    """Create a new PDF set from by selecting replicas from another one.

    Parameters
    -----------
    pdf : validphys.core.PDF
        An existng validphys PDF object from which the indexes will be
        selected.
    indexes : Iterable[int]
        An iterable with integers corresponding to files in the LHAPDF set.
        Note that replica 0 will be calculated for you as the mean of the
        selected replicas.
    set_name : str
        The name of the new PDF set.
    folder : str, bytes, os.PathLike
        The path where the LHAPDF set will be written. Must exsist.
    installgrid: bool, optional, default=``False``.
        Whether to copy the grid to the LHAPDF path.
    use_rep0grid: bool, optional, default=``False``
        Whether to fill the original replica 0 grid when computing replica 0,
        instead of relying that all grids are the same and averaging the
        files directly. It is slower and will call LHAPDF to fill the grids,
        but works for sets where the replicas have different grids.
    """

    if extra_fields is not None:
        raise NotImplementedError()

    if folder is None:
        folder = pathlib.Path()

    set_root = folder/set_name
    if set_root.exists():
        log.warning("Target directory for new PDF already exists %s. "
                    "Deleting contents.", set_root)
        if set_root.is_dir():
            shutil.rmtree(str(set_root))
        else:
            set_root.unlink()

    set_root.mkdir()

    original_info = pdf.infopath
    original_folder = original_info.parent

    new_info = set_root/(set_name + '.info')

    new_len = len(indexes)+1

    with original_info.open() as orig_file, new_info.open('w') as new_file:
        for line in orig_file:
            if line.find('SetDesc') >= 0:
                new_file.write('SetDesc: "Reweighted set from %s"\n' % pdf)
            elif line.find('NumMembers') >=0:
                new_file.write('NumMembers: %d\n' % new_len)
            else:
                new_file.write(line)

    if use_rep0grid:
        _, rep0grid = load_replica(pdf, 0)
    else:
        rep0grid = None

    for newindex,oldindex in enumerate(indexes, 1):
        original_path = _index_to_path(original_folder, pdf, oldindex)
        new_path = _index_to_path(set_root, set_name, newindex)
        shutil.copy(original_path, new_path)

    # Generate replica 0
    oldpaths = lhapdf.paths()
    try:
        lhapdf.setPaths([str(folder)])
        generatedPDF = PDF(set_name)
        generate_replica0(generatedPDF, rep0grid)
    finally:
        lhapdf.setPaths(oldpaths)

    if installgrid:
        newpath = pathlib.Path(lhaindex.get_lha_datapath()) /  set_name
        log.info(f"Installing new PDF set at {newpath}")
        shutil.copytree(set_root, newpath)


def hessian_from_lincomb(pdf, V, set_name=None, folder = None, extra_fields=None):
    """Construct a new LHAPDF grid from a linear combination of members"""

    # preparing output folder
    neig = V.shape[1]

    base = pathlib.Path(lhapdf.paths()[-1])  / pdf.name
    if set_name is None:
        set_name = pdf.name + "_hessian_" + str(neig)
    if folder is None:
        folder = ''
    set_root = pathlib.Path(folder) / set_name
    # In case a Hessian PDF of the same name already exists, we first remove it. Not doing this
    # can lead to the wrong result if Neig is not the same between both PDF sets.
    if os.path.exists(set_root):
        shutil.rmtree(set_root)
        log.warning("Target directory for new PDF, %s, already exists. " "Removing contents.",
                set_root,)
    os.makedirs(os.path.join(set_root))

    # copy replica 0
    shutil.copy(base/f'{pdf}_0000.dat', set_root / f"{set_name }_0000.dat")

    with open(base/f'{pdf}.info', 'r') as inn, \
         open(set_root  / f'{set_name }.info', 'w') as out:

        for l in inn.readlines():
            if l.find("SetDesc:") >= 0:
                out.write(f"SetDesc: \"Hessian {pdf}_hessian\"\n")
            elif l.find("NumMembers:") >= 0:
                out.write(f"NumMembers: {neig+1}\n")
            elif l.find("ErrorType: replicas") >= 0:
                out.write("ErrorType: symmhessian\n")
            else:
                out.write(l)
        if extra_fields is not None:
            yaml.dump(extra_fields, out, default_flow_style=False)

    _headers, grids = load_all_replicas(pdf)
    result  = (big_matrix(grids).dot(V)).add(grids[0], axis=0, )
    hess_header = b"PdfType: error\nFormat: lhagrid1\n"
    for column in result.columns:
        write_replica(column + 1, set_root, hess_header, result[column])
    log.info("Hessian PDF stored at %s", set_root)
    return set_root
