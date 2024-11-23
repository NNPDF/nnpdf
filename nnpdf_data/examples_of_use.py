"""
This file contains examples of use of ``nnpdf_data`` as a library.
This library is currently in pre-alpha form and should not be considered stable.

The functions and examples in this file will be eventually removed but might become
part of the library as an external user-facing interface.

There is currently no user-facing interface so no stability is expected.
"""

from nnpdf_data import path_commondata
from nnpdf_data.commondataparser import parse_new_metadata


def parse_dataset(dataset, variant=None):
    """Given a dataset name, read the observable metadata as a CommonData object.
    A variant can be given.

    The output is a ``ObservableMetaData`` object, with references to all files
    that form the dataset but none of them is loaded.
    This can then be used to _load_ the dataset using load_commondata.

    Example
    -------
    >>> from nnpdf_data.commondataparser import load_commondata
    >>> cd_meta = parse_dataset("LHCB_Z0_7TEV_DIELECTRON_Y")
    >>> cd = load_commondata(cd_meta)
    >>> print(cd)
    CommonData(setname='LHCB_Z0_7TEV_DIELECTRON_Y', ndata=9, commondataproc='DY_Z_Y', nkin=3, nsys=11, legacy=False, legacy_names=['LHCBZ940PB'], kin_variables=['y', 'm_Z2', 'sqrts'])
    """
    setname, observable = dataset.rsplit("_", 1)
    metadata_file = path_commondata / setname / "metadata.yaml"
    metadata = parse_new_metadata(metadata_file, observable, variant=variant)
    return metadata
