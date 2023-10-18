"""
This module contains functions to write commondata and systypes
tables to files
"""


def write_commondata_data(commondata, buffer):
    """
    write commondata table to buffer, this can be a memory map,
    compressed archive or strings (using for instance StringIO)


    Parameters
    ----------

    commondata : validphys.coredata.CommonData

    buffer : memory map, compressed archive or strings
            example: StringIO object


    Example
    -------
    >>> from validphys.loader import Loader
    >>> from io import StringIO

    >>> l = Loader()
    >>> cd = l.check_commondata("NMC").load_commondata_instance()
    >>> sio = StringIO()
    >>> write_commondata_data(cd,sio)
    >>> print(sio.getvalue())

    """
    header = f"{commondata.setname} {commondata.nsys} {commondata.ndata}\n"
    buffer.write(header)
    commondata.commondata_table.to_csv(buffer, sep="\t", header=None)


def write_commondata_to_file(commondata, path):
    """
    write commondata table to file
    """
    with open(path, "w") as file:
        write_commondata_data(commondata, file)


def write_systype_data(commondata, buffer):
    """
    write systype table to buffer, this can be a memory map,
    compressed archive or strings (using for instance StringIO)


    Parameters
    ----------

    commondata : validphys.coredata.CommonData

    buffer : memory map, compressed archive or strings
            example: StringIO object


    Example
    -------
    >>> from validphys.loader import Loader
    >>> from io import StringIO

    >>> l = Loader()
    >>> cd = l.check_commondata("NMC").load_commondata_instance()
    >>> sio = StringIO()
    >>> write_systype_data(cd,sio)
    >>> print(sio.getvalue())

    """
    header = f"{commondata.nsys}\n"
    buffer.write(header)
    commondata.systype_table.to_csv(buffer, sep="\t", header=None)


def write_systype_to_file(commondata, path):
    """
    write systype table to file
    """
    with open(path, "w") as file:
        write_systype_data(commondata, file)
