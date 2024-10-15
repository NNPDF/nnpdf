r"""Implement data form Hepdata reference, data from Table 6.
R is simultaneously determined from the data.
Results are averaged over different $\sqrt{s}$.
"""

import pathlib

HERE = pathlib.Path(__file__).parent

from nnpdf_data.commondata.NMC_NC_NOTFIXED_P_HEPDATA.filter import (
    read_tables,
    write_files,
)

if __name__ == "__main__":
    df = read_tables(HERE / "rawdata", header_line=12)
    write_files(df, HERE)
