r"""Implement data form Hepdata reference, data from Table 5.
R is simultaneously determined from the data.
Results are averaged onver different $\sqrt{s}$.
"""

import pathlib

from nnpdf_data.filter_utils.nmc_hepdata_utils import read_tables, write_files

HERE = pathlib.Path(__file__).parent

if __name__ == "__main__":
    df = read_tables(HERE / "rawdata", header_line=14)
    write_files(df, HERE)
