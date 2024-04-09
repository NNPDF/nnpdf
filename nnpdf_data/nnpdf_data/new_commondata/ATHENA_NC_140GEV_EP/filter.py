from pathlib import Path

import numpy as np

from nnpdf_data.filter_utils.athena_utils import read_excel, write_data

np.random.seed(1234567890)


if __name__ == "__main__":
    BEAMS = (18, 275)
    input_xlsx = Path("./rawdata/ATHENA_ALL_EP.xlsx")
    xdf = read_excel(input_xlsx, beams=BEAMS)
    write_data(xdf)
