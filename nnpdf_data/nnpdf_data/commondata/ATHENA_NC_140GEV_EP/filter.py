from pathlib import Path

import numpy as np
import yaml

from nnpdf_data.filter_utils.poldata_utils import read_excel, write_data
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

np.random.seed(1234567890)


if __name__ == "__main__":
    BEAMS = (18, 275)
    input_xlsx = Path("./rawdata/ATHENA_ALL_EP.xlsx")
    xdf = read_excel(input_xlsx, beams=BEAMS)
    write_data(xdf)
