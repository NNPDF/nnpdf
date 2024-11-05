from pathlib import Path

import numpy as np
import yaml

from nnpdf_data.filter_utils.poldata_utils import (
    fluctuate_data,
    read_central_values,
    read_excel,
    write_data,
)
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

np.random.seed(1234567890)


if __name__ == "__main__":
    BEAMS = (5, 100)
    input_xlsx = Path("./rawdata/ATHENA_ALL_EP.xlsx")
    xdf = read_excel(input_xlsx, beams=BEAMS)
    cv_preds = read_central_values(Path("./rawdata/ATHENA_NC_45GEV_EP.yaml"))
    fluctuated_cv = fluctuate_data(cv_preds, xdf["delta_ALL"].values)
    write_data(xdf, abserr=fluctuated_cv, add_fluctuate=True)
