from pathlib import Path

import numpy as np
import yaml

from nnpdf_data.filter_utils.poldata_utils import (
    fluctuate_data,
    read_central_values,
    read_txt_data,
    write_data,
)
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

np.random.seed(1234567890)


if __name__ == "__main__":
    input_txt = Path("./rawdata/EICC_3_20_A1c_100fb-1.txt")
    df = read_txt_data(input_txt)
    cv_preds = read_central_values(Path("./rawdata/EICC_NC_15GEV_EP.yaml"))
    fluctuated_cv = fluctuate_data(cv_preds, df["abs"].values)
    write_data(df, asym=False, abserr=fluctuated_cv, add_fluctuate=True)
