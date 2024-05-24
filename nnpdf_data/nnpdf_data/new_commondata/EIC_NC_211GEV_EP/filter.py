from pathlib import Path

import numpy as np

from nnpdf_data.filter_utils.eic_utils import (
    fluctuate_data,
    read_central_values,
    read_txt_data,
    write_data,
)

np.random.seed(1234567890)


if __name__ == "__main__":
    input_txt = Path("./rawdata/EIC_18_275_A1c_100fb-1.txt")
    df = read_txt_data(input_txt)
    cv_preds = read_central_values(Path("./rawdata/EIC_NC_211GEV_EP.yaml"))
    fluctuated_cv = fluctuate_data(cv_preds, df["abs"].values)
    # Flat systematic error of 3.2% due to beam polarization
    # from Yuxiang suggestion
    sys_error = 0.032
    write_data(df, asym=False, abserr=fluctuated_cv, add_fluctuate=True, sys_error=sys_error)
