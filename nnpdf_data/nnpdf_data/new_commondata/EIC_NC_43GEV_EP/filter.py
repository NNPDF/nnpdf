from pathlib import Path

import numpy as np

from nnpdf_data.filter_utils.eic_utils import fluctuate_data, read_cvs, read_txt_data, write_data

np.random.seed(1234567890)


if __name__ == "__main__":
    input_txt = Path("./rawdata/EIC_5_41_A1c_100fb-1.txt")
    df = read_txt_data(input_txt)
    cv_preds = read_cvs(Path("./rawdata/EIC_NC_43GEV_EP.yaml"))
    fluctuated_cv = fluctuate_data(cv_preds, df["abs"].values)
    write_data(df, asym=False, abserr=fluctuated_cv, add_fluctuate=True)
