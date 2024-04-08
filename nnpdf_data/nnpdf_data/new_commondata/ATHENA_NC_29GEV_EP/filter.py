from pathlib import Path

import numpy as np
import yaml

from nnpdf_data.athena_utils import read_excel, fluctuate_data, write_data

np.random.seed(1234567890)


def read_cvs() -> np.ndarray:
    cv_preds = Path("./ATHENA_NC_29GEV_EP.yaml")
    cv_yaml = yaml.safe_load(cv_preds.read_text())
    return np.array(cv_yaml["predictions_central"])


if __name__ == "__main__":
    BEAMS = (5, 41)
    input_xlsx = Path("./rawdata/ATHENA_ALL_EP.xlsx")
    xdf = read_excel(input_xlsx, beams=BEAMS)
    cv_preds = read_cvs()
    fluctuated_cv = fluctuate_data(cv_preds, xdf["delta_ALL"].values)
    write_data(xdf, abserr=fluctuated_cv, add_fluctuate=True)
