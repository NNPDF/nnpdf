from pathlib import Path

import numpy as np

from nnpdf_data.filter_utils.eic_utils import fluctuate_data, read_cvs, read_excel, write_data

np.random.seed(1234567890)


if __name__ == "__main__":
    BEAMS = (10, 275)
    input_xlsx = Path("./rawdata/ATHENA_ALL_EP.xlsx")
    xdf = read_excel(input_xlsx, beams=BEAMS)
    cv_preds = read_cvs(Path("./rawdata/ATHENA_NC_105GEV_EP.yaml"))
    fluctuated_cv = fluctuate_data(cv_preds, xdf["delta_ALL"].values)
    write_data(xdf, abserr=fluctuated_cv, add_fluctuate=True)
