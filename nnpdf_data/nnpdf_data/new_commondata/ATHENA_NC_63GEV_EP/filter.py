from pathlib import Path

import numpy as np

from nnpdf_data.filter_utils.eic_utils import (
    fluctuate_data,
    read_central_values,
    read_excel,
    write_data,
)

np.random.seed(1234567890)


if __name__ == "__main__":
    BEAMS = (10, 100)
    input_xlsx = Path("./rawdata/ATHENA_ALL_EP.xlsx")
    xdf = read_excel(input_xlsx, beams=BEAMS)

    for pto in ["NLO", "NNLO"]:
        cv_preds = read_central_values(Path(f"./rawdata/ATHENA_NC_63GEV_EP_{pto}.yaml"))
        fluctuated_cv = fluctuate_data(cv_preds, xdf["delta_ALL"].values)
        write_data(xdf, abserr=fluctuated_cv, add_fluctuate=True, suffix=pto.lower())
