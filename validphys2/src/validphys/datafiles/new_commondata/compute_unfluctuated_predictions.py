import lhapdf
import pineappl
import yaml
import numpy as np

from pathlib import Path
from rich.console import Console

console = Console()


DATASET_NAMES = [
    "ATHENA_NC_105GEV_EP",
    # "ATHENA_NC_140GEV_EP",
    "ATHENA_NC_29GEV_EP",
    "ATHENA_NC_45GEV_EP",
    "ATHENA_NC_63GEV_EP",
    # "EIcC_NC_15GEV_EP",
    # "EIcC_NC_22GEV_EP",
    # "EIC_NC_211GEV_EP",
    # "EIC_NC_43GEV_EP",
    # "EIC_NC_67GEV_EP",
]

PDFSET_NAME = "240319-tr-poldis-nnlo-w2min"


def cfactor_loading(dataname: str, theory_id: int = 822) -> np.ndarray:
    foldername = dataname
    kpath = f"/data/theorie/tanjona/miniconda3/envs/nnpdf/share/NNPDF/theories/theory_{theory_id}/cfactor"
    cpath = f"{kpath}/CF_NRM_{foldername}_G1.dat"
    skiprow = 8 if dataname.startswith("ATHENA") else 9
    cpreds = np.loadtxt(cpath, skiprows=skiprow)
    console.print(f"Reading K-factor {cpath}.", style="bold cyan")
    return cpreds[:, 0]


def compute_predictions(dataname: str, theory_id: int = 822) -> np.ndarray:
    foldername = dataname

    # Prepare & laoad the FK tables
    grids_name = f"{foldername}_G1.pineappl.lz4"
    fkfolder = f"/data/theorie/tanjona/miniconda3/envs/nnpdf/share/NNPDF/theories/theory_{theory_id}/fastkernel"
    fk_paths = Path(f"{fkfolder}/{grids_name}")
    grid = pineappl.grid.Grid.read(fk_paths) 
    console.print(f"Reading pineappl grid {fk_paths}.", style="bold cyan")

    pdf_loaded = lhapdf.mkPDF(PDFSET_NAME)

    predictions = grid.convolute_with_one(
        2212,
        pdf_loaded.xfxQ2,
        pdf_loaded.alphasQ2,
    )
    return predictions


def apply_kfactors(preds: np.ndarray, kfactor: np.ndarray) -> np.ndarray:
    return preds * kfactor


if __name__ == "__main__":
    for dataset in DATASET_NAMES:
        console.print(f"Computing {dataset}", style="bold red")
        cfactor = cfactor_loading(dataset)
        predics = compute_predictions(dataset)

        predics_kfacs = apply_kfactors(predics, cfactor)

        predictions_central_yaml = {"predictions_central": predics_kfacs.tolist()}
        with open(f"{dataset}/{dataset}.yaml", "w") as file:
            yaml.dump(predictions_central_yaml, file, sort_keys=False)
