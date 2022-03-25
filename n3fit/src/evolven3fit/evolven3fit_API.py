from ekobox import evol_pdf, theory_card, operators_card, genpdf, info_file
import pathlib
import numpy as np
import yaml
from . import utils
from validphys.loader import Loader
import os
from ekomark import apply
from eko import basis_rotation as br
from eko import run_dglap

def evolve_fit(conf_folder):
    """
    Evolves all the fitted replica in conf_folder/nnfit 
    """
    initial_PDFs_dict = load_fit(conf_folder)
    eko, theory, op = construct_eko_for_fit(conf_folder)
    for replica in initial_PDFs_dict.keys():
        evolve_exportgrid(
            initial_PDFs_dict[replica],
            eko,
            theory,
            op,
            int(replica.removeprefix("replica_")),
            conf_folder,
        )


def load_fit(conf_folder):
    """
    Loads all the replica pdfs at fitting scale in conf_folder and return the exportgrids

    Parameters
    ----------

        conf_folder: str
            path to the folder containing the fit  
    
    Returns
    -------

            : dict
            exportgrids info
    """
    nnfitpath = pathlib.Path(conf_folder) / "nnfit"
    replica_list = []
    for subdir, dirs, files in os.walk(nnfitpath):
        for dir in dirs:
            replica_list.append(dir)
    replica_list.remove("input")
    pdf_dict = {}
    for replica in replica_list:
        rep_path = pathlib.Path(replica) / (conf_folder + ".exportgrid")
        yaml_file = nnfitpath / rep_path.relative_to(rep_path.anchor)
        with yaml_file.open() as fp:
            data = yaml.safe_load(fp)
        pdf_dict[replica] = data
    return pdf_dict


# Temporary solution. Then it will be loaded from the theory itself
def construct_eko_for_fit(conf_folder):
    """
    Construct the eko operator needed for evolution of fitted pdfs

    Parameters
    ----------
        conf_folder: str
            path to the folder containing the fit  
    """
    # read the runcard
    my_runcard = utils.read_runcard(conf_folder)
    # theory_card construction
    theory = Loader().check_theoryID(my_runcard["theory"]["theoryid"]).get_description()
    theory.pop("FNS")
    t_card = theory_card.generate(theory["PTO"], theory["Q0"], update=theory)

    # construct operator card
    q2_grid = utils.generate_q2grid(theory["Q0"], 1.0e5)
    op_card = operators_card.generate(q2_grid)
    # generate eko operator (temporary because it will be loaded from theory)
    eko = run_dglap(t_card, op_card)
    return eko, t_card, op_card
    #


def evolve_exportgrid(
    exportgrid, eko, theory_card, operator_card, replica, conf_folder
):
    """
    Evolves the provided exportgrid for the desired replica with the eko and dump the replica in conf_folder/nnfit

    Parameters
    ----------
        exportgrid: dict
            exportgrid of pdf at fitting scale
        eko: eko object
            eko operator for evolution
        theory_card: dict
            theory card
        operator_card: dict
            operator card
        replica: int
            num replica
        conf_folder: str
            path to fit folder
    """
    path_where_dump = pathlib.Path(conf_folder) / "nnfit" / conf_folder
    # construct LhapdfLike object
    pdf_grid = np.array(exportgrid["pdfgrid"]).transpose()
    x_grid = np.array(exportgrid["xgrid"]).astype(np.float)
    pdf_to_evolve = utils.LhapdfLike(pdf_grid, exportgrid["q20"], x_grid)
    # evolve pdf
    evolved_pdf = apply.apply_pdf(eko, pdf_to_evolve)
    # construct and dump info file if not already there
    info_path = pathlib.Path(conf_folder) / "nnfit" / (conf_folder + ".info")
    if not info_path.is_file():
        info = info_file.build(
            theory_card, operator_card, 1, info_update={}
        )  # to be changed
        genpdf.export.dump_info(path_where_dump, info)
    # generate block to dump
    targetgrid = operator_card["xgrid"]
    block = genpdf.generate_block(
        lambda pid, x, Q2, evolved_PDF=evolved_pdf: targetgrid[targetgrid.index(x)]
        * evolved_PDF[Q2]["pdfs"][pid][targetgrid.index(x)],
        xgrid=targetgrid,
        Q2grid=operator_card["Q2grid"],
        pids=br.flavor_basis_pids,
    )
    to_write_in_head = "PdfType: replica\nFromMCReplica: " + str(replica) + "\n"
    genpdf.export.dump_blocks(
        path_where_dump, replica, [block], pdf_type=to_write_in_head
    )
