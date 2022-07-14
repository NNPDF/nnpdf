from ekobox import evol_pdf, gen_theory, gen_op, genpdf, gen_info
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

    Parameters
    ----------

        conf_folder: str or pathlib.Path
            path to the folder containing the fit
    """
    usr_path = pathlib.Path(conf_folder)
    initial_PDFs_dict = load_fit(usr_path)
    eko, theory, op = construct_eko_for_fit(usr_path)
    info = gen_info.create_info_file(theory, op, 1, info_update={})  # to be changed
    dump_info_file(usr_path, info)
    for replica in initial_PDFs_dict.keys():
        evolved_block = evolve_exportgrid(initial_PDFs_dict[replica], eko, theory, op)
        dump_evolved_replica(
            evolved_block, usr_path, int(replica.removeprefix("replica_"))
        )


def load_fit(usr_path):
    """
    Loads all the replica pdfs at fitting scale in usr_path and return the exportgrids

    Parameters
    ----------

        usr_path: pathlib.Path
            path to the folder containing the fit  
    
    Returns
    -------

            : dict
            exportgrids info
    """
    nnfitpath = usr_path / "nnfit"
    replica_list = []
    for subdir, dirs, files in os.walk(nnfitpath):
        for dir in dirs:
            replica_list.append(dir)
    replica_list.remove("input")
    # remove the eventual evolution folder
    try:
        replica_list.remove(usr_path.stem)
    except:
        pass
    pdf_dict = {}
    for replica in replica_list:
        rep_path = pathlib.Path(replica) / (usr_path.stem + ".exportgrid")
        yaml_file = nnfitpath / rep_path.relative_to(rep_path.anchor)
        with yaml_file.open() as fp:
            data = yaml.safe_load(fp)
        pdf_dict[replica] = data
    return pdf_dict


# Temporary solution. Then it will be loaded from the theory itself
def construct_eko_for_fit(usr_path):
    """
    Construct the eko operator needed for evolution of fitted pdfs

    Parameters
    ----------
        usr_path: pathlib.Path
            path to the folder containing the fit  
    Returns
    -------
        : eko.output.Output
        eko operator
        : dict
        theory runcard
        : dict
        operator runcard
    """
    # read the runcard
    my_runcard = utils.read_runcard(usr_path)
    # theory_card construction
    theory = Loader().check_theoryID(my_runcard["theory"]["theoryid"]).get_description()
    theory.pop("FNS")
    t_card = gen_theory.gen_theory_card(theory["PTO"], theory["Q0"], update=theory)
    # construct operator card
    q2_grid = utils.generate_q2grid(theory["Q0"], 1.0e5)
    op_card = gen_op.gen_op_card(q2_grid)
    # generate eko operator (temporary because it will be loaded from theory)
    eko_op = run_dglap(t_card, op_card)
    return eko_op, t_card, op_card


def evolve_exportgrid(exportgrid, eko, theory_card, operator_card):
    """
    Evolves the provided exportgrid for the desired replica with the eko and returns the evolved block

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
    Returns
    -------
        : np.array
        evolved block
    """
    # construct LhapdfLike object
    pdf_grid = np.array(exportgrid["pdfgrid"]).transpose()
    x_grid = np.array(exportgrid["xgrid"]).astype(np.float)
    pdf_to_evolve = utils.LhapdfLike(pdf_grid, exportgrid["q20"], x_grid)
    # evolve pdf
    evolved_pdf = apply.apply_pdf(eko, pdf_to_evolve)
    # generate block to dump
    targetgrid = operator_card["interpolation_xgrid"]
    block = genpdf.generate_block(
        lambda pid, x, Q2, evolved_PDF=evolved_pdf: targetgrid[targetgrid.index(x)]
        * evolved_PDF[Q2]["pdfs"][pid][targetgrid.index(x)],
        xgrid=targetgrid,
        Q2grid=operator_card["Q2grid"],
        pids=br.flavor_basis_pids,
    )
    return block


def dump_evolved_replica(evolved_block, usr_path, replica_num):
    path_where_dump = usr_path / "nnfit" / usr_path.stem
    # create folder to dump the evolved replica if it does not exist
    if not os.path.exists(path_where_dump):
        os.makedirs(path_where_dump)
    to_write_in_head = "PdfType: replica\nFromMCReplica: " + str(replica_num) + "\n"
    genpdf.export.dump_blocks(
        path_where_dump, replica_num, [evolved_block], pdf_type=to_write_in_head
    )


def dump_info_file(usr_path, info):
    # dump info file if not already there
    path_where_dump = usr_path / "nnfit" / usr_path.stem
    info_path = path_where_dump / (usr_path.stem + ".info")
    if not info_path.is_file():
        genpdf.export.dump_info(path_where_dump, info)
