

from ekobox import evol_pdf, gen_theory, gen_op, genpdf
import pathlib 
import numpy as np
import yaml
from . import utils
from validphys.loader import Loader


def load_fit(conf_folder, num_replica):
    """
    Loads a replica pdf at fitting scale in conf_folder and return a lhapdf-like object

    Parameters
    ----------

        conf_folder: str
            path to the folder containing the fit 
        num_replicas: int
            number of replica pdf 
    
    Returns
    -------

            : lhapdf-like
            lhapdf-like object
    """
    nnfitpath = pathlib.Path(conf_folder + "/nnfit")
    rep_path = pathlib.Path("/replica_" + str(num_replica) + "/" + conf_folder + ".exportgrid")
    yaml_file = nnfitpath/rep_path.relative_to(rep_path.anchor)
    with yaml_file.open() as fp:
        data = yaml.safe_load(fp)
    pdf_grid = np.array(data["pdfgrid"]).transpose()
    x_grid = np.array(data["xgrid"]).astype(np.float)
    my_pdf = utils.lhapdf_like(pdf_grid,data["q20"],x_grid)
    return my_pdf

def evolve_fit(conf_folder, max_replicas):
    """
    Evolves num_replicas replica pdfs at fitting scale in conf_folder and dump lhapdf pdfs

    Parameters
    ----------
        conf_folder: str
            path to the folder containing the fit 
        num_replicas: int
            number of replica pdfs 
    """
    initial_PDFs_list = []
    for rep in range(2,max_replicas+1):
        initial_PDFs_list.append(load_fit(conf_folder, rep))
    #read the runcard
    my_runcard = utils.read_runcard(conf_folder)
    #theory_card construction
    theory = Loader().check_theoryID(my_runcard["theory"]["theoryid"]).get_description()
    del theory["FNS"]
    theory_card = gen_theory.gen_theory_card(theory["PTO"],theory["Q0"], update = theory )
    #for testing
    #theory_card = gen_theory.gen_theory_card(0,theory["Q0"] ) 
    #construct operator card
    q2_grid = utils.generate_q2grid(theory["Q0"])
    op_card = gen_op.gen_op_card(q2_grid) #Change also x-grid according to pdf one?
    #generate eko operator (temporary because it will be loaded from theory)
    evol_pdf.gen_out(theory_card,op_card,path = conf_folder)
    #evolve and dump
    evol_pdf.evolve_pdfs(initial_PDFs_list, theory_card, op_card,path = conf_folder, name = "Evolved_fit", info_update={"SetDesc": my_runcard["description"], "OrderQCD":theory["PTO"] } ) #name to be changed and also info
