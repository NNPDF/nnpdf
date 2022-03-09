from scipy.interpolate import interp1d
import numpy as np
import math
import pathlib
import yaml

class lhapdf_like():
    pids_dict = {-6:'TBAR', -5:'BBAR', -4:'CBAR', -3:'SBAR', -2:'UBAR', -1:'DBAR', 21:'GLUON', 1:'D', 2:'U', 3:'S', 4:'C', 5:'B', 6:'T', 22:'PHT' }
    pids_order = ['TBAR', 'BBAR', 'CBAR', 'SBAR', 'UBAR', 'DBAR', 'GLUON', 'D', 'U', 'S', 'C', 'B', 'T', 'PHT' ]
    def __init__(self, pdf_grid, q20, x_grid):
        self.pdf_grid = pdf_grid
        self.q20 = q20 
        self.x_grid = x_grid
        self.funcs = [interp1d(self.x_grid, self.pdf_grid[pid]) for pid in range(len(self.pids_order))]
    def xfxQ2(self, pid, x, q2 ):
        if not math.isclose(q2,self.q20,rel_tol = 1e-6):
            raise ValueError("The q2 requested is not the fitting scale of this pdf")
        return x*(self.funcs[self.pids_order.index(self.pids_dict[pid])](x))
    def hasFlavor(self,pid):
        if pid in self.pids_dict.keys():
            return True 
        return False


def read_runcard(conf_folder):
    """
    reads the runcard and returns the relevant information for evolven3fit
    """
    runcard_path = pathlib.Path(conf_folder + "/filter.yml")
    with runcard_path.open() as fp:
        data = yaml.safe_load(fp)
    return data

def generate_q2grid(Q0):
    """
    Generate the q2grid used in the final evolved pdfs (Temporary solution)
    """
    Qfin = 1.e5
    grid = np.geomspace(Q0**2,Qfin**2, num = 100)
    #for testing
    #grid = np.array([100.]) 
    return grid

    
