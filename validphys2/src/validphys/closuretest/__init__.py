"""
closuretest

module containing all actions specific to closure test
"""

from validphys.closuretest.closure_plots import *
from validphys.closuretest.closure_results import *
from validphys.closuretest.multiclosure import *
from validphys.closuretest.multiclosure_output import *
from validphys.closuretest.multiclosure_pdf import *
from validphys.closuretest.multiclosure_pdf_output import *
from validphys.closuretest.multiclosure_preprocessing import *
from validphys.closuretest.multiclosure_pseudodata import *
from validphys.closuretest.inconsistent_closuretest.multiclosure_inconsistent_output import *
from validphys.closuretest.multiclosure_nsigma_helpers import (
    fits_datasets_chi2_nsigma_deviation,
    fits_data,
    is_weighted,
    n_fits,
)
from validphys.closuretest.multiclosure_nsigma import (
    set_1,
    set_2,
    comp_set_1,
    set_3,
    probability_inconsistent,
)
from validphys.closuretest.multiclosure_nsigma_output import (
    plot_all_sets,
    plot_1_minus_all_sets,
    plot_probability_inconsistent,
    plot_probability_consistent,
)
