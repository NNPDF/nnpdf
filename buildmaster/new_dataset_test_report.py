import pandas as pd
import numpy as np
import pathlib

import seaborn as sns

from reportengine.table import table
from reportengine.figure import figure

from validphys.commondataparser import parse_commondata_new, parse_new_metadata
from validphys import covmats
from validphys.calcutils import calc_chi2
from validphys import plotutils



def parsed_metadata(dataset_name, observable_name, variant):
    metadata_file=pathlib.Path(f"/Users/markcostantini/codes/nnpdfgit/nnpdf/buildmaster/{dataset_name}/metadata.yaml")
    if variant == '':
        return parse_new_metadata(metadata_file=metadata_file, observable_name=observable_name, variants=[])
    return parse_new_metadata(metadata_file=metadata_file, observable_name=observable_name, variants=[variant])

def parsed_commondata(parsed_metadata):
    return parse_commondata_new(parsed_metadata)


@table
def table_comparison_old_new_dataset(dataset, pdf, parsed_commondata, variant=''):
    """
    Generates a table to be used to compare a legacy dataset 
    against a newly implemented dataset.


    Parameters
    ----------
    dataset: validphys.core.DataSetSpec

    pdf: validphys.core.PDF

    new_dataset_name: str
                name of new dataset, should end with name observable
                e.g. _PTY

    variant: str
            name of the variant to use, e.g. 'decorrelated'

    Returns
    -------
    pd.DataFrame

    """
    cd_new = parsed_commondata

    cd_old = dataset.load_commondata()

    covmat_new = covmats.covmat_from_systematics(cd_new, dataset, use_weights_in_covmat=False)
    sqrt_covmat_new = covmats.sqrt_covmat(covmat_new)

    covmat_old = covmats.covmat_from_systematics(cd_old, dataset, use_weights_in_covmat=False)
    sqrt_covmat_old = covmats.sqrt_covmat(covmat_old)

    diff = cd_old.central_values - covmats.dataset_t0_predictions(dataset, pdf)

    try:
        chi2_new = calc_chi2(sqrt_covmat_new, diff)
    except:
        import IPython; IPython.embed()
    chi2_old = calc_chi2(sqrt_covmat_old, diff)

    covmats_coincide = np.allclose(covmat_new, covmat_old)
    cv_coincide = np.allclose(cd_new.central_values, cd_old.central_values)

    df = pd.DataFrame(
        {'values':[
            chi2_new/cd_new.ndata,
            chi2_old/cd_old.ndata,
            covmats_coincide,
            cd_new.ndata,
            cd_old.ndata,
            cd_new.nsys,
            cd_old.nsys,
        ]},
        index=[
            'chi2_new',
            'chi2_old',
            'new covmat == old covmat',
            'ndata new',
            'ndata old',
            'nsys new',
            'nsys old',
        ],
    )

    return df


@figure
def figure_comparison_correlation_matrices(dataset, parsed_commondata, variant=''):
    """
    Heat map plot of correlation matrix of legacy dataset vs new dataset.
    """

    cd_new = parsed_commondata

    cd_old = dataset.load_commondata()

    covmat_new = covmats.covmat_from_systematics(cd_new, dataset, use_weights_in_covmat=False)

    covmat_old = covmats.covmat_from_systematics(cd_old, dataset, use_weights_in_covmat=False)


    fig, axs = plotutils.subplots(ncols=2, nrows=1, figsize=(12, 5), sharey=True)

    # Create a shared colorbar axis
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])

    sns.heatmap(
        covmat_new / np.outer(np.sqrt(np.diag(covmat_new)), np.sqrt(np.diag(covmat_new))),
        annot=False,
        cmap="YlGnBu",
        ax=axs[0],
        cbar_ax=cbar_ax,
    )

    sns.heatmap(
        covmat_old / np.outer(np.sqrt(np.diag(covmat_old)), np.sqrt(np.diag(covmat_old))),
        annot=False,
        cmap="YlGnBu",
        ax=axs[1],
        cbar_ax=cbar_ax,
    )

    axs[0].set_title("Legacy Correlation Matrix")
    axs[1].set_title("New Correlation Matrix")
    
    return fig