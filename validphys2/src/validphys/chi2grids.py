"""
chi2grids.py

Compute and store χ² data from replicas, possibly keeping the correlations
between pseudorreplica fluctuations between different fits. This is applied
here to parameter determinations such as those of αs.

This module is severly thwarted by the poor adecuacy of libnnpdf for this use
case. Several pieces of functionality need to be implemented there.
"""
import logging
from collections import namedtuple

import pandas as pd

from reportengine import collect
from reportengine.table import table
from NNPDF import pseudodata, single_replica, RandomGenerator

from validphys.core import PDF
from validphys.results import ThPredictionsResult, DataResult, chi2_breakdown_by_dataset

PseudoReplicaExpChi2Data = namedtuple('PseudoReplicaChi2Data',
    ['experiment', 'dataset', 'ndata', 'chi2', 'nnfit_index'])


log = logging.getLogger(__name__)

def computed_psedorreplicas_chi2(
        experiments, dataseed, pdf, fitted_replica_indexes,
        t0set:(PDF, type(None))):
    """Return a dataframe with the chi² of each replica wirh its corrsponding
    pseudodata (i.e. the one it was fitted with). The chi² is computed for
    both each experiment and each dataset in the experiment. The index of the
    dataframe is

    ['experiment', 'dataset', 'ndata' , 'nnfit_index']

    where 'experiment' is the name of the experiment, 'dataset' is the name of
    the dataset, or "Total" for the total value, 'ndata' is the corresponding
    number of points and 'nnfit_index' is the index specifying the
    pseudorreplica fluctuation.
    """

    #TODO: Everythning about this function is horrible. We need to rewrite
    #experiments.cc from scratch.

    #TODO: Do this somewhere else
    RandomGenerator.InitRNG(0,0)
    if t0set is not None:
        lt0 = t0set.load_t0()
    pdfname = pdf.name
    datas = []

    #No need to save these in the cache, so we call __wrapped__
    original_experiments = [e.load.__wrapped__(e) for e in experiments]
    sqrtcovmat_table = []
    log.debug("Generating dataset covmats")
    for exp in original_experiments:
        if t0set is not None:
            exp.SetT0(lt0)
        #The covariance matrices are currently very expensive to recompute.
        #Store them after computing T0
        sqrtcovmat_table.append([ds.get_sqrtcovmat() for ds in exp.DataSets()])

    for lhapdf_index, nnfit_index in enumerate(fitted_replica_indexes, 1):

        flutuated_experiments = pseudodata(original_experiments, dataseed, nnfit_index)
        lpdf = single_replica(pdfname, lhapdf_index)
        for expspec, exp, mats in zip(experiments, flutuated_experiments, sqrtcovmat_table):
            #We need to manage the memory
            exp.thisown = True

            th = ThPredictionsResult.from_convolution(pdf, expspec,
                loaded_data=exp, loaded_pdf=lpdf)


            results = DataResult(exp), th
            #The experiment already has T0. No need to set it again.
            #TODO: This is a hack. Get rid of this.
            chi2 = chi2_breakdown_by_dataset(results, exp, t0set=None,
                                             datasets_sqrtcovmat=mats)

            for i, (dsname,(value, ndata)) in enumerate(chi2.items()):
                data = PseudoReplicaExpChi2Data(
                    nnfit_index=nnfit_index,
                    experiment=expspec.name,
                    #We set the i so that the sort order is maintaned here.
                    dataset = (i, dsname),
                    ndata = ndata,
                    chi2=value
                    )
                datas.append(data)

    df =  pd.DataFrame(datas, columns=PseudoReplicaExpChi2Data._fields)
    df.set_index(['experiment', 'dataset', 'ndata' , 'nnfit_index'], inplace=True)
    df.sort_index(inplace=True)
    #Now that we have the order we like, we remove the i
    df.index.set_levels([x[1] for x in df.index.levels[1]], level=1, inplace=True)
    return df

#TODO: Probably fitcontext should set all of the variables required to compute
#this. But better setting
#them explicitly than setting some, se we require the user to do that.
fits_computed_psedorreplicas_chi2 = collect(computed_psedorreplicas_chi2, ('fits',))

dataspecs_computed_pseudorreplicas_chi2 = collect(computed_psedorreplicas_chi2, ('dataspecs',))

@table
def export_fits_computed_psedorreplicas_chi2(fits_computed_psedorreplicas_chi2):
    """Hack to force writting the CSV output"""
    return fits_computed_psedorreplicas_chi2
