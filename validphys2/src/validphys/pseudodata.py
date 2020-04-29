# -*- coding: utf-8 -*-
"""
pseudodata.py

Tools to obtain and analyse the pseudodata that was seen by the neural
networks during the fitting of a fit
"""
import logging
import multiprocessing as mp

import numpy as np

log = logging.getLogger(__name__)


def fitted_pseudodata(fit, experiments, num_fitted_replicas, t0pdfset=None, NPROC=None):
    """A function to obtain information about the pseudodata that went
        into an N3FIT fit. Note:
            - this function returns the pseudodata for the replicas
              pre-postfit. Postfit discards some replicas and rearranges
              the order. The correpsondence is done by the get_pseudodata
              function.
            - this code runs in parallel to increase efficiency.

        Parameters
        ----------
        fit: validphys.core.FitSpec
        experiments:
            List of validphys.core.ExeperimentSpec
        num_nnfit_replicas: int
            Provided for by validphys.fitdata. Equal to the number of
            pre-postfit replicas.
        t0pdfset: validphys.core.PDF
        NPROC: int
            Integer specifying how many cores to run on. Default is
            mp.cpu_count()

        Example
        -------
        Create a .yaml file say runcard_for_pseudodata.yaml
            pdf: PN3_DIS_130519
            fit: PN3_DIS_130519

            experiments:
              from_: fit

            theory:
              from_: fit

            t0pdfset:
              from_: datacuts

            datacuts:
              from_: fit

            theoryid:
              from_: theory

            use_cuts: fromfit

        Then run
            >>> with open("./runcard_for_pseudodata.yaml", 'r') as stream:
            ...     from reportengine.compat import yaml
            ...     runcard = yaml.safe_load(stream)
            >>> from validphys.api import API
            >>> API.get_pseudodata(**runcard)
    """
    import n3fit.io.reader as reader
    from n3fit.performfit import initialize_seeds

    if t0pdfset is not None:
        t0pdfset = t0pdfset.load_t0()
    else:
        t0pdfset = None

    # The + 1 coming from the fact that we wish to
    # include the last replica
    replica = range(1, num_fitted_replicas + 1)

    trvlseed, nnseed, mcseed, genrep = [
        fit.as_input().get("fitting").get(i)
        for i in ["trvlseed", "nnseed", "mcseed", "genrep"]
    ]

    seeds = initialize_seeds(replica, trvlseed, nnseed, mcseed, genrep)

    def task(d, mcseeds, trvlseeds, replicas):
        all_exp_infos = [[] for _ in range(len(mcseeds))]
        for exp in experiments:
            all_exp_dicts = reader.common_data_reader(
                exp, t0pdfset, replica_seeds=mcseeds, trval_seeds=trvlseeds
            )
            for i, exp_dict in enumerate(all_exp_dicts):
                all_exp_infos[i].append(exp_dict)
        for i, j in zip(all_exp_infos, replicas):
            d[j] = i

    with mp.Manager() as manager:
        d = manager.dict()

        if NPROC is None:
            NPROC = mp.cpu_count()
            log.warning(
                f"Using all {NPROC} cores available, this may be dangerous "
                "especially for use on a cluster. Consider setting the NPROC "
                "variable to something sensible."
            )
        processes = []
        # XXX: There must be a better way to do this. Note it changes
        # from type int to numpy int and thus require being converted back
        batched_mcseeds = np.array_split(seeds.mcseeds, NPROC)
        batched_trvlseeds = np.array_split(seeds.trvlseeds, NPROC)
        batched_replica_num = np.array_split(replica, NPROC)
        for mc_batch, trvl_batch, replica_batch in zip(
            batched_mcseeds, batched_trvlseeds, batched_replica_num
        ):
            p = mp.Process(
                target=task,
                args=(
                    d,
                    list([int(i) for i in mc_batch]),
                    list([int(i) for i in trvl_batch]),
                    list([int(i) for i in replica_batch]),
                ),
            )
            p.start()
            processes.append(p)
        for p in processes:
            p.join()
        pseudodata_dicts = dict(d)
    return pseudodata_dicts


def get_pseudodata(fitted_pseudodata, fitted_replica_indexes):
    return [fitted_pseudodata[i] for i in fitted_replica_indexes]
