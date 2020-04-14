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


def get_pseudodata(fit, pdf, experiments, NPROC=None):
    """A function to obtain information about the pseudodata that went
        into an N3FIT fit. Note this code runs in parallel to increase efficiency.

        Parameters
        ----------
        fit: validphys.core.FitSpec
        experiments:
            List of validphys.core.ExeperimentSpec
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

    t0pdfset = pdf.load_t0()
    # Note: len(pdf) = num replicas + 1
    # due to rep 0, which we do not consider
    # here.
    replica = range(1, len(pdf))

    trvlseed, nnseed, mcseed, genrep = [
        fit.as_input().get("fitting").get(i)
        for i in ["trvlseed", "nnseed", "mcseed", "genrep"]
    ]

    seeds = initialize_seeds(replica, trvlseed, nnseed, mcseed, genrep)

    def task(L, mcseeds, trvlseeds):
        all_exp_infos = [[] for _ in range(len(mcseeds))]
        for exp in experiments:
            all_exp_dicts = reader.common_data_reader(
                exp, t0pdfset, replica_seeds=mcseeds, trval_seeds=trvlseeds
            )
            for i, exp_dict in enumerate(all_exp_dicts):
                all_exp_infos[i].append(exp_dict)
        L.extend(all_exp_infos)

    with mp.Manager() as manager:
        L = manager.list()

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
        for mc_batch, trvl_batch in zip(batched_mcseeds, batched_trvlseeds):
            p = mp.Process(
                target=task,
                args=(
                    L,
                    list([int(i) for i in mc_batch]),
                    list([int(i) for i in trvl_batch]),
                ),
            )
            p.start()
            processes.append(p)
        for p in processes:
            p.join()
        pseudodata_dicts = [i for i in L]
    return pseudodata_dicts
