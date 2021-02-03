.. _pyobjs:

Python based data objects
=========================

Internal data formats such as PDF sets, CommonData, or :ref:`FKTables
<fktables>` files are currently accessed through the `libnnpdf` C++ code
(interfaced trough the SWIG wrappers). However there is a :ref:`project
<https://github.com/NNPDF/nnpdf/issues?q=label%3Adestroyingc%2B%2B+>` underway
to make these resources available in terms of standard Python containers
(particularly numpy arrays and pandas dataframes). The objectives include
simplifying the codebase, increasing the ease of use and enabling more advanced
computation and storage strategies.

Loading FKTables
----------------

This is implemented
in the :py:mod:`validphys.fkparser` module. For example::

    from validphys.fkparser import load_fktable
    from validphys.loader import Loader
    l = Loader()
    fk = l.check_fktable(setname="ATLASTTBARTOT", theoryID=53, cfac=('QCD',))
    res = load_fktable(fk)

results in an :py:mod:`validphys.coredata.FKTableData` object containing all
the information needed to compute a convolution. In particular the ``sigma``
property contains a dataframe representing the partonic cross-section
(including the cfactors).

Computing theory predictions
----------------------------

The :py:mod:`validphys.convolution` module implements the necessary tooling to
compute theory predictions using Numpy and Pandas. In particular the
:py:func:`validphys.convolution.predictions` function returns predictions in
terms of PDF and dataset objects that can be obtained directly from
``validphys`` runcards::

    from validphys.api import API
    from validphys.convolution import predictions

    inp = {
        'dataset_input': {'dataset': 'ATLASTTBARTOT', 'cfac': ['QCD']},
        'theoryid': 162,
        'use_cuts': 'internal',
        'pdf': 'NNPDF31_nnlo_as_0118'
    }

    preds = predictions(API.dataset(**inp), API.pdf(**inp))

    print(preds.values.mean(axis=1))


The usage of standard scientific Python types opens interesting avenues for
parallelization. For example here is how to compute the mean prediction for all
datasets using the `Dask <https://dask.org/>`_ library::

    from dask.distributed import Client

    from validphys.api import API
    from validphys.convolution import predictions

    c = Client()

    inp = {
        'fit': '181023-001-sc',
        'use_cuts': 'internal',
        'theoryid': 162,
        'pdf': 'NNPDF31_nlo_as_0118',
        'experiments': {'from_': 'fit'}
    }


    all_datasets = [ds for e in API.experiments(**inp) for ds in e.datasets]

    pdf = API.pdf(**inp)

    future_pred = dask.delayed(pure=True)(predictions)
    c.gather(c.compute([np.mean(future_pred(ds, pdf), axis=0) for ds in all_datasets]))

Central predictions
^^^^^^^^^^^^^^^^^^^

The default :py:func:`validphys.convolution.predictions` computes one
prediction for each replica in the PDF set (for Monte Carlo PDF sets). The user
is then supposed to average the replica predictions to get a central value. A
quick approximation is to use the central value directly. This is exact for DIS
observables and a generally very good approximation for hadronic observables.
The :py:func:`validphys.convolution.central_predictions` function may be
appropriate for computations where the PDF error is not required, such as the
central χ².

The previous example can be simpler using ``central_predictions``::


    from validphys.api import API
    from validphys.convolution import central_predictions

    inp = {
        'dataset_input': {'dataset': 'ATLASTTBARTOT', 'cfac': ['QCD']},
        'theoryid': 162,
        'use_cuts': 'internal',
        'pdf': 'NNPDF31_nnlo_as_0118'
    }


    central_preds = central_predictions(API.dataset(**inp), API.pdf(**inp))

    print(central_preds)

Linear predictions
^^^^^^^^^^^^^^^^^^

DIS predictions are linear in the difference between PDF and central value, and
hence in the Hessian error parameters. For hadronic observables this is only
true to a good approximation. The
:py:func:`validphys.convolution.linear_predictions` computes approximate
predictions that are linear in the error parameters, and which may be useful in
specific situations. In particular, for such predictions the prediction of the
central replica is the same as the mean of the replica predictions::

    import numpy as np
    from validphys.loader import Loader
    from validphys.convolution import predictions, linear_predictions, central_predictions

    l = Loader()
    pdf = l.check_pdf('NNPDF31_nnlo_as_0118')
    ds = l.check_dataset('ATLASTTBARTOT', theoryid=53, cfac=('QCD',))

    # "Exact" predictions
    p = predictions(ds, pdf).T
    # Approximate predictions, neglecting the quadratic terms in the
    # differences between each replica and the central value.
    lp = linear_predictions(ds, pdf).T
    # Central predictions
    cp = central_predictions(ds, pdf).T


    assert np.allclose(lp.mean(), cp)
    assert not np.allclose(p.mean(), cp)
    # Compute the size of the differences between approximate and true predictions
    # over the PDF uncertainty. Take the maximum over the three ttbar data points.
    print(((p - lp).std() / p.std()).max())

Loading CommonData
------------------

The underlying functions for loading CommonData can be found in
:py:mod:`validphys.results_providers.commondata_parser`. The data is loaded
as :py:class:`validphys.coredata.CommonData`, which uses the
`dataclasses <https://docs.python.org/3/library/dataclasses.html>`_ module
which automatically generates some special methods for the class. The
underlying data is stored as DataFrames, and so can be used
with the standard pandas machinery::

    import pandas as pd

    from validphys.api import API
    from validphys.results_providers.commondata_parser import load_commondata
    # first get the CommonDataSpec
    cd = API.commondata(dataset_input={"dataset":"NMC"})
    lcd = load_commondata(cd)
    assert isinstance(lcd.central_values, pd.Series)
    assert isinstance(lcd.systematics_table, pd.DataFrame)

The :py:class:`validphys.coredata.CommonData` class has a method which returns
a new instance of the class with cuts applied::

    from validphys.api import API
    from validphys.results_providers.commondata_parser import load_commondata
    inp = {
        "dataset_input": {"dataset":"NMC"},
        "use_cuts": "internal",
        "theoryid": 162
    }
    # first get the CommonDataSpec
    cd = API.commondata(**inp)
    lcd = load_commondata(cd)
    # CommonDataSpec object ndata is always total data points uncut
    assert lcd.ndata == cd.ndata
    cuts = API.cuts(**inp)
    lcd_cut = lcd.with_cuts(cuts)
    # data has been cut, ndata should have changed.
    assert lcd_cut.ndata != cd.ndata

An action already exists which returns the loaded and cut commondata, which is
more convenient than calling the underlying functions::

    api_lcd_cut = API.loaded_commondata_with_cuts(**inp)
    assert api_lcd_cut.ndata == lcd_cut.ndata

Loading Covariance Matrices
---------------------------

Functions which take :py:class:`validphys.coredata.CommonData` s and return
covariance matrices can be found in
:py:mod:`validphys.results_providers.covmat_construction`. As with the commondata
the underlying functions can be accessed directly::

    import numpy as np
    from validphys.api import API
    from validphys.results_providers.covmat_construction import covmat_from_systematics

    inp = {
        "dataset_input": {"dataset":"NMC"},
        "use_cuts": "internal",
        "theoryid": 162
    }
    lcd = API.loaded_commondata_with_cuts(**inp)
    cov = covmat_from_systematics(lcd)
    assert isinstance(cov, np.ndarray)
    assert cov.shape == (lcd.ndata, lcd.ndata)

There exists a similar function which acts upon a list of multiple commondatas
and takes into account correlations between datasets::

    from validphys.results_providers.covmat_construction import datasets_covmat_from_systematics
    inp = {
        "dataset_inputs": [
            {"dataset":"NMC"},
            {"dataset":"NMCPD"},
        ],
        "use_cuts": "internal",
        "theoryid": 162
    }
    lcds = API.dataset_inputs_loaded_cd_with_cuts(**inp)
    total_ndata = np.sum([lcd.ndata for lcd in lcds])
    total_cov = datasets_covmat_from_systematics(lcds)
    assert total_cov.shape == (total_ndata, total_ndata)

These functions are already leveraged by actions, which can be accessed directly
from the API::

    from validphys.api import API

    inp = {
        "dataset_input": {"dataset":"NMC"},
        "use_cuts": "internal",
        "theoryid": 162
    }
    # single dataset covmat
    cov = API.experimental_covmat(**inp)
    inp = {
        "dataset_inputs": [
            {"dataset":"NMC"},
            {"dataset":"NMCPD"},
        ],
        "use_cuts": "internal",
        "theoryid": 162
    }
    total_cov = API.dataset_inputs_experimental_covmat(**inp)
