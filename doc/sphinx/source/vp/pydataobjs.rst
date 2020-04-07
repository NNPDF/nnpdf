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

Currently only FKTables can be directly without C++ code. This is implemented
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
compute theory predictions in pure Python. In particular the
:py:func:`validphys.convolution.predictions` function returns predictions in
terms of PDF and dataset objects that can be obtained directly from `validphys`
runcards::

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
