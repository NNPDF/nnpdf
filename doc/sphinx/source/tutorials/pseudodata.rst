.. _pseudodata:
Obtaining the pseudodata used by an ``n3fit`` fit
=================================================

Suppose one has obtained a fit using the ``n3fit`` framework and wants to do some analysis that requires
knowing exactly the data that the neural networks saw during the fitting procedure. Thankfully, this
information is reproducible due to the various seeds in the fit runcard.

The 3 seeds of interest are ``trvlseed`` which determines the training/validation splitting, ``nnseed``
which concerns the initialization of the neural netowrks themselves, and finally ``mcseed`` which is the
seed used by the pseudodata generation. Clearly, the ones we are interested in are ``trvlseed`` and ``mcseed``.

This functionality is exposed through the API by using
:py:func:`validphys.pseudodata.recreate_fit_pseudodata` which will retrieve the
pseudodata information that we are interested in. The below is a example
usage:

.. code-block :: python

  from validphys.api import API
  API.recreate_fit_pseudodata(fit="pseudodata_test_fit_n3fit")

If instead we wish to account for the ``postfit`` reshuffling of the replicas which make it through
the postfit selection, we must use the closely related :py:func:`validphys.pseudodata.recreate_pdf_pseudodata`
API method:

.. code-block :: python

  from validphys.api import API
  pseudodata = API.recreate_pdf_pseudodata(fit="pseudodata_test_fit_n3fit")

The return type for both these functions is a `list` of :py:class:`validphys.pseudodata.DataTrValSpec`. Which
is a ``namedtuple`` containing the entire dataset, alongside the training and validation indices:

.. code-block :: python

  >>> type(pseudodata)
  list
  >>> type(pseudodata[0])
  validphys.pseudodata.DataTrValSpec
  >>> replica1 = pseudodata[0]
  >>> replica1_tr = replica1.pseudodata.loc[replica1.tr_idx]
  >>> replica1.pseudodata.loc[replica1.tr_idx]
                              replica 1
  group dataset           id
  NMC   NMC               16   0.336004
                          22   0.349966
                          27   0.385452
                          29   0.361615
                          36   0.430297
  ...                               ...
  ATLAS ATLASZPT8TEVMDIST 56  22.123374
                          59   7.284467
                          61   2.204524
                          62   0.671212
                          63   0.023891

  [223 rows x 1 columns]


.. note::

  When running this action from a runcard, it may be worthwhile to use the ``--parallel`` flag when calling validphys.
  This flag parallelizes dependencies which will compute the pseudodata replicas in an asynchronous manner. This is
  advantageous since the MC replica generation is computationally intensive.
