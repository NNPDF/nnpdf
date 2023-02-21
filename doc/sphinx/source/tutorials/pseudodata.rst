.. _pseudodata:
Obtaining the pseudodata used by an ``n3fit`` fit
=================================================

Since version 4.0.7 the Monte Carlo data replicas are saved to disk by default as the fit is performed. 
This can be deactivated by setting the ``savepseudodata`` flag to ``False`` under the ``fitting`` namespace in the fit runcard:

.. code-block :: yaml

   fitting:
    savepseudodata: False

If the ``savepseudodata`` flag is not set to ``False``, the training and validation splits to disk 
under files named ``datacuts_theory_fitting_training_pseudodata.csv`` and similarly for the validation 
split. These can then be loaded within validphys by leveraging the 
:py:func:`validphys.pseudodata.read_fit_pseudodata` action:

.. code-block :: python

   >>> from validphys.api import API
   >>> pseudodata = API.read_fit_pseudodata(fit="pseudodata_test_fit_n3fit")
   >>> replica1_info = pseudodata[0]
   >>> replica1_info.pseudodata.loc[replica1_info.tr_idx]
                                  replica 1
  group dataset           id
  ATLAS ATLASZPT8TEVMDIST 1    29.856281
                          3    14.686290
                          4     8.568288
                          5     2.848544
                          6     0.704977
  ...                                ...
  NMC   NMCPD             247   0.688019
                          249   0.713272
                          255   0.673997
                          256   0.751973
                          259   0.750572

  [223 rows x 1 columns]

With the postfit reshuffling handled instead by :py:func:`validphys.pseudodata.read_pdf_pseudodata`.


Reconstructing pseudodata
--------------------------

.. warning::

  The functionality described here is not guaranteed to work between different versions of the code
  or its dependencies. Specifically, if anything breaks the pseudodata generation between commits, 
  e.g. changes to the theory predictions or settings or the random number generator, it is **not 
  possible to reconstruct previously generated pseudodata** for the code state at such different 
  commits.

Suppose one has obtained a fit using the ``n3fit`` framework and wants to do some analysis that requires
knowing exactly the data that the neural networks saw during the fitting procedure while this has not been stored.
The information is reproducible given the various seeds in the fit runcard.

The 3 random seeds used in the fit are ``trvlseed`` which determines the training/validation splitting, ``nnseed``
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