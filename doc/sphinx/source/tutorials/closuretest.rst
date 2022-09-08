.. code:: {eval-rst}

   .. _tut_closure:

How to run a closure test
=========================

Closure tests are a way to validate methodology by fitting on pseudodata
generated from pre-existing PDFs. There are different levels of closure
tests which aim to validate different components of the fitting
toolchain.

Brief background
----------------

For more detailed information on the conception of closure tests, see
the `NNPDF3.0 paper <https://arxiv.org/abs/1410.8849>`__.

Each closure test defines a ``fakepdf`` in the runcard, which will be
referred to here as the underlying law. For the purpose of the closure
test it can be thought of as being a proxy for the true PDF.

There are three levels of closure test:

1. level 0

   -  central pseudodata is given by central predictions of the
      underlying law
   -  no MC noise is added on top of the central data, each replica is
      fitting the same set of data

2. level 1

   -  central pseudodata is shifted by some noise η which is drawn from
      the experimental covariance matrix and represents ‘real’ central
      values provided by experimentalists which do not sit exactly on
      the underlying law but are consistent with it according to their
      own uncertainty
   -  no MC noise is added, each replica fits a subset of the same
      shifted data. There is however a difference in the
      training/validation split used for stopping, the spread on
      replicas can be thought of as the spread due to this split in
      addition to any methodological uncertainty.

3. level 2

   -  central pseudodata is shifted by level 1 noise η
   -  MC noise is added on top of the level 1 shift
   -  level 2 is a proxy of a real fit, where the underlying law is
      known

The advantage of knowing the underlying law is that we can see how well
the methodology is extracting this from shifted data, using closure test
estimators.

The main obvious disadvantage is that a pre-existing PDF may not be a
suitable proxy for the underlying law.

Preparing the closure test runcard
----------------------------------

To run a closure test we require a standard fit runcard. The main
section which controls closure test specific behaviour can be found
under ``closuretest``. Before you’ve made any changes, a typical
``closuretest`` section will be as follows:

.. code:: yaml

   closuretest:
     filterseed  : 0   # Random seed to be used in filtering data partitions
     fakedata    : False
     fakepdf     : MMHT2014nnlo68cl
     errorsize   : 1.0 # uncertainties rescaling
     fakenoise   : False
     rancutprob  : 1.0 # Fraction of data to be included in the fit
     rancutmethod: 0   # Method to select rancutprob data fraction
     rancuttrnval: False # 0(1) to output training(valiation) chi2 in report
     printpdf4gen: False # To print info on PDFs during minimization

Setting ``fakedata`` to ``True`` will cause closure test pseudodata to
be generated and subsequently fitted. The PDf which the pseudodata will
be generated from is specified by the ``fakepdf`` key. It is strongly
advised to set the ``fakepdf`` and ``t0pdfset``, found under
``datacuts`` to be the same PDF, unless specifically testing the impact
of the t0 procedure.

The ``fakenoise`` key specifies whether or not the level 1 shift η will
be add to the pseudodata during the filtering step, this is require for
**both** level 1 and level 2 closure tests.

An example of a typical level 1 or level 2 ``closuretest`` specification
is given

.. code:: yaml

   closuretest:
     filterseed  : 0   # Random seed to be used in filtering data partitions
     fakedata    : True
     fakepdf     : MMHT2014nnlo68cl
     errorsize   : 1.0 # uncertainties rescaling
     fakenoise   : True
     rancutprob  : 1.0 # Fraction of data to be included in the fit
     rancutmethod: 0   # Method to select rancutprob data fraction
     rancuttrnval: False # 0(1) to output training(valiation) chi2 in report
     printpdf4gen: False # To print info on PDFs during minimization

Note that it is *critical* that two closure tests which are to be
compared have the same ``filterseed``. They should also both have been
run during a time where no major changes were made to data generation.
This is because fits with different level 1 noise produce different
closure test estimators. See for example a
`report <https://vp.nnpdf.science/mbcTUd6-TQmQFvaGd37bkg==/>`__
comparing two level 2 closure tests with identical settings apart from
``filterseed``. Note that setting ``filterseed`` to 0 will use the
default seed for the selected ``fitting::rngalgo``.

There are still some relevant settings to the closure test. For the
above example we would choose that the t0 set was the same as the
underlying law:

.. code:: yaml

   datacuts:
     t0pdfset     : MMHT2014nnlo68cl # PDF set to generate t0 covmat
     ...

Finally we need to specify whether or not MC replicas will be generated
in the fit, differentiating between a level 1 and level 2 closure test.
This can be achieved by setting ``genrep`` under ``fitting`` to be
``True``

.. code:: yaml

   fitting:
     ...
     genrep   : True
     ...

Summary for each level of closure test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See below for the keys which specify each level of closure test, other
keys can be chosen by the user.

Level 0
^^^^^^^

.. code:: yaml

   fitting:
     ...
     genrep   : False
     ...
   closuretest:
     ...
     fakedata    : True
     fakenoise   : False
     ...

Level 1
^^^^^^^

.. code:: yaml

   fitting:
     ...
     genrep   : False
     ...
   closuretest:
     ...
     fakedata    : True
     fakenoise   : True
     ...

Level 2
^^^^^^^

.. code:: yaml

   fitting:
     ...
     genrep   : True
     ...
   closuretest:
     ...
     fakedata    : True
     fakenoise   : True
     ...

Running a closure test with ``n3fit``
-------------------------------------

Running a closure test with ``n3fit`` will require a valid ``n3fit``
runcard, with the closure test settings modified as shown
`above <#preparing-the-closure-test-runcard>`__. The difference between
running a closure fit in ``n3fit`` and a standard fit is that the user
is required to run ``vp-setupfit`` on the runcard before running
``n3fit``. This is because the filtering of the data is required to
generate the pseudodata central values. The filtered data should then be
rebuilt before the fit, so there is no risk of the fit crashing due to
multiple replicas rebuilding the data simultaneously. The workflow is as
follows:

.. code:: bash

   $ vp-setupfit fitname.yml
   $ vp-rebuild-data fitname
   $ n3fit fitname.yml <replica_number>

You will still need to evolve the fit and run ``postfit`` as with a
standard ```n3fit`` <../tutorials/run-fit.md>`__
