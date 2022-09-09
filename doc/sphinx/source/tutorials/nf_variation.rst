.. _howto nf variations:

How to create PDF sets with flavor-number variations
================================================================================

The default PDF fits are generated assuming a ``nf=5`` flavor scheme. It is
possible to use a different flavour scheme for the evolution, given an existing
fit at the initial scale.

After :ref:`performing a fit<n3fit-usage>`, the  ``varflavours`` script can be
used to create a PDF evolved with a different maximum number of flavors.
Instead of running ``evolven3fit`` to evolve the fit with the theory used in
the predictions, one would call `varflavours` as follows

.. code-block:: bash

  varflavors fit_folder max_replicas THEORYID --new_name=[name of output fit]

where ``max replicas`` are the maximum number of replicas evolved by
``evolven3vfit``, ``THEORYID`` is again the theoryID of the theory used to do
the DGLAP evolution, and ``new_name`` is an optional argument setting the new
name for the fit.

To finalize the PDF set, :ref:`postfit<postfit>` should be run in the usual manner.

.. note::
   Note that the resulting PDFs will not correspond to
   the best PDF in that flavour scheme (for which we would need to do an actual
   fit yielding different initial scale PDFs).


Implementation details
----------------------

Under the hood, the ``varflavors`` script uses the ``--theory_id`` flag
of ``evolven3fit`` which allows for the possibility to use a different theory
during evolution than was used during the fit itself.

As such, one can perform a fit using a theory corresponding to the ``nf=5``
flavor scheme, and then perform the evolution by running


.. code-block:: bash

  evolven3fit --theory_id=[THEORYID] [configuration folder] [number of replicas]

where ``THEORYID`` corresponds to a theory with a flavor number scheme that
differs from the ``nf=5`` flavor scheme used during the fit.

.. note::
  In case it is desired to create a PDF with flavor-number differing from the
  fit, while also keeping the original ``nf=5`` version, one can use
  ``vp-fitrename`` or ``vp-pdfrename`` to create two PDF sets with different
  names that are based on the same input fit.

If we do this, the values for  ``AlphaS_MZ`` and ``MZ`` in the ``.info`` file
need to be changed to take the ``alphas`` and ``Qref`` values from the theory
used to do the evolution.
