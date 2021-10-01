.. _howto nf variations:
How to create PDF sets with flavor-number variations
================================================================================

While it is recommended to always perform the PDF fits in an ``nf=5`` flavor
scheme, it is possible to create PDF sets in a variable-flavor-number scheme
where the maximum number of flavors is not equal to five.

After performing a fit in an ``nf=5`` flavor scheme one usually runs
``evolven3fit`` to perform the DGLAP evolution based on the theory that had
also been used during the fit. However, if one instead want to create a PDF set
with a a different maximum number of flavors, the evolution can instead be run
using the ``--theory_id`` flag of ``evolven3fit`` which allows for the
possibility to use a different theory during evolution than was used during the
fit itself.

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

Both, the running of ``evolven3fit`` using a theory id different from the one
used during the fit, and correcting the values for ``AlphaS_MZ`` and ``MZ`` can
be done using the ``varflavors`` script:

.. code-block:: bash

  varflavors [fit folder] [max replicas] [THEORYID] --new_name=[name of output fit]

where ``max replicas`` are the maximum number of replicas evolved by
``evolven3vfit``, ``THEORYID`` is again the theoryID of the theory used to do
the DGLAP evolution, and ``new_name`` is an optional argument that makes the
script copy the original and call it ``new name`` before performing the above
mentioned steps on the ``new name`` fit.

To finalize the PDF set, ``postfit`` should be run in the usual manner.
