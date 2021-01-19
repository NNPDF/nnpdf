How to analyse a closure test
=============================

There are two ``validphys`` reports which are used to analyse closure
tests. The report which compares two different closure tests is
generated using the ``vp-comparefits`` tool. The individual report which
compares a closure test to the PDF used as the underlying law (from which
pseudodata is generated) can be found in the main ``nnpdf`` repository under
``validphys2/examples/closure_templates/``.

Comparison to underlying PDF
----------------------------

This report is very similar to the full validphys fit report produced by
``vp-comparefits``, the main differences are: only the closure fit is
specified and it is compared to the PDF specified under the key
``fakepdf`` in the fit runcard; and the report contains some additional
sections dedicated to the closure test statistical estimators.

These are:

-  bias

   The :math:`\chi^2` between the central predictions of the closure
   test and the underlying PDF predictions, a measure of how well the
   underlying predictions are reproduced by the closure test
-  variance

   The mean across replicas of the :math:`\chi^2` between predictions
   from a single replica and the central replica. Can be
   thought of as the variance of the replica predictions normalised by
   the covariance.

-  :math:`\Delta_{\chi^2}`

   :math:`\chi^2` between central predictions and level 1 replica minus
   the :math:`\chi^2` between level 1 data and underlying PDF, gives
   some information on the direction the central predictions are being
   pulled with respect the the level 1 data

Detailed information on the closure test estimators can be found
`here <https://www.wiki.ed.ac.uk/display/nnpdfwiki/n3fit+closure+results?preview=/407535993/418483968/Statistics.pdf>`__.

Comparison to another closure test
----------------------------------

To compare to another closure test the user can use a special command
line option with ``vp-comparefits``, running

::

    vp-comparefits --closure <other options>

will cause ``vp-comparefits`` to use a special closure comparison report
which compares to closures for the statistical estimators above.

care should be taken that both of the closure tests being compared were fitting
pseudodata which had the same level 1 shift applied to the same underlying law.
Clearly the closure test estimators have some underlying dependence on these.
To ensure this, check that both fits have a matching ``fakepdf`` key and
``filterseed``, additionally the fits should have been ran during a time window
within which the data generation algorithm was not modified.
