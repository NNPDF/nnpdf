.. _datthcomp:

How to do a data theory comparison
==================================

This tutorial explains how to compare the data and theory for a given
data set or list of data sets.

You need to provide:

1. A PDF which includes your data set;
2. A valid theory ID;
3. A choice of cuts policy;
4. A list of data sets to do the comparison for.

Below is an example runcard for a data theory comparison for BCDMSP,
``runcard.yaml``:

.. code:: yaml

   meta:
       title: BCDMSP data/theory comparison
       keywords: [example]
       author: Rosalyn Pearson

   pdfs: 
       - id: NNPDF31_nnlo_as_0118
         label: NNPDF31_nnlo_as_0118

   theoryid: 53

   use_cuts: false

   dataset_inputs:
         - { dataset: BCDMSP}

   template: dthcomparison.md

   actions_:
     - report(main=true)

The corresponding template, ``dthcomparison.md``, looks like this:

.. code::

   %BCDMSP (theory ID 52)

   {@ dataset_inputs plot_fancy @}
   {@ dataset_inputs::pdfs plot_fancy(normalize_to=data)@}
   {@ dataset_inputs::pdfs plot_chi2dist @}
   {@ dataset_inputs::pdfs group_result_table @}

1. ``plot_fancy`` produces data-theory comparison plots for the data.
   This is called twice to produce both normalised and unnormalised sets
   of plots.
2. ``plot_chi2dist`` gives the chi2 distribution between the theory and
   data.
3. ``group_result_table`` gives the numerical values which appear in the
   plots.

Running ``validphys runcard.yaml`` should produce a ``validphys`` report
of the data-theory comparison like the one
`here <https://vp.nnpdf.science/ErmVZEPGT42GCfreWwzalg==/>`__ - see the
`vp-guide <https://data.nnpdf.science/validphys-docs/guide.html#development-installs>`__.
