.. _datthcomp:

How to do a data theory comparison
==================================

This tutorial explains how to compare the data and theory for a given data set or list of data sets.

You need to provide:

1. A PDF which includes your data set;
2. A valid theory ID;
3. A choice of cuts policy;
4. A list of data sets to do the comparison for;
5. Options to shift theoretical predictions according to the correlated part of the experimental uncertainties and/or to normalise the comparison to the central value of the experimental data.

Below is an example runcard for a data theory comparison for BCDMSP, ``runcard.yaml``:

.. code:: yaml

  meta:
      title: BCDMSP data/theory comparison
      keywords: [example]
      author: Rosalyn Pearson

  dataspecs:
    - speclabel: "NNPDF40 (w/o shift)"
      theoryid: 40_000_000
      use_cuts: "internal"
      with_shift: False
      pdf: NNPDF40_nnlo_as_01180
    - speclabel: "NNPDF40 (w/ shift)"
      theoryid: 40_000_000
      use_cuts: "internal"
      with_shift: False
      pdf: NNPDF40_nnlo_as_01180 

  dataset_inputs:
      - { dataset: HERA_NC_318GEV_EAVG_BOTTOM-SIGMARED, variant: legacy}
      - { dataset: ATLAS_1JET_8TEV_R06, variant: legacy}
      - { dataset: BCDMS_NC_NOTFIXED_P_EM-F2, variant: legacy}

  template_text: |
     # Data theory comparison with and without shifts
     {@ with dataset_inputs @}
     {@ plot_fancy_dataspecs(normalize_to=data) @}
     {@ endwith @}

  actions_:
    - report(main=true)

The function ``plot_fancy_dataspecs`` produces data-theory comparison plots for the specified list of data for all of the data specifications ``dataspecs``. The code can be run as :code:`validphys runcard.yaml` which will produce a ``validphys`` report with the desired plots. See the runcard ``data_theory_comparison.yaml`` in the validphys ``examples`` folder for details.
