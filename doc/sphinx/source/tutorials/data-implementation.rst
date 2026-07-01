.. _data-implementation:

How to Implement New Data
=========================

This tutorial covers the steps and some points that need to
be taken into account when implementing new datasets into the
``CommonData`` format (which you can read more about :ref:`here <commondata>`).

Dataset implementation refers to transfering the information from
experimental datasets into the format required for PDF fitting.
You can read more about data storage and handling in
:ref:`this section <data-storage-th-pred>`.


Step 1: Identify what to implement and download your data
---------------------------------------------------------

Datasets to implement are usually associated to a corresponding
release paper. You should identify which figures or what sections
of the paper describe the observables and experimental results to
be implemented.

The raw data is usually available in the `HEPData website <https://www.hepdata.net/>`_.
You should create a new direcory for your dataset in::

    nnpdf/nnpdf_data/nnpdf_data/commondata/

where the name of your new directory should be equivalent to that
of your dataset, and should follow the :ref:`naming convention <dataset-naming-convention>`.

You can then create a directory called ``rawdata`` within your dataset
directory and download the HEPData entries corresponding to your observable
as ``.yaml`` files in this new directory.


Step 2: Create a ``metadata.yaml`` file
---------------------------------------

In your dataset directory, create a file called ``metadata.yaml``, and
populate it with the :ref:`relevant information <metadata-spec>`. 

**Pro tip:** just copy a metadata file from a previously implemented
dataset with the same process and/or experiment and change all the entries
to correspond to your dataset.


Step 3: Filter the data
-----------------------

You are now ready to filter your raw data so that the information is broken
down into the following files:

- ``data_<observable>.yaml``: central values of data points.
- ``kinematics_<observable>.yaml``: the kinematic values for each bin of the data, eg. pseudorapidity, transverse energy...
- ``uncertainties_<observable>.yaml``: list of all the uncertainties for each data point.

Note that, in some cases, you will only have one observable per dataset, so
you won't need to specify the observable in your file names. 

To filter your data, you should create a ``filter.py`` script. The best way to do this
is to have a look at other such scripts in other previously implemented datasets,
ideally for one made for the same or a similar process and/or observable.

You should also save your ``filter.py`` script in your dataset directory for review
and as a record. 

Things to watch out for:

- Some data points may have asymmetric uncertainties. Your filter script should symmetrise these points by taking the average between the maximum and minimum uncertainties and shifting the central value accordingly.
- The units of implemented data should be *fb*. If your data is in a different unit, your filter script should convert it accordingly.

Step 4: Classify the uncertainties
----------------------------------

The ``filter.py`` script will assign a ``treatment`` and ``type`` to each uncertainty.
For example, the first lines of your script might look like this:

.. code-block::

    definitions:
        stat:
            description: Uncorrelated statistical uncertainties
            treatment: ADD
            type: UNCORR
        sys_corr_1:
            description: sysPhotonID
            treatment: ADD
            type: UNCORR
        sys_corr_2:
            description: RESOLUTION_AF2
            treatment: MULT
            type: CORR
        sys_corr_3:
            description: RESOLUTION_MATERIALCALO
            treatment: MULT
            type: CORR
        sys_corr_4:
            description: RESOLUTION_MATERIALCRYO
            treatment: MULT
            type: CORR
        sys_corr_5:
            description: RESOLUTION_MATERIALGAP
            treatment: MULT
            type: CORR

You should check one by one to make sure that the ``treatment`` and ``type`` for
every source of uncertainty is correct. The paper where your dataset was presented
should describe the nature of each uncertainty.

If the information on the treatment or type is not available, here are some rules
of thumb that may assist you, baring in mind that they don't always work and should
be applied on a case-by-case basis:

- Statistical uncertainties are usually **uncorrelated**,
- Uncertainties are most likely to be **multiplicative** and **uncorrelated**.