How to implement a new experiment in buildmaster
================================================

Buildmaster is the code that allows the user to generate the ``DATA``
and ``SYSTYPE`` files that contain, respectively, the experimental data
and the information pertaining to the treatment of systematic errors.
Data made available by experimental collaborations comes in a variety of
formats: for use in a fitting code, this data must be converted into a
common format, that contains all the required information for use in PDF
fitting. Such a conversion is realised by the buildmaster code according
to the layout described in :ref:`exp_data_files`.

The user is strongly encouraged to go through that section with care, in
order to familiarise himself with the features of the experimental data,
in general, and the nomenclature of the ``NNPDF`` code, in particular.

To implement a new experiment in buildmaster the first thing to do is to
find the relevant experimental information. As mentioned above, this can
come in a variety of formats. Usually, this is made available from the
`hepdata`_ repository as soon as the corresponding preprint is accepted
for publication. Additional useful resources are the public pages of the
(LHC) experimental collaborations:

- `ATLAS`_
- `CMS`_
- `LHCb`_

A careful reading of the experimental paper is strongly recommended to
understand the information provided, in particular concerning the origin
and the treatment of uncertainties.

Once the details of the experimental measurement are clear, one should
assign the corresponding experiment a name. Such a name must follow the
convention

::

   <name_exp>_<name_obs>_[<extra_info>]

where is the ``<name_exp>`` is name of the experiment in full
(e.g. ATLAS, CMS, LHCB, …), ``<name_obs>``
is the name of the observable (e.g. 1JET, SINGLETOP, TTB, …), and
``[<extra_info>]`` (optional) is a set of strings, separated by underscore, that
encapsulate additional information needed to univocally identify the
measurement (e.g. the c.m. energy, the final state, the luminosity, the
jet radius, …).

The experimental information retrieved from the above must be collected
(ideally with minimal editing and in plain text format) in a new
directory

::

   buildmaster/rawdata/<name_exp>_<name_obs>_[<extra_info>]

A metadata file has to be created in the ``.yaml`` format as

::

   buildmaster/meta/<name_exp>_<name_obs>_[<extra_info>].yaml

with the following structure

.. code-block:: yaml

   ndata:    <number of datapoints>
   nsys:     <number of systematic errors>
   setname:  <setname in double quotes, i.e. "<name_exp>_<name_obs>_[<extra_info>]">
   proctype: <process type> in double quotes)

A list of the available process types can be found at :ref:`process_type_label`.
If the process type corresponding to the experiment under consideration
is not contained in that list, a new process type should be defined and
implemented.

Then the user has to create the header for a new class with the dataset
name in

::

   /buildmaster/inc/<name_exp>_<name_obs>_[<extra_info>].h

as follows

::

   class MY_NEW_DATASET_CLASSFilter: public CommonData {
   public: MY_NEW_DATASET_CLASSFilter("MY_NEW_DATASET_NAME") { ReadData(); }
   private:
       void ReadData();
   }

and implement the ``ReadData()`` function in

::

   /buildmaster/filter/<name_exp>_<name_obs>_[<extra_info>].cc

Such a function should read from the rawdata file

- the kinematic variables required for the specific process under consideration:
  ``fKin1``, ``fKin2``, ``fKin3``
- the data: ``fData``
- the statistical uncertainty: ``fStat``
- the systematic uncertainties: ``fSys``

Important remarks.

1. The relevant information regarding uncertainty correlations must be
   consistently implemented. Depending on the specific experiment one is
   considering, this may be provided either as a full breakdown of
   correlated systematics or through a covariance (or correlation)
   matrix. In the latter case, if the dataset is made by ``N`` data,
   ``N`` systematics have to be produced from the decomposition of the
   covariance matrix, using the function ``genArtSys`` (in
   ``buildmaster/src/buildmaster_utils.cc``). Sometimes a covariance
   matrix is provided also for the statistical uncertainties. In such
   cases the ``fStat`` variable should be set to zero, and the
   statistical uncertainty should be implemented as a set of ``N``
   additional artificial systematics obtained from the decomposition of
   the systematic covariance matrix through ``genArtSys``.

2. Uncertainties are sometimes provided as sets of independent (left and
   right) asymmetric values. They are usually estimated, data point by
   data point, by varying upwards and downwards the nuisance parameters
   in the experimental model used for their determination. Note that an
   upwards (downwards) variation of the nuisance parameters does not
   necessarily generate a positive (negative) variation of the data
   point expectation value. Therefore, left and right uncertainties can
   be both positive, both negative, positive and negative, or negative
   and positive. If the left uncertainty is negative and the right
   uncertainty is positive (i.e. a downwards shift of the nuisance
   parameters generates a decrease of the data point expectation value
   and an upwards shift of the nuisance parameters generates an increase
   of the data point expectation value), they can be symmetrised using
   the `D’Agostini`_ rule, as implemented in the ``symmetriseErrors``
   function (in ``buildmaster/src/buildmaster_utils.cc``). The data
   point expectation value should be shifted accordingly. If the signs
   of the left and right asymmetric uncertainties are mixed, other
   prescriptions (to preserve correlations/anticorrelations) must be
   adopted, see the implement

3. Consider testing that the additive and multiplicative columns of the
   commondata are self-consistent. The multiplicative columns should be related
   to the additive columns (schematically) by
   ``add_columns = mult_columns * central_values * 1e-2``. The easiest
   way to test this is to add the newly implemented dataset to the list
   of datasets tested in :py:mod:`validphys.tests.test_commondata_columns`.
   If you commit this change to the repo then the :ref:`CI <CI>` will always check this is
   the case, in case somebody edits the dataset in the future.

.. _hepdata: https://www.hepdata.net/
.. _ATLAS: https://twiki.cern.ch/twiki/bin/view/AtlasPublic
.. _CMS: https://home.cern/news
.. _LHCb: http://lhcbproject.web.cern.ch/lhcbproject/Publications/LHCbProjectPublic/Summary_all.html
.. _D’Agostini: https://arxiv.org/abs/physics/0403086
