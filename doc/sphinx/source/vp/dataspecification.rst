.. _data_specification:

Data specification
==================

.. _datasetspec-core-label:

``DataSetSpec`` - Core dataset object
-------------------------------------

The declaration of dataset specifications within the ``validphys`` framework
is handled using the ``dataset_input``  key, or ``dataset_inputs`` namespace
list for a collection of datasets. Through this keyword the user
is provided a granular degree of customizability for each dataset considered
in the runcard; in particular the handling of K-factors, systematic uncertainties,
training fraction, or dataset weight can be modified in this declaration. Moreover,
the ``metadata_group`` keyword allows for a flexible grouping of a collection of
datasets, organizing them into disjoint subsets depending on, for example, the
experiment class to which they belong or their process type.

The core dataset object in ``validphys`` is the
:py:mod:`validphys.core.DataSetSpec` which is responsible for loading the
dataset, covariance matrix and applying cuts.

Specifying a dataset
~~~~~~~~~~~~~~~~~~~~

In a validphys runcard the settings for a single dataset are specified using a
``dataset_input``. This is a dictionary which minimally specifies the name of
the dataset, but can also control behaviour such as contributions to the
covariance matrix for the dataset and C-factors.

Here is an example dataset input:

.. code:: yaml

    dataset_input:
        dataset: CMSZDIFF12
        cfac: [QCD,NRM]
        sys: 10

This particular example is for the ``CMSZDIFF12`` dataset, the user has
specified to use some C-factors given by ``cfac`` as well as ``sys: 10``, which
corresponds to an additonal contribution to the covariance matrix accounting for
statistical fluctuations in the C-factors. These settings correspond to NNLO
predictions and so presumably elsewhere in the runcard the user would have
specified an NNLO theory - such as theory 53.

We can use the API to return an instance of ``DataSetSpec`` in a development
environment using the settings above

.. code:: python

    >>> from validphys.api import API
    >>> ds_spec = API.dataset(
    ...     dataset_input={"dataset": "CMSZDIFF12", "cfac": ["QCD", "NRM"], "sys": 10},
    ...     use_cuts="internal",
    ...     theoryid=53
    ... )
    >>> type(ds_spec)
    <class 'validphys.core.DataSetSpec'>

Here we are obtaining the result from the production rule
:py:mod:`validphys.config.CoreConfig.produce_dataset`, the required arguments
are ``dataset_input``, ``cuts`` and ``theoryid``.

.. note::

    It seems odd to require theory settings such as a ``theoryid`` in the
    ``dataset_input`` in order to load data. However, this is a relic of the
    underlying C++ code that performs the loading of data, which intrinsically
    groups together the commondata (CSVs containing data central values and
    uncertainties) and :ref:`fktables`.

    Clearly there is a big margin for error when manually entering
    ``dataset_input`` and so there is a
    `project <https://github.com/NNPDF/nnpdf/issues/226>`_ that aims to have a
    stable way of filling many of these settings with correct default values.

The ``DataSetSpec`` contains all of the information used to construct it, e.g.

.. code:: python

    >>> ds_spec.thspec
    TheoryIDSpec(id=53, path=PosixPath('/Users/michael/conda/envs/nnpdf-dev/share/NNPDF/data/theory_53'))
    >>> ds_spec.name
    'CMSZDIFF12'

but also importantly has a ``load`` method, which returns an instance of the
``DataSet`` that is generated from the C++ code using SWIG. This new object
contains numpy arrays of data central values and experimental covariance
matrices, e.g:

.. code:: python

    >>> ds_libnnpdf = ds_spec.load()
    >>> ds_libnnpdf.get_cv() # get central values of dataset
    array([2917.  , 1074.  ,  460.5 ,  222.6 ,  109.8 ,   61.84,   30.19,
           2863.  , 1047.  ,  446.1 ,  214.5 ,  110.  ,   58.13,   29.85,
           2588.  ,  935.5 ,  416.3 ,  199.  ,  103.1 ,   54.06,   28.45,
           1933.  ,  719.5 ,  320.7 ,  161.1 ,   84.62,   47.57,   24.13])

In practice, actions that require experimental data and/or covariance matrices
will make use of the :py:mod:`validphys.results.results` provider, which is a
tuple of :py:mod:`validphys.results.DataResult` and
:py:mod:`validphys.results.ThPredictionsResult`. Since in this case we are
additionally generating theory predictions, we are additionally required to
specify a PDF

.. code:: python

    >>> results = API.results(
    ...     dataset_input={"dataset": "CMSZDIFF12", "cfac": ["QCD", "NRM"], "sys": 10},
    ...     use_cuts="internal",
    ...     theoryid=53,
    ...     pdf="NNPDF31_nnlo_as_0118"
    ... )
    PDF: NNPDF31_nnlo_as_0118  ErrorType: Monte Carlo booked
    LHAPDF 6.2.3 loading all 101 PDFs in set NNPDF31_nnlo_as_0118
    NNPDF31_nnlo_as_0118, version 1; 101 PDF members
    NNPDF31_nnlo_as_0118 Initialised with 100 members and errorType replicas
    >>> results
    (<validphys.results.DataResult object at 0x1518528350>, <validphys.results.ThPredictionsResult object at 0x1a19a4da50>)

The covariance matrix associated with the ``DataResult`` in this tuple was
constructed by :py:mod:`validphys.results.covmat`, which allows the user to
change the behaviour of the covariance matrix - such as adding theory
uncertainties computed from scale variations or using a t0 PDF to calculate the
multiplicative contributions to the covariance matrix - for more detail see
:py:mod:`validphys.results.covmat` itself.

``DataGroupSpec`` - core object for multiple datasets
-----------------------------------------------------

The core object for multiple datasets is :py:mod:`validphys.core.DataGroupSpec`,
which is similar in many regards to the `DataSetSpec`, but additionally handles
the loading of multiple datasets. In particular, when constructing the
covariance matrix, it takes into account any uncertainties which are correlated
across the different datasets.

Specifying multiple datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Multiple datasets are specified using the ``dataset_inputs`` key, which is a
list in which each element of the list is a valid ``dataset_input``. For
example:

.. code:: yaml

    dataset_inputs:
        - { dataset: NMC }
        - { dataset: ATLASTTBARTOT, cfac: [QCD] }
        - { dataset: CMSZDIFF12, cfac: [QCD,NRM], sys: 10 }

We see that multiple datasets are inputted as a flat list, i.e. that there is no
hierarchy to the datasets which splits them into experiments or process types.
The grouping of datasets is done internally according to the metadata of
datasets and is controlled by the ``metadata_group`` key. This can be any key
which is present in the ``PLOTTING`` file of each dataset - for example
``experiment`` or ``nnpdf31_process``. The default value for ``metadata_group``
is ``experiment``. Other groupings might be relevant, for example when
contructing a theory covariance matrix, in which case you want to group datasets
according to process type rather than experiment. The grouping is performed by
the production rule
:py:mod:`validphys.config.CoreConfig.produce_group_dataset_inputs_by_metadata`,
which returns a list with length equal to the number of distinct groups. Each
element is a namespace with the ``group_name`` and list of ``dataset_input`` s
for that specific group e.g:

.. code:: python

    >>> API.group_dataset_inputs_by_metadata(
    ...    dataset_inputs=[
    ...        {"dataset":"NMC"},
    ...        {"dataset": "ATLASTTBARTOT", "cfac": ["QCD"]},
    ...        {"dataset": "CMSZDIFF12", "cfac": ["QCD","NRM"], "sys": 10 }],
    ...    metadata_group="experiment"
    ... )
    [
        {'data_input': [DataSetInput(name='NMC', sys=None, cfac=(), frac=1, weight=1)], 'group_name': 'NMC'},
        {'data_input': [DataSetInput(name='ATLASTTBARTOT', sys=None, cfac=['QCD'], frac=1, weight=1)], 'group_name': 'ATLAS'},
        {'data_input': [DataSetInput(name='CMSZDIFF12', sys=10, cfac=['QCD', 'NRM'], frac=1, weight=1)], 'group_name': 'CMS'}
    ]

Here we see that the namespace key is ``data_input`` rather than
``dataset_inputs``, which is because ``data_input`` bridges the gap between the
current way of specifying data (with ``dataset_inputs``) and a deprecated
specification using the ``experiments`` key. The production rule that returns a
``DataGroupSpec`` is :py:mod:`validphys.config.CoreConfig.produce_data` through
the following pipeline

.. code::

    dataset_inputs or experiments -> data_input -> data

For example, the following runcard produces a single-column table with a row
containing the 𝞆² of the specificed datasets, grouped by ``experiment``:

.. code:: yaml

    dataset_inputs:
        - { dataset: NMC }
        - { dataset: ATLASTTBARTOT, cfac: [QCD] }
        - { dataset: CMSZDIFF12, cfac: [QCD,NRM], sys: 10 }

    theoryid: 53

    dataspecs:
     - pdf: NNPDF31_nnlo_as_0118
     - speclabel: "3.1 NNLO"

    use_cuts: internal

    actions_:
     - dataspecs_groups_chi2_table

If we specify a ``metadata_group`` in the runcard, like so

.. code:: yaml

    metadata_group: nnpdf31_process

    dataset_inputs:
        - { dataset: NMC }
        - { dataset: ATLASTTBARTOT, cfac: [QCD] }
        - { dataset: CMSZDIFF12, cfac: [QCD,NRM], sys: 10 }

    theoryid: 53

    dataspecs:
     - pdf: NNPDF31_nnlo_as_0118
       speclabel: "3.1 NNLO"

    use_cuts: internal

    actions_:
     - dataspecs_groups_chi2_table

then we instead get a single-column table, but with the datasets grouped by
process type, according the
`theory uncertainties paper <https://arxiv.org/abs/1906.10698>`__.


Note that actions which rely on grouping use a fallback value of
``metadata_group`` which gets set in the production rule for
``processed_metadata_group``. It may be useful to use the
namespace key ``processed_metadata_group`` in reports and
actions alike to make use of this. Here is an example giving
sensible titles/section headings e.g.:

.. code:: yaml

    template_text: |
     # chi2 grouped by {processed_metadata_group}
     {@dataspecs_groups_chi2_table@}

    actions_:
     - report(main=True)

Custom grouping
---------------

It is possible to define a custom grouping at the level of the runcard, which
is useful for temporary groupings or testing out a new group which may
eventually be added the the metadata. The user can use custom groupings
by setting ``metadata_group=custom_group`` in the runcard and then adding
the ``custom_group`` key to each dataset_input as follows

.. code:: yaml

  metadata_group: custom_group

  dataset_inputs:
    - { dataset: NMC, custom_group: traca }
    - { dataset: NMCPD, custom_group: traco }
    - { dataset: LHCBWZMU7TEV, cfac: [NRM], custom_group: pepe }
    - { dataset: LHCBWZMU8TEV, cfac: [NRM], custom_group: pepa }
    - { dataset: ATLASWZRAP36PB}

Note that we didn't set any group for ``ATLASWZRAP36PB``, but that's ok: any
datasets which are not explicitly given a ``custom_group`` get put into the
``unset`` group.

For more information on how to immortalise your custom grouping in the
metadata and call that grouping as in the previous examples
(i.e with ``nnpdf31_process``) see :ref:`add_special_label`.

Action naming conventions
-------------------------

There are some general rules that should be observed when adding new actions to
``validphys``. Firstly, try to indicate the required runcard input for an
action in the name of the function. Take for example the provider
``dataset_inputs_results``. The returned object is a ``results`` object: a
tuple of data and theory predictions, which is used by a wide range of other
actions, notably when calculating a 𝞆². The first part of the name
``dataset_inputs`` refers to the runcard input required to process the action.
This is especially useful for actions which use a group of datasets or
``data``, because the dependency tree for these actions is not neccessarily
obvious to somebody who is unfamiliar with the code. As explained above,
``dataset_inputs -> data_input -> data`` and so the action name serves to guide
the user to creating a working runcard as easily as possible.

The second general rule is that if your action makes use of ``collect``
somewhere in the dependency graph, then consider prepending what is collected
over to the action name. For example: ``dataspecs_groups_chi2_table``, which
depends on

.. code:: python

    dataspecs_groups_chi2_data = collect("groups_chi2", ("dataspecs",))

and in turn

.. code:: python

    groups_chi2 = collect("dataset_inputs_abs_chi2_data", ("group_dataset_inputs_by_metadata",))

Without having to find these specific lines in the code we would be able to
guess that the 𝞆² is collected first over groups of data (``groups_chi2``), and
then over ``dataspecs``. Naming functions according to these rules helps make
the general workings of the underlying code more transparent to an end user.

.. _backwards-compatibility:

Backwards compatibility
-----------------------

Where possible, backwards compatibility with runcards which use the
``experiments`` key has been preserved. For example, with the
``dataspecs_groups_chi2_table`` example above we could also use the following
input

.. code:: yaml

    experiments:
     - experiment: NMC
       datasets:
        - { dataset: NMC }
     - experiment: ATLAS
       datasets:
        - { dataset: ATLASTTBARTOT, cfac: [QCD] }
     - experiment: CMS
       datasets:
        - { dataset: CMSZDIFF12, cfac: [QCD,NRM], sys: 10 }

    theoryid: 53

    dataspecs:
     - pdf: NNPDF31_nnlo_as_0118
       speclabel: "3.1 NNLO"

    use_cuts: internal

    actions_:
     - dataspecs_groups_chi2_table

The user should be aware, however, that any grouping introduced in this way is
purely superficial and will be ignored in favour of the experiments defined by
the metadata of the datasets.

Runcards that request actions that have been renamed will not work anymore.
Generally, actions that were previously named ``experiments_*`` have been
renamed to highlight the fact that they work with more general groupings.

If you are writing a runcard whereby you want to take the data from a ``fit``,
and either do not know whether the fit uses the new or old data specification or
require the runcard to be agnostic to the data specification in the fit,
there are a couple of options.

First and foremost try using the ``fitinputcontext`` production rule to extract
the data from the fit. This production rule handles both styles of runcard
out of the box:

.. code:: yaml

    metadata_group: nnpdf31_process

    fit: NNPDF31_nnlo_as_0118_DISonly

    dataspecs:
     - pdf: NNPDF31_nnlo_as_0118
       speclabel: "3.1 NNLO"

    use_cuts: internal

    actions_:
     - fitinputcontext dataspecs_groups_chi2_table

The production rule sets the ``theoryid`` and ``data_input`` based on the
runcard for the specified ``fit``. Note that you can also use ``fitcontext``
which does all of the above, and additionally sets the ``pdf`` to be the
``fitpdf``.

In many cases where an action is prefixed with ``dataspecs``, indicating that
a table or plot will contain some results collected over the ``dataspecs``,
there will be a similar action prefixed with ``fits``, where instead the
results in the table or plot will have been collected over ``fits`` with
``fitcontext`` taken into account.

.. warning::
  Whilst it is possible to specify ``data_input: {from_: fitinputcontext}``
  directly in the runcard, it is highly recommended **not** to do this where
  possible. Instead take ``dataset_inputs`` directly ``from_: fit``
  irrespective of whether the fit uses new or old data specification; since
  the conversion from the old style data specification is handled internally
  using :py:func:`validphys.utils.experiments_to_dataset_inputs` in
  conjunction with :py:meth:`validphys.core.FitSpec.as_input`.  (See below for
  a detailed explanation).

Currently the ``pseudodata`` and ``chi2grids`` modules have not been updated to
use ``dataset_inputs`` and so require ``experiments`` to be specified in the
runcard.

.. seealso:: Why not to use ``data_input: {from_: fitinputcontext}``?

  Taking a key ``from_`` a production rule causes that key to be
  overwritten in inner namespaces. The grouping function essentially returns a
  namespace list with each item in the list specifying a different namespace,
  where ``data_input`` is defined as the datasets within that group. If
  the user specifies ``data_input: {from_: fitinputcontext}`` in the runcard,
  the inner ``data_input`` for each group will be overwritten and instead each
  group will contain all of the datasets from the fit - which is incorrect.
  This is regarded as a bug, the relevant issue is:
  https://github.com/NNPDF/reportengine/issues/38

What do I need to change in my runcards?
-----------------------------------------------------

Efforts have been made to ensure a degree of backwards compatibility, however
there are two main things which may need to be changed in old runcards.

1. For theorycovariance runcards, you must add a line with
``metadata_group: nnpdf31_process``, or else the prescriptions for scale
variations will not vary scales coherently for data
within the same process type, as usually desired, but rather for data within
the same experiment.

2. Many actions which were based on experiments have changed names as they are
now based on arbitrary groupings given by ``metadata_group``. The table below gives
the old name alongside the new one. These need to be updated for the runcards to 
continue to work.

.. list-table:: Updated names for old actions
   :widths: 25 25
   :header-rows: 1

   * - Old name
     - New name
   * - plot_fits_experiments_phi
     - plot_fits_groups_data_phi
   * - plot_phi_experiment_dist
     - plot_dataset_inputs_phi_dist
   * - plot_experiments_chi2
     - plot_groups_data_chi2
   * - plot_fits_experiments_chi2
     - plot_fits_groups_data_chi2
   * - experiment_result_table
     - group_result_table
   * - experiment_result_table_68cl
     - group_result_table_68cl
   * - experiments_covmat
     - groups_covmat
   * - experiments_sqrtcovmat
     - groups_sqrtcovmat
   * - experiments_invcovmat
     - groups_invcovmat
   * - experiments_normcovmat
     - groups_normcovmat
   * - experiments_corrmat
     - groups_corrmat
   * - experiment_results
     - dataset_inputs_results
   * - experiments_chi2_table
     - groups_chi2_table
   * - fits_experiments_chi2_table
     - fits_groups_chi2_table
   * - fits_experiments_phi_table
     - fits_groups_phi_table
   * - dataspecs_experiments_chi2_table
     - dataspecs_groups_chi2_table
   * - experiments_central_values
     - groups_central_values
   * - experiments_chi2_table_theory
     - groups_chi2_table_theory
   * - experiments_chi2_table_diagtheory
     - groups_chi2_table_diagtheory
