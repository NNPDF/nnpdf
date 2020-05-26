Data specification
==================

``DataSetSpec`` - Core dataset object
-------------------------------------

The core dataset object in ``validphys`` is the :py:mod:`validphys.core.DataSetSpec`
which is responsible for loading the dataset, covariance matrix and
applying cuts.

Specifying a dataset
~~~~~~~~~~~~~~~~~~~~

In a validphys runcard the settings for a single dataset are specified
using a ``dataset_input``. This is a dictionary which minimally
specifies the name of the dataset, but can also control behaviour such
as contributions to the covariance matrix for the dataset and
C-factors.

Here is an example dataset input:

.. code:: yaml

    dataset_input:
        dataset: CMSZDIFF12
        cfac: [QCD,NRM]
        sys: 10

This particular example is for ``CMSZDIFF12`` dataset, the user has
specified to use some C-factors ``cfac`` and ``sys: 10`` which corresponds
to an additonal contribution to the covariance matrix accounting for
statistical fluctuations in the C-factors. These settings correspond to
NNLO predictions and so presumably elsewhere in the runcard the user
would have specified a NNLO theory - such as theory 53.

We can use the API to return an instance of ``DataSetSpec`` in a
development environment using the settings above

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
:py:mod:`validphys.config.CoreConfig.produce_dataset`, the required arguments are
``dataset_input``, ``cuts`` and ``theoryid``.

.. note::
    It seems odd to require a `theoryid`
    and parameters in the `dataset_input` which refer to theory settings
    in order to load data. However, this is a relic of the underlying c++ code
    which performs the loading of data, which intrinsically groups together the
    commondata (CSVs containing data central values and uncertainties) and :ref:`fktables`.

    Clearly there is a big margin for error when manually entering `dataset_input`
    and so there is a [project](https://github.com/NNPDF/nnpdf/issues/226) which
    aims to have a stable way of filling many of these settings with correct
    default values.

The ``DataSetSpec`` contains all of the information used to construct
it, e.g

.. code:: python

    >>> ds_spec.thspec
    TheoryIDSpec(id=53, path=PosixPath('/Users/michael/conda/envs/nnpdf-dev/share/NNPDF/data/theory_53'))
    >>> ds_spec.name
    'CMSZDIFF12'

but also importantly has a ``load`` method, which returns an instance of
the ``DataSet``, which was generated from the c++ code using SWIG. This
new object contains numpy arrays of data central values and experimental
covariance matrices, e.g:

.. code:: python

    >>> ds_libnnpdf.get_cv() # get central values of dataset
    array([2917.  , 1074.  ,  460.5 ,  222.6 ,  109.8 ,   61.84,   30.19,
           2863.  , 1047.  ,  446.1 ,  214.5 ,  110.  ,   58.13,   29.85,
           2588.  ,  935.5 ,  416.3 ,  199.  ,  103.1 ,   54.06,   28.45,
           1933.  ,  719.5 ,  320.7 ,  161.1 ,   84.62,   47.57,   24.13])

In practice actions which require experimental data and/or covariance
matrices will make use of the :py:mod:`validphys.results.results`
provider which is a tuple of :py:mod:`validphys.results.DataResult`
and :py:mod:`validphys.results.ThPredictionsResult`. Since we are additionally
generating theory predictions we additionally are required to specify a
PDF

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

The covariance matrix associated with the ``DataResult`` in this tuple
was constructed by :py:mod:`validphys.results.covmat` which allows the
user to change the behaviour of the covariance matrix - such as adding
theory uncertainties from scale variation or using a t0 pdf to calculate
the multiplicative contributions to the covariance matrix - for more
detail see :py:mod:`validphys.results.covmat`.

``DataGroupSpec`` - core object for multiple datasets
-----------------------------------------------------

The core object for multiple datasets is :py:mod:`validphys.core.DataGroupSpec`
which is similar in many regards to the DataSetSpec, but handles the loading
of multiple datasets. In particular, when constructing the covariance matrix,
it takes into account any uncertainties which are correlated across different
datasets.

Specifying multiple datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Multiple datasets are specified using ``dataset_inputs`` key: a list
where each element of the list is a valid ``dataset_input``. For
example:

.. code:: yaml

    dataset_inputs:
        - { dataset: NMC }
        - { dataset: ATLASTTBARTOT, cfac: [QCD] }
        - { dataset: CMSZDIFF12, cfac: [QCD,NRM], sys: 10 }

We see that multiple datasets are inputted as a flat list and there is
no hierarchy to the datasets, splitting them into experiments or process
types. The grouping of datasets is done internally according to the
metadata of datasets and is controlled by ``metadata_group`` key. This
can be any key which is present in the ``PLOTTING`` file of each dataset
- for example ``experiment`` or ``nnpdf31_process``.
The default value for ``metadata_group`` is ``experiment``. Other groupings
might be relevant for example when contructing a theory covariance matrix, where
you want to group datasets according to process type rather than experiment.
The grouping is performed by the production rule
:py:mod:`validphys.config.CoreConfig.produce_group_dataset_inputs_by_metadata`
which returns a list with length equal to number of distinct groups. Each element
is a namespace with the ``group_name`` and list of ``dataset_input`` s for that
specific group e.g:

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

Here we see that the namespace key is ``data_input`` rather than ``dataset_inputs``
which is because ``data_input`` bridges the gap between the current way of specifying
data (with ``dataset_inputs``) and a deprecated specification using the ``experiments``
key. The production rule which returns a ``DataGroupSpec`` is
:py:mod:`validphys.config.CoreConfig.produce_data` through the following
pipeline

.. code::

    dataset_inputs or experiments -> data_input -> data

For example the following runcard produces a single column table with a
row containing the 𝞆² of the specificed datasets, grouped by
``experiment``

.. code:: yaml

    dataset_inputs:
        - { dataset: NMC }
        - { dataset: ATLASTTBARTOT, cfac: [QCD] }
        - { dataset: CMSZDIFF12, cfac: [QCD,NRM], sys: 10 }

    theoryid: 53

    dataspecs:
     - pdf: NNPDF31_nnlo_as_0118

    use_cuts: internal

    actions_:
     - dataspecs_groups_chi2_table

If we specify a grouping in the runcard:

.. code:: yaml

    metadata_group: nnpdf31_process

    dataset_inputs:
        - { dataset: NMC }
        - { dataset: ATLASTTBARTOT, cfac: [QCD] }
        - { dataset: CMSZDIFF12, cfac: [QCD,NRM], sys: 10 }

    theoryid: 53

    dataspecs:
     - pdf: NNPDF31_nnlo_as_0118

    use_cuts: internal

    actions_:
     - dataspecs_groups_chi2_table

then we instead get a single column table, but with the datasets grouped
by process type, according the `theory uncertainties
paper <https://arxiv.org/abs/1906.10698>`__.

To expose the result of checking if ``metadata_group`` was specified in the
runcard and applying a default value, use the namespace key
``processed_metadata_group``. We can use this key in reports and actions alike
for example to give sensible titles/section headings e.g:

`` code:: yaml
    template_text: |
     # chi2 grouped by {processed_metadata_group}
     {@dataspecs_groups_chi2_table@}

    actions_:
     - report(main=True)

Action naming conventions
-------------------------

There are some general rules which should be observed when adding
new actions to ``validphys``. Firstly try to indicate the required runcard input
for an action in the name of the function. Take for example the provider
``dataset_inputs_results``. The returned object is a ``results`` object: a tuple
of data and theory which is used by a wide range of other actions, notably
when calculating 𝞆². The first part of the name ``dataset_inputs`` refers
to the runcard input required to process that action. This is especially
useful with actions for a group of datasets or ``data``, because the dependency
tree for these actions is not neccessarily obvious to somebody who is unfamiliar
with the code. As explained above,
``dataset_inputs -> data_input -> data`` and so the action name serves to
guide the user to creating a working runcard as easily as possible.

The second general rule is if your action has makes use of ``collect`` somewhere
in the dependency graph, then consider prepending what is collected over to
the action name. For example: ``dataspecs_groups_chi2_table``, which depends on

.. code:: python
    dataspecs_groups_chi2_data = collect("groups_chi2", ("dataspecs",))

and in turn

.. code:: python
    groups_chi2 = collect("dataset_inputs_abs_chi2_data", ("group_dataset_inputs_by_metadata",))

Without having to find these specific lines in the code we were able to guess
that 𝞆² was collected first over groups of data (``groups_chi2``), and then
over ``dataspecs``. Naming functions according to these rules helps make the
general workings of the underlying code more transparent to an end user.

Backwards compatibility
-----------------------

Where possible, backwards compatibility with runcards which use the ``experiments``
key has been preserved. For example with the ``dataspecs_groups_chi2_table``
example above we could also use the following input

.. code:: yaml

    metadata_group: nnpdf31_process

    dataset_inputs:
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

    use_cuts: internal

    actions_:
     - dataspecs_groups_chi2_table

The user should be aware, however, that any grouping introduced in this way
is purely superficial and will be ignored in favour of the experiments defined
by the metadata of the datasets.

Note that some theory uncertainties runcards will need to be updated to explicitly
set the dataset grouping to ``experiment``.

Runcards which request actions that have been renamed won't work anymore,
generally actions which were previously named ``experiments_*`` have been
renamed to highlight that they work with more general groupings.

Currently ``n3fit``, and the ``pseudodata``, ``closuretest`` and ``chi2grids``
modules have not been updated to use ``dataset_inputs`` yet and so require
``experiments`` to be specified in the runcard. The c++ fitting code ``nnfit``
is not scheduled to be updated to use ``dataset_inputs`` and so will always
require ``experiments`` to be specified in the runcard.
