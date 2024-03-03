.. _commondata:

=======================
Experimental data files
=======================

Data made available by experimental collaborations comes in a variety of
formats. For use in a fitting code, this data must be converted into a common
format that contains all the required information for use in PDF fitting.
Existing formats commonly used by the community, such as in `HepData <https://www.hepdata.net/>`_,
are generally unsuitable. Principally as they often do not fully describe the
breakdown of systematic uncertainties. Therefore over several years an NNPDF
standard data format has been iteratively developed, now denoted ``CommonData``.

This documentation describes the ``CommonData`` format
used in NNPDF starting from code version 4.0.10 and compatible with releases beyond 4.0.


Naming convention and organization of the datasets
--------------------------------------------------

All datasets in the new data format follow the exact same naming convention::

    <setname>_<observable>

where the setname is defined by::

    <experiment>_<process>_<energy>{_<extras>}

The naming convention for the set names is defined in the :ref:`naming convention documentation<dataset-naming-convention>`.

Each ``<setname>`` defines a folder in which the data is contained.
While the separation of data in different folders can be arbitrary,
a folder cannot contain more than one hepdata entry
or datasets that mix different processes, energies or experiment.
Due to historical reasons and for backwards compatibility the special energy ``NOTFIXED`` is used
for datasets where more than one center of mass energy is used.
When in doubt, it is preferable to utilize two different folders.
The ``<extras>`` string is free and can be used to disambiguate.

The data downloaded or parsed from hepdata or other sources is kept in the
``<setname>/<rawdata>`` folder and it is not installed with the rest of the code.
Each folder must contain a ``<setname>/metadata.yaml`` file which will define
all datasets implemented within the folder and that will be described below.
Only ``.yaml`` file are allowed to be installed together with the ``nnpdf`` code.

In order to keep backward compatibility and allow the reproducibility of the 4.0 family of fits
a ``dataset_names.yml`` file keeps a mapping of the datasets that were used in 4.0.
When using the old names in a runcard, ``validphys`` will automatically translate
them using this file.
The format of this mapping is as follow:

..  code-block:: yaml

    old_name_1:
        dataset: new_name_1
        variant: legacy

                                 
CommonData Metadata specification
---------------------------------

The ``metadata.yaml`` file defines unequivocally the datasets implemented within a folder.
The general structure is a first portion of general information (references, name of the set)
and a list of ``implemented_observables`` which define separate datasets.


Shared information
==================


..  code-block:: yaml

    setname: "EXPERIMENT_PROCESS_ENERGY{_EXTRA}"

    version: 1
    version_comment: "A comment about this version"

    # References
    arXiv:
        url: "https://arxiv.org/abs/XYZ.ABC"
    iNSPIRE:
        url: "https://inspirehep.net/literature/XYZ"
    hepdata:
        url: "https://www.hepdata.net/record/insXYZ"
        version: 1

    nnpdf_metadata:
        nnpdf31_process: "PROCESS"
        experiment: "EXPERIMENT_NAME"

    implemented_observables:
      - observable_metadata_1
      - observable_metadata_2


The header of the ``metadata.yaml`` file contains information shared among different datasets.

Setname
~~~~~~~

Correspond to the name of the set and must be equal to the folder. It acts a s a sanity check.

Versioning
~~~~~~~~~~

The initial version of a dataset should be set to ``version: 1``.
Any change on a dataset should be *always* accompanied of a version bump and a ``version_comment`` explaining the update.
This will allow to keep an exact tracking of all changes to every dataset even if they change over time due to bugs, updates in hepdata, etc.

References
~~~~~~~~~~

References to the original source of the data. 
This can be ``arXiv``, ``iNSPIRE`` or ``hepdata``.
All information must be provided unless it is explicitly missing.

nnpdf_metadata
~~~~~~~~~~~~~~

Grouping information used internally by ``validphys`` up to NNPDF4.0.
It accepts the keys ``experiment``, which should in general coincide
with the ``EXPERIMENT`` key in the ``<setname>`` and the key ``nnpdf31_process``
which is the process grouping information used in the 3.1 and 4.0 MHOU papers.

Observable specific information
===============================

Within a ``metadata.yaml`` we can find one or more implemented datasets.
These correspond to different observables of a single measurement.
For instance, the LHCB publication of Z rapidity measurements at 13 TeV
(``setname: LHCB_Z0_13TEV``) contains two observables: Z decay into two electrons
and Z decay into 2 muons.
This setname contain two datasets: ``LHCB_Z0_13TEV_DIELECTRON-Y`` and ``LHCB_Z0_13TEV_DIMUON-Y``.

In the following we describe the metadata corresponding to the observable within the ``metadata.yaml`` file.


..  code-block:: yaml
    
   implemented_observables:
    - observable_name: "DIMUON-Y"
      process_type: "EWK_RAP"
      tables: [5]
      ndata: 18
      observable:
        description: "Differential cross-section of Z-->µµ as a function of Z-rapidity"
        label: r"$d\sigma / d|y|$"
        units: "[fb]"
      kinematics:
        file: kinematics_dimuon.yaml
        variables:
          y: {description: "Z boson rapidity", label: "$y$", units: ""}
          M2: {description: "Z boson Mass", label: "$M^2$", units: "$GeV^2$"}
          sqrts: {description: "Center of Mass Energy", label: '$\sqrt{s}$', units: "$GeV$"}
      kinematic_coverage: [y, M2, sqrts]
      data_central: data_dimuon.yaml
      data_uncertainties:
        - uncertainties_dimuon.yaml
      variants:
        - example_variant:
            data_uncertainties:
              - uncertainties_different_treatment.yaml
      theory:
        FK_tables:
          - - LHCB_DY_13TEV_DIMUON
        operation: 'null'
        conversion_factor: 1000.0
      # Plotting information
      plotting:
        dataset_label: "LHCb $Z\\to µµ$"
        plot_x: y
        y_label: '$d\sigma_{Z}/dy$ (fb)'

``observable_name``
~~~~~~~~~~~~~~~~~~~
The observable name is used to construct the full name of the dataset ``<setname>_<observable_name>``.
It must be unique within a set and contain no ``_`` (as it could lead to confusion).

``process_name``
~~~~~~~~~~~~~~~~
One of the processes defined in the ``process_options`` module at
``validphys/src/validphys2/process_options.py``.
This is used internally by validphys to describe the combination of observable
and process in various plots, to check that the kinematic variables utilized by the
dataset are sensible and to generate derived plots such as the ``x-q2`` kinematic coverage plots.

``tables``
~~~~~~~~~~
Tables from the hepdata entries that have been used to construct the dataset

``ndata``
~~~~~~~~~
Number of datapoints in the dataset.
While this quantity could be derived from the data itself,
many other pieces (crucially backwards compatibility with cuts and theories) requires
the number of datapoints to be set in stone.
If an update requires to change the number of datapoint,
it should be added as a separate observable.

``observable``
~~~~~~~~~~~~~~
This is a dictionary with the entries ``description``, ``label`` and ``units``.
All entries must be latex-compilable as they are used by various plotting routines in ``validphys``.

``kinematics::file``
~~~~~~~~~~~~~~~~~~~~
A reference to a ``.yaml`` file containing all kinematic information.
The file contain a list of ``ndata`` ``bins`` for which information about all variables
is included for all bins.
When ``mid`` is not given, it will be automatically filled with the midpoint between min and max.
Only ``mid`` is used for cuts, while ``min`` and ``max`` may be used for plotting routines.

..  code-block:: yaml

    bins:
        - var_1:
            min: 0
            max: 1
            mid: 0.5
          var_2:
            min: 0
            max: 1
            mid: 0.5

``kinematics::variables``
~~~~~~~~~~~~~~~~~~~~~~~~~
Metadata for each of the variables contained in the ``kinematics::file``
and which can be ``description``, ``label`` and ``units``.
Latex syntax is accepted and encouraged since they will be used by plotting routines.

..  code-block:: yaml

    variables:
      var_1: {description: "my var 1", label: "$m$", "units: "GeV"}


``kinematic_coverage``
~~~~~~~~~~~~~~~~~~~~~~
A list of the variables within the kinematic files


``data_central``
~~~~~~~~~~~~~~~~
A reference to a ``yaml`` file containing the measurement central data.
The format of the data is a ``yaml`` file with an entry ``data_central`` which
list for all values for all bins.

..  code-block:: yaml

    data_central:
        - val1
        - val2
        - val3

``data_uncertainties``
~~~~~~~~~~~~~~~~~~~~~~
A list of ``.yaml`` file containing the uncertainty information for the measurement.
When using more than one uncertainty file they will be concatenated. 
This allows the user the flexibility of creating variants
where only a subset of the uncertainties are modified.

The format of the uncertainty files is of two fields, a ``definitions`` field that contains
metadata about all the uncertainties: name, treatment (``ADD`` or ``MULT``) and type
and a second field ``bins`` which is a list of mappings with ``ndata`` entries
with the named uncertainties.

Note that, regardless of their treatment, uncertainties should always be written as absolute values
and not relative to the data values. If the data should be updated, the uncertainties should be too.

..  code-block:: yaml

    definitions:
        stat:
            description:
            treatment:
            type:
        error_name:
            description:
            treatment:
            type:
        error_name_2:
            description:
            treatment:
            type:
    bins:
        - stat:
          error_name:
          error_name_2:




``variants``
~~~~~~~~~~~~

In some occasions we might want to maintain two variations of the same observable.
For instance, we might have two incompatible sources of uncertainties. In such case a variant can be added.
These variants can overwrite certain keys if necessary.
When a variant is used, the key under the variant will be used instead of the key defined in the observable.

A ``variant`` can only overwrite the entries ``data_central``, ``theory`` and ``data_uncertainties``.
Example:

..  code-block:: yaml

    data_uncertainties:
        - uncertainties.yaml

    variants:
        name_of_the_variant:
            data_uncertainties:
                - uncertainties.yaml
                - extra_uncertainties.yaml
        another_variant:
            data_central: different_data.yaml
            data_uncertainties:
                - different_uncertainties.yaml
              
When loading this dataset with no variant only the ``uncertainties.yaml`` file will be read.
Instead, when choosing ``variant: name_of_the_variant``, both ``uncertainties.yaml`` and  ``extra_uncertainties.yaml`` will be loaded.
If we select ``variant: another_variant`` both the ``data_uncertainties`` and the ``data_central`` keys will be substituted.
Note that if we want to substitute the default set of uncertainties we just need to not include it in the variant (as done in ``another_variant``).

``theory``
~~~~~~~~~~

The theory field defines how predictions for the dataset are to be computed.
It includes two entries:

- ``FK_tables``: this is a list of lists which defines the FK Tables to be loaded. The outermost list are the operands (in case an operation is needed to recover the observable, more on that below). The innermost list are the grids that are to be concatenated in order to form the operands.
- ``operaton``: operation to be applied in order to compute the observable

Example:

..  code-block:: yaml
  
  theory: 
  FK_tables:
      - - Z_contribution
        - Wp_contribution
        - Wm_total
      - - total_xs
  operation: 'ratio'

In this case the ``fktables`` for the Z, W+ and W- contributions will be concatenated (the dataset might include predictions for all three contributions).
After that, the final observable will be computed by taking the ratio of the concatenation of all those observables and the total cross section (``total_xs``).

``plotting``
~~~~~~~~~~~~

The ``plotting`` section defines the plotting style inside ``validphys``
and is described in detail in :ref:`plotting-format`.

Note that name of the variables need to be the same in the plotting and kinematics.
