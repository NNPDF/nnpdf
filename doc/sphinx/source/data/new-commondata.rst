Naming convention and organization of the datasets
--------------------------------------------------

All datasets in the new data format follow the exact same naming convention::

    <experiment>_<process>_<energy>{_<extras>}_<observable>

The data is contained in folders, each folder containing one single hepdata publication. 
In all cases one can reconstruct the name of the folder by separating the observable name on the last ``_``, i.e., the folder will always be named::

    <experiment>_<process>_<energy>{_<extras>}

Where all observables contained in one hepdata entry are separated by their observable name.

Each folder will contain one single metadata file named ``metadata.yaml`` which defines all observables implemented for a given dataset.

In order to keep backward compatibility and ease the comparison between new and old commondata, the ``buildmaster/dataset_names.yml`` file keeps a mapping of the datasets implemented in both formats.
When a ``legacy`` variant is available, the usage of the old name automatically enables such variants. The format of this mapping is as follow (which enables using variants):

..  code-block:: yaml

    old_name_1: new_name_1
    old_name_2:
        dataset: new_name_2
        variant: this_particular_variant


Metadata Format
---------------

This ``metadata.yaml`` file contains a first portion of general information which might be shared by several sets and a list of ``implemented_observables`` which define the separate observables.


..  code-block:: yaml

    setname: "EXPERIMENT_PROCESS_ENERGY{_EXTRA}"

    version: 1
    version_comment: "Initial implementation"

    # References
    arXiv:
        url: ""
    iNSPIRE:
        url: "https://inspirehep.net/literature/302822"
    hepdata:
        url: "https://www.hepdata.net/record/ins302822"
        version: 1

    nnpdf_metadata:
        nnpdf31_process: "PROCESS"
        experiment: "EXPERIMENT_NAME"

    implemented_observables:
      - observable_name: "OBS"
        observable:
            description: "Description of the observable"
            label: "Latex label for the observable"
            units: "[u]"
        ndata: n_of_datapoints
        tables: [n, j, k] # (optional) corresponding tables in the hepdata entry
        npoints: [n, j, k] # (optional) number of points per table
        process_type: INC # for instance, INC, JET, DIJET, etc

        # Plotting information (for instance, the kinematics variable could be pt, mt, q2)
        plotting:
            dataset_label: "Label to be used in reports"
            kinematics_override: identity
            x_scale: log
            plot_x: var_1
            figure_by:
                - var_2

        kinematic_coverage: [var_1, var_2, var_3]

        kinematics:
            variables:
                var_1: {description: "Description of var", label: "latex", units: "u"}
                var_2: {description: "Description of var", label: "latex", units: "u"}
                var_3: {description: "Description of var", label: "latex", units: "u"}
            file: kinematics.yaml

        data_central: data.yaml
        data_uncertainties:
            - uncertainties.yaml
            - uncertainties_2.yaml

        # Having variants is optional
        # variants can overwrite the data_uncertainties 
        variants:
            different_errors:
                data_uncertainties:
                    - uncertainties.yaml
                    - uncertainties_3.yaml

        # The theory field is always optional
        theory: 
            FK_tables:
                - - DYE605
            operation: 'null'




Versioning
~~~~~~~~~~

The initial version of a dataset should be set to ``version: 1``.
Any change on a dataset should be *always* accompanied of a version bump and a ``version_comment`` explaining the update.
This will allow to keep an exact tracking of all changes to every dataset even if they change over time.

Variants
~~~~~~~~

In some occasions we might want to maintain two variations of the same observable.
For instance, we might have two incompatible sources of uncertainties. In such case a variant can be added.
The syntax of the ``variants`` is.

Theory
~~~~~~

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


..  code-block:: yaml

    data_uncertainties:
        - uncertainties.yaml

    variants:
        name_of_the_variant:
            data_uncertainties:
                - uncertainties.yaml
                - extra_uncertainties.yaml
        another_variant:
            data_uncertainties:
                - different_uncertainties.yaml


When loading this dataset with no variant only the ``uncertainties.yaml`` file will be read.
Instead, when choosing ``variant: name_of_the_variant``, both ``uncertainties.yaml`` and  ``extra_uncertainties.yaml`` will be loaded.
Note that if we want to substitute the default set of uncertainties we just need to not include it in the variant (as done in ``another_variant``).


Data
----

The format of the data is a ``yaml`` file with an entry ```data_central``` which is a list for all values for all bins.

..  code-block:: yaml

    data_central:
        - val1
        - val2
        - val3

Uncertainties
-------------

The uncertainties are (also) ``.yaml`` files. 
Note that in the ``metadata.yaml`` the ``data_uncertainties`` entry is given as a list. 
When using more than one uncertainty file they will be concatenated. 
This allows the user the flexibility of creating variants where only a subset of the uncertainties are modified.

The format of the uncertainty files is of two fields, a ``definitions`` field that contains metadata about all the uncertainties (their name, their treatment (``ADD`` or ``MULT``) and their type) and a second field ``bins`` which is a list of mappings with as many entries as the `data_central` with the named uncertainties.

Note that, regardless of their treatment type, the uncertainties should always be written as absolute values and not relative to the data values.

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

Kinematics:
-----------
The kinematics file follow a convention very similar to the uncertainties file, where the ``definitions`` field is skipped since that information is already contained in the parent ``metadata.yaml`` file.

Therefore, we have a list of ``bins`` (of the same size as the list for `data_central`) and for each entry we have the information of all the variables.

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

Plotting
~~~~~~~~

The ``plotting`` section defines the plotting style inside ``validphys``.
In previous implementations there were per-process options that defined plotting options for family of processes.
In the commondata format defined in this page every plotting option must be defined in the ``plotting`` section of each observable.

Internally within ``validphys`` only 3 kinematic variables are taken into account. The 3 selected variables (and their order) is defined by ``plotting::kinematic_coverage``.

The name of the variables (which in this example are `var_1`, `var_2`, `var_3`) need to be the same in the plotting and kinematics.
