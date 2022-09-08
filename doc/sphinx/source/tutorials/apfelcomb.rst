.. _tutorialfktables:

How to generate and implement FK tables
=======================================

APFELcomb is the project that allows the user to generate ``FK`` tables.
These are lookup tables that contain the relevant information to compute
theoretical predicitons in the NNPDF framework. Broadly speaking, this
is achieved by taking DGLAP evolution kernels from ``APFEL`` and
combining them with interpolated parton-level observable kernels in the
APPLgrid or FastNLO format (see `How to generate APPLgrid and fastNLO
tables <../tutorials/APPLgrids>`__). The various data formats used in
APFELcomb are described in `Experimental data
files <../data/exp-data-files.rst#exp-data-files>`__.

The user is strongly encouraged to go through that note with care, in
order to familiarise himself with the features and the structure of the
APFELcomb project.

Generate a FK table
-------------------

Each ``FK`` table is generated piecewise in one or more ``subgrids``.
The ``subgrids`` implemented in APFELcomb can be displayed by running
the script

::

   ./scripts/disp_grids.py

The generation of each subgrid can by achieved with the following
command

::

   ./apfel_comb <source=app/dis/dyp> <subgrid id> <theory id>

where ``<app/dis/dyp>`` specifies whether the subgrid is in the APP, DIS
or DYP subgrid categories in the database (``db/apfelcomb.dat``), where:
- APP: refers to applgrids, partonic cross sections produced externally
by a MonteCarlo generator. - DIS: Deep Inelastic Scatting, coefficient
fucnctions computed by ``APFEL``. - DYP: Drell-Yan, partonic cross
sections computed by ``APFEL``.

``<subgrid id>`` is the corresponding ID in that database (visible in
the ``disp\_grids`` script) and ``<theory id>`` specifies the desired
NNPDF theory index (the entry in nnpdf/nnpdfcpp/data/theory.db). As an
example:

.. code:: shell

   ./apfel_comb app 500 53

will generate the subgrid for CDFZRAP and theory 53 (NNPDF3.1 NNLO
fitted charm). The resulting FK subgrid will be written out to

::

   $RESULTS_PATH/theory_<theory id>/subgrids/FK_<setname>_<subgrid id>.dat.

APPLgrids and FastNLO tables should be properly stored in the
``applgrids`` folder by means of `Git
LFS <https://git-lfs.github.com/>`__ (see `here <storage>`__ for
details).

Once all the relevant subgrids for the desired dataset(s) are generated,
one should run

::

   ./merge_allgrids.py <theory id>

which will loop over all datasets and attempt to merge their subgrids
into a complete ``FK`` table. The resulting final ``FK`` table should be
stored at

::

   $RESULTS_PATH/theory_<theory id>/fastkernel/FK_<setname>.dat.

Implement a new FK table
------------------------

Whenever a new dataset is implemented, it should be accompanied by the
corresponding ``FK`` table. To implement a new ``FK`` table, one must
first add a corresponding entry into the apfelcomb database (by editing
the ``./db/apfelcomb.dat`` file) under the ``grids`` table. These
entries are comprised of the following fields. - **id** - The primary
key identifier of the FK table. - **setname** - The COMMONDATA set name
of the corresponding dataset. - **name** - The name of the FK table. -
**description** - A one-line description of the FK table. - **nx** - The
number of x-grid interpolation points. - **positivity** - A flag
specifying if the FK table is a positivity set. - **source** - Specifies
if the corresponding subgrids are [APP/DIS/DYP].

Note that **setname** and **name** may be different in the case of
compound observables such as ratios, where multiple FK tables are
required to compute predictions for a single dataset. The ``nx``
parameter specifies the interpolation accuracy of the dataset (this must
currently be tuned by hand, e.g. by making sure that the native applgrid
and the generated FK tables lead to numerically equivalent results once
they are convolved with the same PDF set). The ``positivity`` parameter
restricts the observable to NLO matrix elements and disables target-mass
corrections. Once this entry is complete, one must move on to adding
entries in the corresponding subgrid table.

Implementing a new APPLgrid/FastNLO subgrid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To add a new APPLgrid- or FastNLO–based subgrid, one must add a
corresponding entry into the ``app\_subgrids`` table of the apfelcomb
database. One entry should be added for each APPLgrid making up the
final target ``FK`` table. The entries have the following fields: -
**id** - The primary key identifier of the subgrid. - **fktarget** - The
name of the FK table this subgrid belongs to. - **applgrid** - The
filename of the corresponding APPLgrid. - **fnlobin** - The fastNLO
index if the table is a fastNLO grid, or -1 if not. - **ptmin** - The
minimum perturbative order (1 when the LO is zero, 0 if not). -
**pdfwgt** - A boolean flag, 1 if the APPLgrid has PDF weighting, 0 if
not (depending on how the native applgrid was generated). - **ppbar** -
A boolean flag, 1 if the APPLgrid should be transformed to *ppbar*
beams, 0 if not. - **mask** - A boolean mask, specifying which APPLgrid
entries should be considered data points. - **operators** - A list of
operators to handle certain special cases (see below). The mask should
have as many entries as APPLgrid bins and each boolean value should be
separated by a space. For example, for an applgrid with five bins where
we want to exclude the penultimate bin, the mask would be:

::

   1 1 1 0 1

Note that there is no way to know a priori whether ``pdfwgt`` should be
set to 0 or to 1, that is whether the grid is unweighted or weighted.
However, this can easily be checked a posteriori, since setting
``pdfwgt`` to the wrong value should lead to ``./apfel_comb`` failing
due to a large relative error between the value in the APPLgrid and that
in the FK table.

The applgrid filename assumes that the grid can be found at

::

   $APPL_PATH/<setname>/<applgrid>

where ``APPL_PATH`` is defined in Makefile.am, ``<setname>`` is the
corresponding ``COMMONDATA`` set name specified in the grids table (that
should match the name used in the
`buildmaster <../tutorials/buildmaster.md>`__ implementation), and
``<applgrid>`` is specified in the field described above.

Implementing a new DIS or DYP subgrid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

New DIS or DYP subgrids should be entered respectively into the
``dis_subgrids`` or ``dyp_subgrids`` tables of the apfelcomb database.
Typically only one subgrid is needed per DIS or DYP FK table. Each
subgrid entry has the following fields: - **id** - The primary key
identifier of the subgrid - **fktarget** - The name of the FK table this
subgrid belongs to - **operators** - A list of operators to handle
certain special cases (see Subgrid operators). For DIS there is one
additional field: - **process** - The process string of the observable
(e.g DIS_F2P, see DIS Processes in APFEL below)

DIS Processes in APFEL
~~~~~~~~~~~~~~~~~~~~~~

For DIS processes and since the coefficient functions are computed
solely with APFEL, one needs to specify the process of the observable,
in ``dis_subgrids`` following ``APFEL``\ ’s nomenclature. The list of
processes below can be found in ``apfel/src/DIS/FKObservables.f`` in the
headers corresponding to the different observables called.

**Deep Inelastic Scattering Structure Functions**: - DIS_F2L: [EM] Light
structure function F2light (electron-proton) - DIS_F2U: [EM] Up
structure function F2u (electron-proton[up]) - DIS_F2d: [EM] Down
structure function F2d (electron-proton[down]) - DIS_F2S: [EM] Strange
structure function F2s (electron-proton[strange]) - DIS_F2C: [EM] Charm
structure function F2charm (electron-proton) - DIS_F2B: [EM] Bottom
structure function F2bottom (electron-proton) - DIS_F2T: [EM] Top
structure function F2top (electron-proton) - DIS_F2D: [EM] Deuteron
structure function F2 (electron-isoscalar) - DIS_FLL: [EM] Light
structure function FLlight (electron-proton) - DIS_FLC: [EM] Charm
structure function FLcharm (electron-proton) - DIS_FLB: [EM] Bottom
structure function FLbottom (electron-proton) - DIS_FLT: [EM] Top
structure function FLtop (electron-proton) - DIS_FLD: [EM] Deuteron
structure function FL (electron-isoscalar) - DIS_F2P_NC: [NC] Proton
structure function F2 (electron-isoscalar) - DIS_F2P: [EM] Proton
structure function F2 (electron-proton) - DIS_FLP_NC: [NC] Proton
structure function FL (electron-proton) - DIS_FLP_CON_NC: [NC] Proton
structure function FL (electron-proton) - DIS_FLP: [EM] Proton structure
function FL (electron-proton) - DIS_F3P_NC: [NC] F3 structure function
(electron-proton)

**Deep Inelastic Scattering Reduced Cross-Sections**: - DIS_NCE_L: [NC]
Electron scattering Reduced Cross-Section, light (electron-proton) -
DIS_NCP_L: [NC] Positron scattering Reduced Cross-Section, light
(positron-proton) - DIS_NCE_CH: [NC] Electron scattering Reduced
Cross-Section, charm (electron-proton) - DIS_NCP_CH: [NC] Positron
scattering Reduced Cross-Section, charm (positron-proton) - DIS_NCE_BT:
[NC] Electron scattering Reduced Cross-Section, bottom (electron-proton)
- DIS_NCP_BT: [NC] Positron scattering Reduced Cross-Section, bottom
(positron-proton) - DIS_NCE_TP: [NC] Electron scattering Reduced
Cross-Section, top (electron-proton) - DIS_NCP_TP: [NC] Positron
scattering Reduced Cross-Section, top (positron-proton) - DIS_NCE_D:
[NC] Electron scattering Reduced Cross-Section on deuteron, inclusive
(electron-isosclar) - DIS_NCP_D: [NC] Positron scattering Reduced
Cross-Section on deuteron, inclusive (positron-isoscalar) - DIS_NCE:
[NC] Electron scattering Reduced Cross-Section, inclusive
(electron-proton) - DIS_NCP: [NC] Positron scattering Reduced
Cross-Section, inclusive (positron-proton) - DIS_CCE_L: [CC] Electron
scattering Reduced Cross-Section, light (electron-proton) - DIS_CCP_L:
[CC] Positron scattering Reduced Cross-Section, light (positron-proton)
- DIS_CCE_C: [CC] Electron scattering Reduced Cross-Section, charm
(electron-proton) - DIS_CCP_C: [CC] Positron scattering Reduced
Cross-Section, charm (positron-proton) - DIS_CCE: [CC] Electron
scattering Reduced Cross-Section, inclusive (electron-proton) - DIS_CCP:
[CC] Positron scattering Reduced Cross-Section, inclusive
(positron-proton)

**Deep Inelastic Scattering Reduced Cross-Sections (heavy-ion)**: -
DIS_SNU_L_Pb: [CC] Neutrino scattering Reduced Cross-Section, light
(neutrino-lead) - DIS_SNB_L_Pb: [CC] Antineutrino scattering Reduced
Cross-Section, light (antineutrino-lead) - DIS_SNU_C_Pb: [CC] Neutrino
scattering Reduced Cross-Section, charm (neutrino-lead) - DIS_SNB_C_Pb:
[CC] Antineutrino scattering Reduced Cross-Section, charm
(antineutrino-lead) - DIS_SNU_Pb: [CC] Neutrino scattering Reduced
Cross-Section, inclusive (neutrino-lead) - DIS_SNB_Pb: [CC] Antineutrino
scattering Reduced Cross-Section, inclusive (antineutrino-lead) -
DIS_SNU_L: [CC] Neutrino scattering Reduced Cross-Section, light
(neutrino-isoscalar) - DIS_SNB_L: [CC] Antineutrino scattering Reduced
Cross-Section, light (antineutrino-isoscalar) - DIS_SNU_C: [CC] Neutrino
scattering Reduced Cross-Section, charm (neutrino-isoscalar) -
DIS_SNB_C: [CC] Antineutrino scattering Reduced Cross-Section, charm
(antineutrino-isoscalar) - DIS_SNU: [CC] Neutrino scattering Reduced
Cross-Section, inclusive (neutrino-isoscalar) - DIS_SNB: [CC]
Antineutrino scattering Reduced Cross-Section, inclusive
(antineutrino-isoscalar) - DIS_DM_NU: [CC] Dimuon neutrino cross section
(neutrino-iron) - DIS_DM_NB: [CC] Dimuon anti-neutrino cross section
(antineutrino-iron)

**Single-Inclusive electron-positron annihilation, Time-Like Evolution
(SIA)**: - SIA_F2: [NC] SIA structure function F2 = FT + FL
(electron-proton) - SIA_FL: [NC] SIA structure function FL
(electron-proton) - SIA_FA: [NC] SIA structure function FA
(electron-proton) - SIA_XSEC_NF4: [NC] SIA absolute cross section (nf=4)
(electron-proton) - SIA_XSEC: [NC] SIA absolute cross section
(electron-proton) - SIA_NORM_XSEC_LONG_L: [NC] SIA normalized light
longitudinal cross section (electron-proton) - SIA_NORM_XSEC_LONG_BT:
[NC] SIA normalized bottom longitudinal cross section (electron-proton)
- SIA_NORM_XSEC_LONG: [NC] SIA normalized total longitudinal cross
section (electron-proton) - SIA_NORM_XSEC_L: [NC] SIA normalized light
cross section (electron-proton) - SIA_NORM_XSEC_CH: [NC] SIA normalized
charm cross section (electron-proton) - SIA_NORM_XSEC_BT: [NC] SIA
normalized bottom cross section (electron-proton) - SIA_NORM_XSEC_TP:
[NC] SIA normalized top cross section (electron-proton) -
SIA_NORM_XSEC_NF4: [NC] SIA normalized total cross section (nf=4)
(electron-proton) - SIA_NORM_XSEC: [NC] SIA normalized total cross
section (electron-proton)

Subgrid operators
~~~~~~~~~~~~~~~~~

Subgrid operators are used to provide certain subgrid-wide
transformations that can be useful in certain circumstances. They are
formed by a key-value pair with syntax:

::

   <KEY>:<VALUE>

If using multiple operators, they should be comma-separated. Currently
these operators are implemented: - \*:*V* - Duplicate the subgrid data
point (there must be only one for this operator) *V* times. - +:*V* -
Increment the starting data point index of this subgrid by *V*. - N:*V*
- Normalise all data points in this subgrid by *V*.

The \* operator is typically used for normalised cross-sections, where
the total cross-section computation (a single data point) must be
duplicated *N_dat* times to correspond to the size of the ``COMMONDATA``
file. The + operator is typically used to compensate for missing
subgrids, for example when a ``COMMONDATA`` file begins with several
data points that cannot yet be computed from theory, the + operator can
be used to skip those points. The N operator is used to perform unit
conversions or the like.

Compound files and C-factors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the new dataset is a compound observable (that is, theory predictions
are a function of more than one FK-product), then one should write a
corresponding ``COMPOUND`` file as described in `Theory data
files <../data/th-data-files.rst#compound-file-format>`__. This compound
file should be stored in the APFELcomb repository under the ``compound``
directory.

C-factors should be in the format specified in `Theory data
files <../data/th-data-files.rst#cfactor-file-format>`__ and stored in
the nnpdfcpp repository under

::

   nnpdf/nnpdfcpp/data/N*LOCFAC/

directory.

Important note on subgrid ordering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the FK table consists of more than one subgrid to be merged into a
single table, then the ordering of the subgrids in their subgrid **id**
is vital. The ``merge_allgrids.py`` script will merge the subgrids in
order of their **id**. So if one is constructing an FK table for a
merged W+/W-/Z dataset, it is crucial that the ordering of the
corresponding W+/W-/Z subgrids in id matches the ordering in
``COMMONDATA``.

Important note on committing changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If one makes a modification to the ``apfelcomb.db`` database, once he is
happy with it one *must* export it to the plain-text dump file at
``db/apfelcomb.dat``. This file must then be committed. It is important
to note that the binary sqlite database is not stored in the repository.

A helper script is provided to do this. If you want to convert your
binary database to the text dump, run ``db/generate_dump.sh`` and then
commit the resulting ``apfelcomb.dat`` file.

Also, note that, if one conversely modifies the ``apfelcomb.dat`` file,
one has to delete and re-generate the sqlite database ``apfelcomb.db``
This is easily done by running ``db/generate_database.sh``.

Helper scripts
--------------

Several helper scripts are provided to make using APFELcomb easier
(particularly when generating a full set of FK tables for a particular
theory). - ``scripts/disp_grids.py`` displays a full list of
APPLgrid/FastNLO, DIS or DYP subgrids implemented in APFELcomb. -
``run_allgrids.py [theoryID] [job script]`` scans the results directory
and submits jobs for all missing subgrids for the specified theory. -
``test_submit.py`` is an example [job script] to be used for
``run\_allgrids.py``. These scripts specify how jobs are launched on a
given cluster. - ``hydra_submit.py`` is the [job script] for the HYDRA
cluster in Oxford. - ``merge_allgrids.py [theoryID]`` merges all
subgrids in the results directory for a specified theory into final FK
tables. This does not delete subgrids. - ``finalise.sh [theoryID]`` runs
C-factor scaling, copies ``COMPOUND`` files, deletes the subgrids, and
finally compresses the result into a theory.tgz file ready for upload. -
``results/upload_theories`` automatically upload to the server all the
theory.tgz files that have been generated.

Generating a complete theory
----------------------------

The general workflow for generating a complete version of a given theory
(on a cluster) cluster is then: \``\` ./run_allgrids.py
./hydra_submit.sh # Submit all APFELcomb subgrid-jobs # Once all subgrid
jobs have successfully finished ./merge_allgrids.py # Merge subgrids
into FK tables # If merging is successful ./finalise.sh # Results in a
final theory at ./results/theory\_.tgz
