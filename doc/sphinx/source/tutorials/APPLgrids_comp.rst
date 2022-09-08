How to benchmark and store APPLgrid/FastNLO tables
==================================================

APPLgrid and fastNLO tables produced according to :ref:`How to generate
APPLgrid and FastNLO tables <applgrids>` must be
benchmarked and properly stored before they can be used to produce FK
tables, as outlined in :ref:`How to generate and implement FK
tables <tutorialfktables>`.

Benchmark
---------

The easiest way to benchmark an APPLgrid or a FastNLO table is to
convolve it with a given set of PDFs and check that the ensuing numbers
are equivalent (within the statistical precision) to results obtained in
a completely independent way (e.g. the ones usually provided by
collaborators).

The APPLgrid software does not provide a built-in script to convolve
APPLgrid tables with a set of PDFs. Nevertheless, a simple script that
does this task is available in the

::

   external/APPLgrid_check

folder. The code can be easily compiled by doing

::

   make

or, if using conda in an environment where the user already has
``cern-root`` and ``applgrid``) by doing

.. code:: text

   conda install meson
   mkdir bld
   cd bld
   meson
   ninja

The code can be executed as

::

   ./ratio_check NAME_OF_APPLGRID.root PDF_SET

where ``PDF_SET`` must be the one used in the production of the
APPLgrids. The code will write to screen the results for all the PDF
replicas.

The FastNLO software comes with a built-in function to perform the
convolution between a FastNLO table and a PDF set. This function can be
run as

::

   fnlo-tk-cppread NAME_OF_FASTNLO.dat PDF_SET <options>

where ``options`` denote additional arguments that select the PDF member
and the scale of the FastNLO grid (see the FastNLO
`manual <https://fastnlo.hepforge.org/>`__ for details).

Important note: the numerical precision of an APPLgrid and a FastNLO
table must be sufficiently high to make it negligible in comparison to
the data and/or theoretical uncertainties. Such a precision depends on
the number of Monte Carlo events generated for each process, and must be
checked case by case. In order to increase the Monte Carlo statistics in
a reasonable amount of time, it is customary to run different
APPLgrid/FastNLO tables for the same process starting from a different
seed. The ensuing tables can then be combined with appropriate built-in
functions: - for APPLgrid

::

   applgrid-combine

-  for fastNLO

::

   fnlo-tk-merge

see the `APPLgrid <https://applgrid.hepforge.org/>`__ and
`FastNLO <https://fastnlo.hepforge.org/>`__ documentation for further
details.

.. _storage:

Storage
-------

Once the APPLgrid and/or FastNLO tables have been generated, they must
be stored in

::

   applgrids

in a folder that matches the name of the dataset used in the
``buildmaster`` implementation. Each folder contains a specific README
file with the summary information about the grid’s origin and usage.

APPLgrid and FastNLO grids are stored using `Git
LFS <https://git-lfs.github.com/>`__, which allows users to handle large
files efficiently. Git LFS can be installed with conda using:

::

   conda install git-lfs

Git LFS can then be linked to the user’s Git account and used for a
particular repository by following the instructions under `‘Getting
Started’ <https://git-lfs.github.com/>`__. If the ``applgrids``
repository already exists on the user’s system, it may need to be
re-cloned once Git LFS has been set-up to benefit from its installation.
