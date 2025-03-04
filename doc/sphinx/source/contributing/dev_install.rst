.. _source:

Development installation on Linux or macOS
===========================================

To contribute to the NNPDF code,
you will need to install a development version of the code.

Code development is carried out using Github, see :ref:`gitsection`.

While we recommend that you install ``nnpdf`` in the most comfortable way for your workflow,
installing some of the dependencies such as LHAPDF might be challenging.
Below we describe how to prepare a development environment using conda.

1. Setup a conda environment with a recent python version supported by NNPDF (see `here <https://github.com/NNPDF/nnpdf/blob/master/pyproject.toml>`_).

.. code::

   conda create -n nnpdf_dev lhapdf pandoc python==3.12 -c conda-forge


2. Clone the repository

.. code::

    git clone https://github.com/NNPDF/nnpdf.git
    cd nnpdf


Now choose how to install the dependencies of the code.
For a pip-only installation, go to 3., if you prefer to use conda go to 4.


3. Install NNPDF packages and its dependencies (make sure the conda environment
   is activated)

.. code::

  conda activate nnpdf_dev
  python -m pip install -e .

.. note::

  Following the installation steps above will set up a development
  environment that makes it possible to run and work on the nnpdf code. One
  may wish to install additional, optional, dependencies. Depending on the
  specific needs for an environment this can be dependencies that enable the
  running of the CI test, building the documentations, or performing a QED fit.

    .. code::

      python -m pip install -e .[tests,docs,qed]

4. Install the dependencies of NNPDF via conda and install the code without dependencies.

.. code::

   conda activate nnpdf_env
   conda install nnpdf --only-deps -c conda-forge -c https://packages.nnpdf.science/public
   python -m pip install -e .[tests,docs,qed] --no-deps

.. note::

   While this method of installation will ensure that the environment is self-consistent with
   respect to the packages available in conda-forge (and which might be different in pip),
   it won't keep track of changes made to the pyproject.toml file.
