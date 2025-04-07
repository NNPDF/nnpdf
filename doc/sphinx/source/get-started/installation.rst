Installing the code on Linux or macOS
=====================================

Installing the NNPDF code requires a system with a recent version of python and
either a Linux-based operating system or macOS.
At least 4 GB of storage are needed in order to download all data and theory ingredients
necessary for a fit or most analyses.
Some specific features such as theory uncertainties or QED fits might require more space.

The code can be installed either with :ref:`conda <condainstall>` or with :ref:`pip <pip>`,
the former offering a "batteries included" approach where non-python software
such as LHAPDF and pandoc will be installed automatically.

If you plan to contribute to the development of NNPDF, please see :ref:`Source`.


.. _condainstall:

Installation using conda
------------------------

The ``nnpdf`` package is available from `conda-forge <https://anaconda.org/conda-forge/nnpdf>`_,
therefore it is possible to install the most recent release of ``nnpdf`` as well as its dependencies
by running the following command (optionally, with ``--override-channels`` to ensure all dependencies are from ``conda-forge``.

::

  conda create -n environment_nnpdf nnpdf -c conda-forge [--override-channels]
  conda activate environment_nnpdf


This will create a new conda environment ``environment_nnpdf`` in which the conda package for ``nnpdf``,
downloaded from our public repository, will be installed.
You are now ready to use the NNPDF code! Check out the :ref:`Tutorials`.

.. note::

   Make sure to be using a recent version of `conda <https://docs.anaconda.com/miniconda/install/>`_. These instructions have been tested with conda 25.



.. _pip:

Installation using pip
----------------------

Most NNPDF packages and its dependencies are available in the `PyPI repository <https://pypi.org>`_.
While the fitting code is currently not available, it can be installed directly from the git repository:

::

  python -m venv environment_nnpdf
  . environment_nnpdf/bin/activate
  python -m pip install git+https://github.com/NNPDF/nnpdf.git@4.0.10


.. warning::

   When you install using pip, non-python codes such as LHAPDF and pandoc won't be installed automatically and neeed to be manually installed in the environment. If using python 3.9, make sure it is newer than ``3.9.2`` (see issue `here <https://github.com/NNPDF/reportengine/pull/69>`_)


Shared data
-----------

By default, shared data in NNPDF will be saved to the enviroment's ``share`` path.
In the case of conda, it defaults to ``${CONDA_PREFIX}/share/NNPDF``.
It is possible to configure where to download theory and results using a ``nnprofile`` file as described in :ref:`nnprofile`.




Using the code with docker
--------------------------

We provide docker images for tag release of the code using GitHub Packages. The
docker images contain a pre-configured linux environment with the NNPDF
framework installed with the specific tag version. The code is installed using
miniconda3.

Please refer to the download and authentication instructions from the `NNPDF GitHub Packages`_.

In order to start the docker image in interactive mode please use docker
standard syntax, for example:

.. code::

    docker run -it ghcr.io/nnpdf/nnpdf:<tag_version> bash

This will open a bash shell with the ``nnpdf`` environment already activated, with
all binaries and scripts from the NNPDF framework.

.. _NNPDF GitHub Packages: https://github.com/NNPDF/nnpdf/pkgs/container/nnpdf


.. _dependencies:

Dependencies and requirements
-----------------------------

The NNPDF framework would not be possible with a number of dependencies.
These are automatically installed when using conda and the full list can be consulted in the
`conda recipe <https://github.com/NNPDF/nnpdf/blob/master/conda-recipe/meta.yaml#L20>`_ available in the repository.

Below we list some of the most relevant external libraries than underpin the different aspects of this framework.

::

  LHAPDF

  keras
  tensorflow

  numpy
  pandas
  scipy
  matplotlib
  numba
