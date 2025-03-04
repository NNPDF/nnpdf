Installing the code on Linux or macOS
=====================================

Installing the NNPDF code requires a system with a recent version of python and
either a Linux-based operating system or macOS.
The code can be installed either with :ref:`conda <condainstall>` or with :ref:`pip <pip>`,
the former offering a "batteries included" approach where non-python software
such as LHAPDF and pandoc will be installed automatically.

If you don't have a recent version of conda installed, please use the
following :ref:`bootstrap-installation`
to install conda as well as set up relevant channels.

If you plan to contribute to the development of NNPDF, please see :ref:`Source`.


.. _condainstall:

Installation using conda
------------------------

It is possible to install the most recent release of ``nnpdf`` as well as its dependencies
by running the following command:

::

  conda create -n environment_nnpdf nnpdf -c conda-forge -c https://packages.nnpdf.science/public [--override-channels]
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
  python -m pip install git+https://github.com/NNPDF/nnpdf.git@4.0.9


.. warning::

   When you install using pip, non-python codes such as LHAPDF and pandoc won't be installed automatically and neeed to be manually installed in the environment.


.. _bootstrap-installation:

Installation using the bootstrap script
---------------------------------------


If you don't have a recent version of conda installed,
a helper script exists to aid the configuration using conda. To obtain it use:

::

       mkdir nnpdfgit
       cd nnpdfgit
       git clone git@github.com:NNPDF/binary-bootstrap.git

-  Execute the script

   ::

        ./binary-bootstrap/bootstrap.sh

-  **Path**: the conda installer will ask to add the conda bin path to
   the default *$PATH* environment variable (by editing your .bashrc
   file). Confirm this unless you know that you have a specific reason
   not to. Note that newer versions of conda give the option of having
   conda available, but not any environment (which you have to enable
   explicitly by either having conda activate in .bashrc or typing it
   each time you want to use the environment). On remote machines, the
   addition to .bashrc should read as follows

   ::

        if shopt -q login_shell; then
            . <path-to-conda>/etc/profile.d/conda.sh
            conda activate
        fi

the if condition is important because conda activate prints to the
standard output, which interferes with commands like scp and rsync.

-  Note that the script may ask you to perform some actions manually (
   e.g. it will not overwrite your existing conda configuration). Please
   pay attention to the output text of the script.

Installing the NNPDF code
~~~~~~~~~~~~~~~~~~~~~~~~~

After the helper script has run, navigate to the miniconda3 installation
directory, by default this is ``~/miniconda3``, and run the command

.. code::

       . ./etc/profile.d/conda.sh
       conda activate
       conda install nnpdf

**Note:** The installer will set up its own version of the LHAPDF code,
with its own path for storing PDFs, which can be seen running ``lhapdf --help``.
If you have an existing directory with LHAPDF grids, you may want to
either move, symlink or copy them to the new path (depending on whether
you want to keep around the older installation). The command for
symlinking would be something like:

.. code::

   ln -s <old path>/share/LHAPDF/* <new path>/miniconda3/share/LHAPDF

This will avoid symlinking the existing LHAPDF configuration, which may
be corrupted or incompatible. You should make sure only the grid directories
are transferred if you copy or move instead.


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
