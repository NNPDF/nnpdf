Installing the code on Linux or macOS
=====================================

Installing the NNPDF code requires a system with a recent version of either a
Linux-based operating system or macOS. There are two methods for installing the
code, both of which require conda. You can either install the code entirely with
conda or install the code from source, with the dependencies still being
installed via conda. :ref:`conda` is preferable if you simply want to run the
code, while the :ref:`source` is necessary if you want to develop the code.

.. _conda:

Installation using conda
------------------------

A helper script exists to aid the configuration using conda. To obtain it use:

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
with its own path for storing PDFs, which can be seen running ``lhapdf -–help``.
If you have an existing directory with LHAPDF grids, you may want to
either move, symlink or copy them to the new path (depending on whether
you want to keep around the older installation). The command for
symlinking would be something like:

.. code::

   ln -s <old path>/share/LHAPDF/* <new path>/miniconda3/share/LHAPDF

This will avoid symlinking the existing LHAPDF configuration, which may
be corrupted or incompatible. You should make sure only the grid directories
are transferred if you copy or move instead.


.. _source:

Installation from source
------------------------

If you intend to work on the NNPDF code, then installing from source is the
recommended installation procedure. Note that the ``binary-bootstrap.sh`` should
be downloaded and run as explained above, if the user has not already done so.

1. Setup a conda environment with a recent Python version and, if you don't have
   them yet, install ``lhapdf`` and ``pandoc`` since they are not provided
   through PyPI. Any recent python version should work, but to get an exact
   range of the supported Python versions see the `Github action workflow files
   <https://github.com/NNPDF/nnpdf/tree/master/.github/workflows>`_.

   .. code::

      conda create -n nnpdf_dev lhapdf pandoc

2. Clone the repository

   .. code::

      git clone https://github.com/NNPDF/nnpdf.git
      cd nnpdf

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

This will open a bash shell with the nnpdf environment already activated, with
all binaries and scripts from the NNPDF framework.

.. _NNPDF GitHub Packages: https://github.com/NNPDF/nnpdf/pkgs/container/nnpdf
