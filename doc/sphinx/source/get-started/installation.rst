Installing the code
===================

There are two methods for installing the code, both of which require
conda. You can either install the code entirely with conda or install
the code from source, with the dependencies still being installed via
conda. :ref:`conda` is preferable if you simply want to run the
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

This will provide functioning C++ and Python executables.

**Note:** The installer will set up its own version of the LHAPDF code,
with its own path for storing PDFs, which can be seen running ``lhapdf –help``.
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

If you intend to work on the NNPDF code, then building from source is
the recommended installation procedure. However, you can still use conda
to get all the dependencies and setup the validphys and C++ development
environment. Note
that the ``binary-bootstrap.sh`` should be downloaded and run as
explained above, if the user has not already done so.

**NOTE:** For installation on M1/M2 Macs, please see the :ref:`M1` section.

1. Create an NNPDF developer environment ``nnpdf-dev`` and install all
   relevant dependencies using

   .. code::

       conda create -n nnpdf-dev
       conda activate nnpdf-dev
       conda install --only-deps nnpdf

   Note that the user should be in the conda environment ``nnpdf-dev``
   whenever they wish to work on NNPDF code. The conda environment can
   be exited using ``conda deactivate``.

   .. note::

        If you are a macOS user, you will need to download the `Mac Software
        Development Kit`_ or SDK for short. This is necessary to get the
        correct C compiler. The `anconda documentation`_ explains in more
        detail why you need this file, and why they cannot include it with
        the compilers by default.

        You can check which version of SDK is currently being used by the
        :ref:`CI` system by checking the ``MACOS_SDK_URL``
        and ``MACOS_SDK_FILE`` variables inside ``.travis.yml``. At the time
        of writing this documentation, the version used is 10.9 but the user
        is advised to check in case the documentation has become out of sync
        with the CI configuration. Once you know the URL of the SDK file, you
        can download it from the commandline using ``curl``, e.g.:

        .. code::

            curl -L -O https://github.com/phracker/MacOSX-SDKs/releases/download/11.3/MacOSX10.9.sdk.tar.xz

        You can then unpack it into your root conda directory by running

        .. code::

            tar xfz MacOSX10.9.sdk.tar.xz -C <path_to_root

        where you can find ``<path_to_root_conda_directory>`` by typing
        ``echo $CONDA_PREFIX`` when your base conda environment is activated. You
        should then export the following path

        .. code::

            export CONDA_BUILD_SYSROOT=<path_to_root_conda_directory>/MacOSX10.9.sdk

        which you may wish to write to one of your ``~/.bashrc`` or
        ``~/.bash_profile`` scripts so that the SDK is easily accessible from the
        shell.

2. Install the appropriate C++ compilers using

   .. code::

       conda install gxx_linux-64

   macOS users should replace ``gxx_linux-64`` with ``clangxx_osx-64``.

3. Ensure that the NNPDF repositories ``nnpdf`` and ``apfel`` are in the
   ``nnpdfgit`` directory. These are required to be able to run fits and
   can be obtained respectively by

   .. code::

       git clone git@github.com:NNPDF/nnpdf.git
       git clone https://github.com/scarrazza/apfel.git

4. Obtain the dependencies of the code you want to build. Where to find
   those depends on the particular code. For example, something linking
   to ``libnnpdf`` will likely require ``pkg-config``. Projects based on
   ``autotools`` (those that have a ``./configure`` script) will
   additionally require ``automake`` and ``libtool``. Similarly projects
   based on ``cmake`` will require installing the ``cmake`` package. In
   the case of ``nnpdf`` itself, the build dependencies can be found in
   ``<nnpdf git root>/conda-recipe/meta.yaml``. We have to install the
   remaining ones manually:

   .. code::

       conda install pkg-config swig cmake
      
   When working on a Linux system it is `currently 
   <https://github.com/NNPDF/nnpdf/pull/1280>`_ also needed to run
   
   .. code::
   
       conda install sysroot_linux-64=2.17

5. We now need to make the installation prefix point to our
   ``nnpdf-dev`` environment. Fortunately, when you activate the environment,
   the location is saved to the environment variable ``CONDA_PREFIX``, e.g.

   .. code::

       $ conda activate nnpdf-dev
       $ echo $CONDA_PREFIX
       /home/miniconda3/envs/nnpdf-dev/

6. Navigate to the ``nnpdf`` directory obtained from the Github
   repository and create a directory.

   .. note::

        The directory name is unimportant,
        its role is to contain all of the build files, separately from the source
        files - we will refer to it as the build directory. A clean
        build and installation can always be achieved by
        deleting the contents of the build directory (or even creating a new one)
        and re-running ``cmake``.

   For this example we have created a directory called ``conda-bld`` by

   .. code::

       nnpdf$ mkdir conda-bld
       nnpdf$ cd conda-bld

   Note that it is important that for the following step to be executed
   while the user is in the ``nnpdf-dev`` conda environment. The project
   can be built using:

   .. code::

       nnpdf/conda-bld$ cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX

7. When the user wishes to work on the NNPDF code, they should do so in,
   for example, ``/nnpdfgit/nnpdf/libnnpdf``. To compile the code
   navigate to the build directory created above and run

   .. code::

       rm -r ./*
       cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX
       make
       make install

   If you are reinstalling the code using the same build directory, it is
   recommended to remove all files from the build directory as is shown
   in the example above.

.. _here: https://github.com/settings/keys
.. _Mac Software Development Kit: https://github.com/phracker/MacOSX-SDKs
.. _anconda documentation: https://docs.conda.io/projects/conda-build/en/latest/resources/compiler-tools.html#macos-sdk


.. _M1:

Installation from source on M1/M2 Macs
--------------------------------------

Installation on M1/M2 Macs is not directly supported, so everything needs to be
built manually. The following steps are required:

1. Clone the repositories

  .. code::

      mkdir nnpdfgit
      cd nnpdfgit
      git clone git@github.com:NNPDF/nnpdf.git
      git clone git@github.com:NNPDF/binary-bootstrap.git
      git clone https://github.com/scarrazza/apfel.git

2. Execute binary bootstrap to set the channels in ``.condarc``

   .. code::

      ./binary-bootstrap/bootstrap.sh

3. Setup conda environment using python 3.9

   .. code::

      conda create -n nnpdf-dev python=3.9
      conda activate nnpdf-dev

4. Install ARM compiler

   .. code::

      conda install clangxx_osx-arm64

5. LHAPDF

   Download version 6.4.0 and decompress

   .. code::

       wget -O LHAPDF-6.4.0.tar.gz https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.4.0.tar.gz
       tar -xzvf LHAPDF-6.4.0.tar.gz
       rm LHAPDF-6.4.0.tar.gz
       cd LHAPDF-6.4.0

   Regenerate the configuration files, configure the build with python disabled, compile and 
   install. You may need to `brew install automake` first:

   .. code::

      autoreconf -f -i
      ./configure --prefix=$CONDA_PREFIX --disable-python
      make -j
      make install

   Install the python wrapper

   .. code::

      cd wrappers/python
      pip install -e .

   Test

   .. code::

      lhapdf install CT18NNLO
      python -c "import lhapdf"

6. Apfel


   First, we need to install some dependencies:

   .. code::

      conda install pkg-config swig cmake

   Then build it

   .. code::

      cd ../../../apfel
      autoreconf -f -i
      PYTHON=$(which python) ./configure --prefix=$CONDA_PREFIX 
      make clean
      make -j
      make install

7. validphys

   First install reportengine and validobj, then validphys itself:

   .. code::

      pip install reportengine validobj
      cd ../nnpdf/validphys2
      pip install -e .

8. nnpdf

   Install other dependencies

   .. code::

      conda install libarchive sqlite gsl yaml-cpp

   Run cmake in nnpdf/build directory:

   .. code::

      cd ..
      mkdir build
      cd build
      cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX

   Edit the file ``nnpdfgit/nnpdf/CMakeLists.txt`` :

    - on line 8 change the option to true, so it says:

      .. code::

         SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

    - comment out line 58 (:code:`set(LIBNNPDF_HAVE_SSE "#define SSE_CONV")`)

    - line 104 should read:

      .. code::

         set(DEFAULT_CXX_OPTIONS "-Wall -Wextra -fvisibility-inlines-hidden -fmessage-length=0 -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe")

      (so delete ``-march=nocona -mtune=haswell``).

   Then make:

   .. code::

      make -j
      make install

   Install remaining packages

   .. code::

      pip install seaborn prompt_toolkit scipy psutil hyperopt

9. Install tensorflow

   Not specifying versions will install at the time of writing macos 2.12.0 and metal 0.8.0, which both work.
   They only give warnings on the optimizers, that the legacy versions are faster.
   If you want an older version, macos 2.9.2 and metal 0.5.0 are also tested to work.

   .. code::

      conda install -c apple tensorflow-deps
      pip install tensorflow-macos==2.9.2
      pip install tensorflow-metal==0.5.0

10. Test

   .. code::

      cd nnpdf/n3fit/runcards/examples
      vp-setupfit Basic_runcard.yml
      n3fit Basic_runcard.yml 1
      evolven3fit Basic_runcard 1

   With these settings tensorflow will run by default on GPU which makes
   the fit run very slow. To disable the GPU, type the following command:
   
   .. code::

      export CUDA_VISIBLE_DEVICES=0
      
   or insert the following line in the ``set_initial_state`` function in ``n3fit/src/n3fit/backends/keras_backend/internal_state.py``:

   .. code::

      tf.config.set_visible_devices([], 'GPU')

   And to use legacy optimizers, you only need to change one line in ``n3fit/src/n3fit/backends/keras_backend/MetaModel.py``:

   .. code::

      # from tensorflow.keras import optimizers as Kopt
      import tensorflow.keras.optimizers.legacy as Kopt

   With both these tweaks, and the latest tensorflow versions, the basic runcard with 1 replica should take about 30 seconds.

.. _docker:

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
