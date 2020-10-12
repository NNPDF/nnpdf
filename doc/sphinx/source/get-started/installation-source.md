```eval_rst
.. _source:
```
## Installation from source

If you intend to work on the NNPDF code, then building from source is the recommended installation procedure. However, you can still use conda to get all the dependencies and setup the validphys and C++ development environment. Further information is available in the [vp-guide](https://data.nnpdf.science/validphys-docs/guide.html#development-installs). Note that the `binary-bootstrap.sh` should be downloaded and run as explained above, if the user has not already done so.

1. Create an NNPDF developer environment `nnpdf-dev` and install all relevant dependencies using

		conda create -n nnpdf-dev
		conda activate nnpdf-dev
		conda install --only-deps nnpdf

	Note that the user should be in the conda environment `nnpdf-dev` whenever they wish to work on NNPDF code. The conda environment can be exited using `conda deactivate`.

	<details>
      <summary>Important note for Mac users</summary>

      If you are a macOS user, you will need to download the
      [Mac Software Development Kit](https://github.com/phracker/MacOSX-SDKs) or
      SDK for short. This is necessary to get the correct C compiler. Assuming
      that you already have access to the [server](NNPDF-server), you can
      download the version of the SDK used by the [Continuous Integration](CI)
      system by doing

		curl -L -O https://data.nnpdf.science/MacOSX10.9.sdk.tar.xz

	  You can then unpack it into your root conda directory by running

		tar xfz MacOSX10.9.sdk.tar.xz -C <path_to_root_conda_directory>

	  where you can find `<path_to_root_conda_directory>` by typing
	  `echo $CONDA_PREFIX` when your base conda environment is activated. You
	  should then export the following path

		export CONDA_BUILD_SYSROOT=<path_to_root_conda_directory>/MacOSX10.9.sdk

	  which you may wish to write to one of your `~/.bashrc` or
	  `~/.bash_profile` scripts so that the SDK is easily accessible from the
	  shell.

    </details>

2. Install the appropriate C++ compilers using
		
		conda install gxx_linux-64 

	macOS users should replace `gxx_linux-64` with `clangxx_osx-64`.

3. Ensure that the NNPDF repositories `nnpdf` and `apfel` are in the `nnpdfgit` directory. These are required to be able to run fits and can be obtained respectively by

		git clone git@github.com:NNPDF/nnpdf.git
		git clone https://github.com/scarrazza/apfel.git

4. Obtain the dependencies of the code you want to build. Where to find those depends on the particular code. For example, something linking to `libnnpdf` will likely require `pkg-config`. Projects based on `autotools` (those that have a `./configure` script) will additionally require `automake` and `libtool`. Similarly projects based on `cmake` will require installing the `cmake` package. In the case of `nnpdf` itself, the build dependencies can be found in  `<nnpdf git root>/conda-recipe/meta.yaml`. We have to install the remaining ones manually:

		conda install pkg-config swig=3.0.10 cmake

5. We now need to make the installation prefix point to our `nnpdf-dev` environment, this can be done using:

		$CONDA_PREFIX=~/miniconda3/envs/nnpdf-dev/

	this assumes `miniconda3` is installed in the default place which is the home directory.

6. Navigate to the `nnpdf` directory obtained from the Github repository and create a folder called `conda-bld` by
		
		nnpdf$ mkdir conda-bld
		nnpdf$ cd conda-bld

	Note that it is important that for the following step to be executed while the user is in the `nnpdf-dev` conda environment. The project can be built using:

		nnpdf/conda-bld$ cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX

7. When the user wishes to work on the NNPDF code, they should do so in, for example, `'/nnpdfgit/nnpdf/libnnpdf'`. To compile the code navigate to the `conda-bld` directory created above and run

		make
		make install
