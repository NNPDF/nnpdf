# Installation

## Installation using conda
### Adding RSA key to Github

To clone the NNPDF Github repositories, you must add a public RSA key to Github.

* Checking for an existing RSA key

	Firstly, check if your machine has already generated an RSA key. This can be done by the following:

		cd ~/.ssh
		ls -a

	and verify the files `id_rsa` and `id_rsa.pub` exist. An alternative is to perform `find ~ -name 'id_rsa.*'` and verify `id_rsa.pub` is found.

* Generating an RSA key

	If no RSA key exists use the command 

		ssh-keygen

	and follow the instructions to generate an RSA key.

* Adding the RSA key to Github

	Copy the public RSA key stored in `id_rsa.pub` to your clipboard using

		cat ~/.ssh/id_rsa.pub

	and copying the output.

	The public RSA key must be pasted to Github [here](https://github.com/settings/keys) and clicking on the 'New SSH key' button. Enter an appropriate title and paste the RSA Key. 

### Obtain the helper script

A helper script exists to aid the configuration. To obtain it use:

		mkdir nnpdfgit
		cd nnpdfgit
		git clone git@github.com:NNPDF/binary-bootstrap.git

* Execute the script
	
		./binary-bootstrap/bootstrap.sh

	The script will ask you for the password of the NNPDF private repository, which can be found [here](https://www.wiki.ed.ac.uk/pages/viewpage.action?pageId=292165461)

* **Path**: the conda installer will ask to add the conda bin path to the default *$PATH* environment variable (by editing your  .bashrc file). Confirm this unless you know that you have a specific reason not to. Note that newer versions of conda give the option of having conda available, but not any environment (which you have to enable explicitly by either having  conda activate in .bashrc or typing it each time you want to use the environment). On remote machines, the addition to .bashrc should read as follows

		if shopt -q login_shell; then
			. <path-to-conda>/etc/profile.d/conda.sh
			conda activate
		fi

the if condition is important because conda activate prints to the standard output, which interferes with commands like scp and rsync.

* Note that the script may ask you to perform some actions manually ( e.g. it will not overwrite your existing conda configuration). Please pay attention to the output text of the script.

### Installing the NNPDF code

After the helper script has run, navigate to the miniconda3 installation directory, by default this is `~/miniconda3`, and run the command

		. ./etc/profile.d/conda.sh
		conda activate
		conda install nnpdf

This will provide functioning C++ and Python executables.

**Note:** The installer will set up it’s own version of the LHAPDF code, with its own path for storing PDFs, which can be seen running lhapdf --help. If you have an exiting folder with LHAPDF grids, you may want to either move, symlink or copy them to the new path (depending on whether you want to keep around the older installation). The command for symlinking would be something like:

	ln -s <old path>/share/LHAPDF/* <new path>/miniconda3/share/LHAPDF

This will avoid symlinking the existing LHAPDF configuration, which may be corrupted or incompatible. You should make sure only the grid folders are transferred if you copy or move instead.

## Installation from source
If you intend to work on the NNPDF code, then building from source is the recommended installation procedure. However, you can still use conda to get all the dependecies and setup the validphys and C++ development environment. Further information is available in the [vp-guide](https://data.nnpdf.science/validphys-docs/guide.html#development-installs). Note that the `binary-bootstrap.sh` should be downloaded and run as explained above, if the user has not already done so.

1. Create an NNPDF developer environment `nnpdf-dev` and install all relevant dependencies using
		
		conda create -n nnpdf-dev 
		conda activate nnpdf-dev
		conda install --only-deps nnpdf

	Note that the user should be in the conda environment `nnpdf-dev` whenever they wish to work on NNPDF code. The conda environment can be exited using `conda deactivate`.

2. Install the appropriate C++ compilers using
		
		conda install gxx_linux-64 

	macOS users should replace `gxx_linux-64` with `clangxx_osx-64`.

3. Ensure that the NNPDF repositories `nnpdf` and `apfel` are in the `nnpdfgit` directory. These are required to be able to run fits and can be obtained respectively by

		git clone git@github.com:NNPDF/nnpdf.git
		git clone https://github.com/scarrazza/apfel.git

4. Obtain the dependencies of the code you want to build. Where to find those depends on the particular code. For example, something linking to `libnnpdf` will likely require `pkg-config`. Projects based on `autotools` (those that have a `./configure` script) will additionally require `automake` and `libtool`. Similarly projects based on `cmake` will require installing the `cmake` package. In the case of `nnpdf` itself, the build dependencies can be found in  `<nnpdf git root>/conda-recipe/meta.yaml`. We have to install the remaining ones manually:

		conda install pkg-config swig=3.0.10 cmake

5. We now need to make the installation prefix point to our `nnpdf-dev` environment, this can be done using:

		$CONDA_PREFIX==~/miniconda3/envs/nnpdf-dev/

	this assumes `miniconda3` is installed in the default place which is the home directory.

6. Navigate to the `nnpdf` directory obtained from the Github repository and create a folder called `conda-bld` by
		
		nnpdf$ mkdir conda-bld
		nnpdf$ cd conda-bld

	Note that it is important that for the following step to be executed while the user is in the `nnpdf-dev` conda environment. The project can be built using:

		nnpdf/conda-bld$ cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX

7. When the user wishes to work on the NNPDF code, they should do so in, for example, `'/nnpdfgit/nnpdf/libnnpdf'`. To compile the code navigate to the `conda-bld` directory created above and run

		make
		make install
