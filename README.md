# NNPDF/nnpdf 
The main repository with the NNPDF framework. 

**Table of Contents**
  * [Project summary and aims](#project-summary-and-aims)
  * [Installation](#installation)
    * [Binary packages](#binary-packages)
    * [From source](#from-source)    
  * [Using the code](#using-the-code)
 
## Project summary and aim

This project is contains the following components:
- libnnpdf: common core objects shared among projects
- nnpdfcpp: provides a set of programs used in the NNPDF fitting framework. 
- validphys2: the postfit analysis framework

### Release and Tag policy

The library is tagged and released when a major and stable status is achieved. 
Tags and releases follows the NNPDF releases. This project requires compatible versions of `libnnpdf` see binary installation for details.

### Code development policy/rules

Developers must never commit code structure modifications to master. The development pattern should follow these rules:
- Open an issue explaining your bug or feature request. If you report a bug, post information to reproduce it.
- The resolution of issues must be performed in a new branch through a pull request.
- If you have already a local version of the code that you would like to merge in the master, open a pull request.
- The pull request must be reviewed by at least 2 core developers.

### Code style

Originally the code of this library was written at the beginning of 2012 so at that time C++14 was not used as default. 
During the last months the library is receiving constant feedback and requests to improve the data types, so a code modernization is in place.

### Testing

Testing is actually not implemented in the current repository, but there are plans for that.


## Installation

We provide two installation methods: binary packages and from source.
If do not plan to apport modifications to the library the recommended installation system is through binary packages.

### Binary packages

The master version of `libnnpdf` and its dependencies can be obtained in binary format, as a conda package. A bootstrap script exist to aid the configuration. Simply clone its repository and execute it:
```Shell
git clone git@github.com:NNPDF/binary-bootstrap.git
./binary-botstrap/bootstrap.sh
```
The script will ask for the password of the private NNPDF repositories. It is:
```
BifaSali9
```
Once configured, it is possible to install libnnpdf or apfel by simply:
```Shell
conda install nnpdf validphys2
```

Which will pull also LHAPDF, libnnpdf, apfel and all the other dependencies.

When the packages are installed, the necessary binaries are added to
the `bin/` directory of the corresponding conda environment (which is
typically in the `PATH`). Users can simply run `filter`, `nnfit`, or
`postfit2` from any directory. Once you have the binary package
installed, you should not need the git repositories to run a fit.

By default, the data from theory and experiment is located under:
```
<conda root>/share/NNPDF/data
```

and the fits will be stored in:
```
<conda root>/share/NNPDF/results
```

These paths can be changed by tweaking `nnprofile.yaml` as
described in [NNPDF paths and URLS](#nnpdf-paths-and-urls).

A detailed validphys2 guide including conda installation instructions can be found here:

http://pcteserver.mi.infn.it/~nnpdf/validphys-docs/guide.html

### From source

`nnpdfcpp` depends on the following libraries in order to compile all programs:

- pkg-config
- lhapdf
- gsl
- libarchive
- sqlite
- yaml-cpp
- libnnpdf
- CERN-ROOT (for validphys)
- APFEL
- fiatlux

please ensure to have the dependencies correctly installed and in your PATH before compiling nnpdfcpp.
The exact or minimal version requirements for each package is summarized in https://github.com/NNPDF/nnpdfcpp/blob/master/conda-recipe/meta.yaml.

#### Compiling the code

Several options are available, please use ccmake:
```Shell
mkdir build
cd build
cmake ..
ccmake .
```
Make sure you have set correctly the `CMAKE_INSTALL_PREFIX`, this path is used when `make install` is invoked.

After running cmake proceed with `make` followed by `make install`. The later will copy binaries and scripts to the `CMAKE_INSTALL_PREFIX` that you selected, while the content of `nnpdfcpp/data` folder will be copied to the `data_path` set in nnprofile.yaml.

### NNPDF paths and URLS

The paths that various codes (such as `nnfit` and `validphys`) will
use to find and write resources, as well as the URLS to upload and
download them are defined in a `nnprofile.yaml` file. By default, it
is stored in the `libnnpdf` install prefix, under `<libnnpdf install
prefix>/share/NNPDF/nnprofile.yaml`. For binary packages, the
`libnnpdf` install prefix is simply the path of the conda environment
where the packages is installed.
The paths and URLs can be
modified: This can be useful to make the code work under specific
cluster configurations, for example to avoid excessive I/O in NFS
mounts. However, do not do it by editing the `nnprofile.yaml` file in the
default location, since it will be overwritten every time that
`libnnpdf` is installed. 
Instead copy it to some other location, make
the changes you wish, and define an `NNPDF_PROFILE_PATH` environment
variable pointing to your modified file. For example, you could write
in your `.bashrc`:
```shell
export NNPDF_PROFILE_PATH=/home/user/mynnprofile.yaml
```

## Using the code

### Runcard

The runcard is written in YAML. The interpretation of the entries should be self-explanatory. The runcard is the unique identifier of a fit, it is also the only required configuration input required for many programs of this repository.

### Workflow

0. compile the code and install

1. Create a runcard by taking as template one of the files in `nnpdfcpp/config`.

2. Go to `<data_path>` and download the specific theory folder using `disp_theory.py` and passing the `theoryid` of your runcard as theory number.

3. Filter the data:
```Shell
filter <runcard>.yml
```
this command will generate a `<runcard>` folder in the current directory with a copy of the original YAML runcard and its md5 key.

4. All programs take the `<runcard>` folder as input, e.g.
```Shell
nnfit <replica_number> <runcard_folder>
```
where replica_number goes from 1-n.

### Code documentation

THe code is documented with Doxygen, if you find methods or classes not fully documented open a issue request.

### Layout documentation

For specifications about data please check `data/doc`.
For specifications about the code design see Chapter 3 of http://arxiv.org/pdf/1509.00209.pdf
