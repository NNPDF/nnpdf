# libnnpdf
Library for core NNPDF utilities.

## Project summary and aim

The aim of `libnnpdf` is to provide a set of common tools shared between multiple
projects for the NNPDF Collaboration. The output of this repository is a C++ library
which can be imported and shared to other native C++ programs and python codes through
the SWIG wrapper. 

The library implements basically the following components:
- Data I/O
- NN parameterization
- FK tables parser
- Theoretical prediction convolutions
- PDF set handlings

### Release and Tag policy

The library is tagged and released when a major and stable status is achieved. 
Tags and releases do not necessarily follows the NNPDF releases.

### Code development policy/rules

Developers must never commit code structure modifications to master. The development pattern should follow these rules:
- Open an issue explaining your bug or feature request. If you report a bug, post information to reproduce it.
- The resolution of issues must be performed in a new branch through a pull request.
- If you have already a local version of the code that you would like to merge in the master, open a pull request.
- The pull request must be reviewed by at least 2 core developers.

### Code style

Originally the code of this library was written at the beginning of 2012 so at that time C++11 was not used as default. 
During the last months the library is receiving constant feedback and requests to improve the data types, so a code modernization is in place.

### Continuous integration (CI)

We implement CI at CERN's gitlab. After each commit `libnnpdf` is compiled at the CI nodes. The CI is enabled for all branches.
Zahari has a private TravisCI instance used to build binary packages (conda packages) for this project.

### Testing

Basic testing is implemented through two interfaces:
- c++ using catch: it tests c++ specific features, like copy constructors memory allocation, etc.
- swig+vp2 using pytests: it tests the python wrapper, provides utilities to dump results to disk and use a reference in future tests.

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

A detailed validphys2 guide including conda installation instructions can be found here:

http://pcteserver.mi.infn.it/~nnpdf/validphys-docs/guide.html

### From source

`libnnpdf` depends on the following libraries:

- pkg-config
- lhapdf
- gsl
- libarchive
- sqlite
- yaml-cpp
- python 3.x (optional)

please ensure to have the dependencies correctly installed and in your PATH before installing libnnpdf.
The exact or minimal version requirements for each package is summarized in https://github.com/NNPDF/libnnpdf/blob/master/conda-recipe/meta.yaml.

#### Configurations

Possible configurations:

```Shell
cmake .

```
or (recommended):

```Shell
mkdir build
cd build
cmake ..

```
You can control the optional flags with `ccmake` or from cmd line, the most relevant flags are:

```Shell
CMAKE_INSTALL_PREFIX
ENABLE_OPENMP
ENABLE_OPENMPI
ENABLE_PYWRAP
ENABLE_TESTS
PROFILE_PREFIX
```

On the command line, options are controlled appending a `-D` flag. For
example:

```
cmake .. -DENABLE_PYWRAP=on
```

To compile the Python wrappers, set ENABLE_PYWRAP=ON as above, and
after libnnpdf is installed, then compile with:

```Shell
make wrapper
make wrapper-clean
```

This will also install the wrapper in the current Python environment.

If `PROFILE_PREFIX` sets paths and urls required to access the nnpdfcpp data, theory and results resources.
If this flag is empty the paths will be set to the libnnpdf installation folder.

## Documentation

### Code documentation

THe code is documented with Doxygen, if you find methods or classes not fully documented open a issue request.

### Layout documentation

For specifications about data please check the `nnpdfcpp` repository in `data/doc`.
For specifications about the code design see Chapter 3 of http://arxiv.org/pdf/1509.00209.pdf
