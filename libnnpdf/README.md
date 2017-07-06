# libNNPDF
Library for core NNPDF utilities.

## Project summary and aims

The aim of `libNNPDF` is to provide a set of common tools shared between multiple
projects for the NNPDF Collaboration. The output of this repository is a C++ library
which can be imported and shared to other native C++ programs and python codes through
the SWIG wrapper.

### Release policy

### Code development policy/rules

### Code style

### Continuous integration (CI)

### Testing

## Installation

### Binary packages

### From source

libnnpdf depends on the following libraries:

- pkg-config
- lhapdf
- gsl
- libarchive
- sqlite
- yaml-cpp

please ensure to have the dependencies correctly installed and in your PATH before installing libnnpdf.

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
You can control the optional flags with ccmake or from cmd line, the most relevant flags are:

```Shell
CMAKE_INSTALL_PREFIX
ENABLE_OPENMP
ENABLE_OPENMPI
ENABLE_PYWRAP
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

## Documentation
