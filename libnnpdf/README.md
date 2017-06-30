libnnpdf 1.2.0b2
----------------------------
----------------------------

Library for core NNPDF utilities.

Dependencies
----------------------------

libnnpdf depends on the following libraries:

- pkg-config
- lhapdf
- gsl
- libarchive
- sqlite 
- yaml-cpp

please ensure to have the dependencies correctly installed and in your PATH before installing libnnpdf.

Configuration
----------------------------

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
