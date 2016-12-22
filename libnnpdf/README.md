libnnpdf 1.1.1
----------------------------

Libary for core NNPDF utilities.

Possible configurations:
```Shell
cmake .
```
or
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
ENABLE_SAFEMODE
```
Activate python wrapper set ENABLE_PYWRAP=ON, then compile with:
```Shell
make wrapper
make wrapper-clean
```
install the libnnpdf and then make wrapper
