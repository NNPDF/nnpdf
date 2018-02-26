# libnnpdf
Library for core NNPDF utilities.

## Project summary and aim

The aim of `libnnpdf` is to provide a set of common tools shared between multiple
projects for the NNPDF Collaboration. The output of this repository is a C++ library
which can be imported and shared to other native C++ programs and python codes through
the SWIG wrapper. 

The library implements the following principle components:
- Data I/O
- PDF parametrisation
- FK tables parser
- Theoretical prediction convolutions
- PDF set handling

### Testing

Basic testing is implemented through two interfaces:
- c++ using catch: it tests c++ specific features, like copy constructors memory allocation, etc.
- swig+vp2 using pytests: it tests the python wrapper, provides utilities to dump results to disk and use a reference in future tests.

### Layout documentation

For specifications about data please check the `NNPDF/nnpdfcpp` repository in `data/doc`.
For specifications about the code design see Chapter 3 of http://arxiv.org/pdf/1509.00209.pdf
