#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT(libnnpdf, 1.1.0b, nathan.hartland@physics.ox.ac.uk)
AC_CONFIG_MACRO_DIR([m4])

LT_INIT()

AC_CONFIG_SRCDIR([src/common.cc])
AC_CONFIG_HEADERS([src/NNPDF/config.h])

AC_CONFIG_FILES([Makefile src/NNPDF/common.h])
AC_CONFIG_FILES([bin/Makefile bin/nnpdf-config])
AM_INIT_AUTOMAKE([subdir-objects])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
AC_SEARCH_GSL
AC_SEARCH_LHAPDF
AC_SEARCH_PYTHON
AC_SEARCH_SSE
AC_SEARCH_OPENMP
AC_SEARCH_MPI

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([gettimeofday memset mkdir pow sqrt strstr])

AC_OUTPUT
