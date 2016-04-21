AC_PREREQ([2.59])

define([gitversion], [esyscmd([sh -c "git rev-parse --short HEAD | tr -d '\n'
"])])

AC_INIT([nnpdfsrc],[gitversion],[stefano.carrazza@mi.infn.it])
AC_CONFIG_MACRO_DIR([m4])
AM_MAINTAINER_MODE

AC_CONFIG_AUX_DIR([autoconfig])
AC_CONFIG_HEADERS([autoconfig/config.h])
AC_CONFIG_FILES([common/inc/svn.h])

AM_INIT_AUTOMAKE([subdir-objects])
AC_PREFIX_DEFAULT(`pwd`/../nnpdfbuild)

LT_INIT

# checks for programs
AC_PROG_CXX
AC_PROG_CC

# check for root
AC_SEARCH_GSL
AC_SEARCH_ROOT
AC_SEARCH_LHAPDF
AC_SEARCH_PYTHON
AC_SEARCH_NNPDF
AC_SEARCH_SQLITE3
AC_SEARCH_APFEL
AC_SEARCH_YAML

# Checks for header files
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([projects/Makefile])
AC_CONFIG_FILES([testcode/Makefile])

AC_OUTPUT
