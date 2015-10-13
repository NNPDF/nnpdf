AC_DEFUN([AC_SEARCH_MPI],[

AC_ARG_ENABLE([openmpi],
  [AC_HELP_STRING(--enable-openmpi, [Enable MPI convolution.])],
  [enable_openmpi=yes], [enable_openmpi=no])
if test x$enable_openmpi == xyes; then
  AX_EXT
  AC_SUBST(LIBNNPDF_HAVE_MPI, ["#define OPENMPI"])
fi

])
