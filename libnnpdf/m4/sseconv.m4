AC_DEFUN([AC_SEARCH_SSE],[

AC_ARG_ENABLE([safemode],
  [AC_HELP_STRING(--enable-safemode, [Disable SSE convolution (for debug).])],
  [enable_safemode=no], [enable_safemode=yes])
if test x$enable_safemode == xyes; then
  AX_EXT
  AC_SUBST(LIBNNPDF_HAVE_SSE, ["#define SSE_CONV"])
fi

])
