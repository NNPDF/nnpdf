AC_DEFUN([AC_SEARCH_PYWRAP],[

AC_ARG_ENABLE([pywrap],
  [AC_HELP_STRING(--enable-pywrap, [Creates python wrapper.])],
  [enable_pywrap=yes], [enable_pywrap=no])
  SWIG_FOUND=no
if test x$enable_pywrap == xyes; then
  AX_PKG_SWIG(1.3.17, [], [ AC_MSG_ERROR([SWIG is required to build..]) ])
  AM_PATH_PYTHON([3.0])
  SWIG_FOUND=yes
fi

])
