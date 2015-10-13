#AC_SEARCH_NNPDF(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_NNPDF],[
  AC_ARG_WITH([nnpdf], AC_HELP_STRING(--with-nnpdf, [path to NNPDF library and header files]))

  AC_SUBST(NNPDF_SVN_VERSION, ["#define SVN_REV \"$VERSION\""])
  
  ## Use a specified --with-nnpdf arg to set basic paths, if provided
  NNPDFCONFIG_PATH=$PATH
  if test -e "$with_nnpdf"; then
    NNPDFCONFIG_PATH=$with_nnpdf/bin:$NNPDFCONFIG_PATH
    NNPDFPATH="$with_nnpdf"
    NNPDFINCPATH="$NNPDFPATH/include"
    NNPDF_CPPFLAGS="-I$NNPDFINCPATH"
    NNPDF_CXXFLAGS=""
    NNPDF_LDFLAGS="-L$NNPDFPATH/lib -lnnpdf"
  fi

  ## Try to do better, using the nnpdf-config script
  AC_PATH_PROG(NNPDFCONFIG, nnpdf-config, [], [$NNPDFCONFIG_PATH])
  if test -x "$NNPDFCONFIG"; then
    AC_MSG_NOTICE(Using $NNPDFCONFIG to find NNPDF flags)
    NNPDFPATH=`$NNPDFCONFIG --prefix`
    NNPDFINCPATH="$NNPDFCONFIG --incdir"
    NNPDF_CPPFLAGS=`$NNPDFCONFIG --cppflags`
    NNPDF_CXXFLAGS=`$NNPDFCONFIG --cppflags`
    NNPDF_LDFLAGS=`$NNPDFCONFIG --ldflags`
  fi

  ## If it's worked, propagate the variables and execute success arg
  if test -e "$NNPDFPATH"; then
    ## Otherwise  execute the fail arg
    AC_SUBST([NNPDFPATH])
    AC_SUBST([NNPDFINCPATH])
    AC_SUBST([NNPDF_CPPFLAGS])
    AC_SUBST([NNPDF_CXXFLAGS])
    AC_SUBST([NNPDF_LDFLAGS])
    AM_CONDITIONAL([WITH_NNPDF], true)
    AM_CONDITIONAL([WITH_NNPDFLIB], true)
    AM_CONDITIONAL([WITH_NNPDFINC], true)
    AM_CONDITIONAL([WITHOUT_NNPDF], false)
    AM_CONDITIONAL([WITHOUT_NNPDFLIB], false)
    AM_CONDITIONAL([WITHOUT_NNPDFINC], false)
    AC_MSG_NOTICE([NNPDF include path is $NNPDFINCPATH])
    AC_MSG_NOTICE([NNPDF CPPFLAGS is $NNPDF_CPPFLAGS])
    AC_MSG_NOTICE([NNPDF CXXFLAGS is $NNPDF_CXXFLAGS])
    AC_MSG_NOTICE([NNPDF LDFLAGS is $NNPDF_LDFLAGS])
    $1
  else
    AM_CONDITIONAL([WITH_NNPDF], false)
    AM_CONDITIONAL([WITH_NNPDFLIB], false)
    AM_CONDITIONAL([WITH_NNPDFINC], false)
    AM_CONDITIONAL([WITHOUT_NNPDF], true)
    AM_CONDITIONAL([WITHOUT_NNPDFLIB], true)
    AM_CONDITIONAL([WITHOUT_NNPDFINC], true)
    $2
  fi
])
