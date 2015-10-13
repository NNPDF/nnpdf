#AC_SEARCH_LHAPDF(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_LHAPDF],[
  AC_ARG_WITH([lhapdf], AC_HELP_STRING(--with-lhapdf, [path to LHAPDF library and header files]))

  ## Use a specified --with-lhapdf arg to set basic paths, if provided
  LHAPDFCONFIG_PATH=$PATH
  if test -e "$with_lhapdf"; then
    LHAPDFCONFIG_PATH=$with_lhapdf/bin:$LHAPDFCONFIG_PATH
    LHAPDFPATH="$with_lhapdf"
    LHAPDFINCPATH="$LHAPDFPATH/include"
    LHAPDF_CPPFLAGS="-I$LHAPDFINCPATH"
    LHAPDF_CXXFLAGS=""
    LHAPDF_LDFLAGS="-L$LHAPDFPATH/lib -llhapdf -lm"
  fi

  ## Try to do better, using the lhapdf-config script
  AC_PATH_PROG(LHAPDFCONFIG, lhapdf-config, [], [$LHAPDFCONFIG_PATH])
  if test -x "$LHAPDFCONFIG"; then
    AC_MSG_NOTICE(Using $LHAPDFCONFIG to find LHAPDF flags)
    LHAPDFPATH=`$LHAPDFCONFIG --prefix`
    LHAPDFINCPATH="$LHAPDFPATH/include"
    LHAPDF_CPPFLAGS=`$LHAPDFCONFIG --cflags`
    LHAPDF_CXXFLAGS=`$LHAPDFCONFIG --cflags`
    LHAPDF_LDFLAGS=`$LHAPDFCONFIG --libs`
  fi

  ## If it's worked, propagate the variables and execute success arg
  if test -e "$LHAPDFPATH"; then
    ## Otherwise  execute the fail arg
    AC_SUBST([LHAPDFPATH])
    AC_SUBST([LHAPDFINCPATH])
    AC_SUBST([LHAPDF_CPPFLAGS])
    AC_SUBST([LHAPDF_CXXFLAGS])
    AC_SUBST([LHAPDF_LDFLAGS])
    AM_CONDITIONAL([WITH_LHAPDF], true)
    AM_CONDITIONAL([WITH_LHAPDFLIB], true)
    AM_CONDITIONAL([WITH_LHAPDFINC], true)
    AM_CONDITIONAL([WITHOUT_LHAPDF], false)
    AM_CONDITIONAL([WITHOUT_LHAPDFLIB], false)
    AM_CONDITIONAL([WITHOUT_LHAPDFINC], false)
    AC_MSG_NOTICE([LHAPDF include path is $LHAPDFINCPATH])
    AC_MSG_NOTICE([LHAPDF CPPFLAGS is $LHAPDF_CPPFLAGS])
    AC_MSG_NOTICE([LHAPDF CXXFLAGS is $LHAPDF_CXXFLAGS])
    AC_MSG_NOTICE([LHAPDF LDFLAGS is $LHAPDF_LDFLAGS])
    $1
  else
    AM_CONDITIONAL([WITH_LHAPDF], false)
    AM_CONDITIONAL([WITH_LHAPDFLIB], false)
    AM_CONDITIONAL([WITH_LHAPDFINC], false)
    AM_CONDITIONAL([WITHOUT_LHAPDF], true)
    AM_CONDITIONAL([WITHOUT_LHAPDFLIB], true)
    AM_CONDITIONAL([WITHOUT_LHAPDFINC], true)
    AC_MSG_ERROR('LHAPDF not installed!')
    $2
  fi
])
