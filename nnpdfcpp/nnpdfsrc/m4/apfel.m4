#AC_SEARCH_LHAPDF(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_APFEL],[
  AC_ARG_WITH([apfel], AC_HELP_STRING(--with-apfel, [path to APFEL library and header files]))

  ## Use a specified --with-lhapdf arg to set basic paths, if provided
  APFELCONFIG_PATH=$PATH
  if test -e "$with_apfel"; then
    APFELCONFIG_PATH=$with_apfel/bin:$APFELCONFIG_PATH
    APFELPATH="$with_apfel"
    APFELINCPATH="$APFELPATH/include"
    APFEL_CPPFLAGS="-I$APFELINCPATH"
    APFEL_CXXFLAGS=""
    APFEL_LDFLAGS="-L$APFELPATH/lib -lAPFEL"
  fi

  ## Try to do better, using the lhapdf-config script
  AC_PATH_PROG(APFELCONFIG, apfel-config, [], [$APFELCONFIG_PATH])
  if test -x "$APFELCONFIG"; then
    AC_MSG_NOTICE(Using $APFELCONFIG to find APFEL flags)
    APFELPATH=`$APFELCONFIG --prefix`
    APFELINCPATH="$APFELPATH/include"
    APFEL_CPPFLAGS=`$APFELCONFIG --cppflags`
    APFEL_CXXFLAGS=`$APFELCONFIG --cppflags`
    APFEL_LDFLAGS=`$APFELCONFIG --ldflags`
  fi

  ## If it's worked, propagate the variables and execute success arg
  if test -e "$APFELPATH"; then
    ## Otherwise  execute the fail arg
    AC_SUBST([APFELPATH])
    AC_SUBST([APFELINCPATH])
    AC_SUBST([APFEL_CPPFLAGS])
    AC_SUBST([APFEL_CXXFLAGS])
    AC_SUBST([APFEL_LDFLAGS])
    AM_CONDITIONAL([WITH_APFEL], true)
    AM_CONDITIONAL([WITH_APFELLIB], true)
    AM_CONDITIONAL([WITH_APFELINC], true)
    AM_CONDITIONAL([WITHOUT_APFEL], false)
    AM_CONDITIONAL([WITHOUT_APFELLIB], false)
    AM_CONDITIONAL([WITHOUT_APFELINC], false)
    AC_MSG_NOTICE([APFEL include path is $APFELINCPATH])
    AC_MSG_NOTICE([APFEL CPPFLAGS is $APFEL_CPPFLAGS])
    AC_MSG_NOTICE([APFEL CXXFLAGS is $APFEL_CXXFLAGS])
    AC_MSG_NOTICE([APFEL LDFLAGS is $APFEL_LDFLAGS])
    $1
  else
    AM_CONDITIONAL([WITH_APFEL], false)
    AM_CONDITIONAL([WITH_APFELLIB], false)
    AM_CONDITIONAL([WITH_APFELINC], false)
    AM_CONDITIONAL([WITHOUT_APFEL], true)
    AM_CONDITIONAL([WITHOUT_APFELLIB], true)
    AM_CONDITIONAL([WITHOUT_APFELINC], true)
    AC_MSG_ERROR('APFEL not installed!')
    $2
  fi
])
