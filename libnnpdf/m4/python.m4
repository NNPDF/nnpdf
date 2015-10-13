#AC_SEARCH_PYTHON(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_PYTHON],[
  AC_ARG_WITH([python], AC_HELP_STRING(--with-python, [path to PYTHON library and header files]))

  ## Use a specified --with-python arg to set basic paths, if provided
  PYTHONCONFIG_PATH=$PATH
  if test -e "$with_python"; then
    PYTHONCONFIG_PATH=$with_python/bin:$PYTHONCONFIG_PATH
    PYTHONPATH="$with_python"
    PYTHONINCPATH="$PYTHONPATH/include"
    PYTHON_CPPFLAGS="-I$PYTHONINCPATH"
    PYTHON_CXXFLAGS=""
    PYTHON_LDFLAGS="-L$PYTHONPATH/lib -lpython -lm"
  fi

  ## Try to do better, using the python-config script
  AC_PATH_PROG(PYTHONCONFIG, python-config, [], [$PYTHONCONFIG_PATH])
  if test -x "$PYTHONCONFIG"; then
    AC_MSG_NOTICE(Using $PYTHONCONFIG to find PYTHON flags)
    PYTHONPATH=`$PYTHONCONFIG --prefix`
    PYTHONINCPATH="$PYTHONPATH/include"
    PYTHON_CPPFLAGS=`$PYTHONCONFIG --cflags`
    PYTHON_CXXFLAGS=`$PYTHONCONFIG --cflags`
    PYTHON_LDFLAGS=`$PYTHONCONFIG --libs`
  fi

  ## If it's worked, propagate the variables and execute success arg
  if test -e "$PYTHONPATH"; then
    ## Otherwise  execute the fail arg
    AC_SUBST([PYTHONPATH])
    AC_SUBST([PYTHONINCPATH])
    AC_SUBST([PYTHON_CPPFLAGS])
    AC_SUBST([PYTHON_CXXFLAGS])
    AC_SUBST([PYTHON_LDFLAGS])
    AM_CONDITIONAL([WITH_PYTHON], true)
    AM_CONDITIONAL([WITH_PYTHONLIB], true)
    AM_CONDITIONAL([WITH_PYTHONINC], true)
    AM_CONDITIONAL([WITHOUT_PYTHON], false)
    AM_CONDITIONAL([WITHOUT_PYTHONLIB], false)
    AM_CONDITIONAL([WITHOUT_PYTHONINC], false)
    AC_MSG_NOTICE([PYTHON include path is $PYTHONINCPATH])
    AC_MSG_NOTICE([PYTHON CPPFLAGS is $PYTHON_CPPFLAGS])
    AC_MSG_NOTICE([PYTHON CXXFLAGS is $PYTHON_CXXFLAGS])
    AC_MSG_NOTICE([PYTHON LDFLAGS is $PYTHON_LDFLAGS])
    $1
  else
    AM_CONDITIONAL([WITH_PYTHON], false)
    AM_CONDITIONAL([WITH_PYTHONLIB], false)
    AM_CONDITIONAL([WITH_PYTHONINC], false)
    AM_CONDITIONAL([WITHOUT_PYTHON], true)
    AM_CONDITIONAL([WITHOUT_PYTHONLIB], true)
    AM_CONDITIONAL([WITHOUT_PYTHONINC], true)
    $2
  fi
])
