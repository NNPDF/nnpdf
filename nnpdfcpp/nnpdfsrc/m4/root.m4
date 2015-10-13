#AC_SEARCH_ROOT(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_ROOT],[ AC_ARG_WITH([root],
AC_HELP_STRING(--with-root, [path to ROOT library and header files]))

  ## Use a specified --with-root arg to set basic paths, if provided
  ROOTCONFIG_PATH=$PATH
  if test -e "$with_root"; then
    ROOTCONFIG_PATH=$with_root/bin:$ROOTCONFIG_PATH
    ROOTPATH="$with_root"
    ROOTINCPATH="$ROOTPATH/include"
    ROOT_CPPFLAGS="-I$ROOTINCPATH"
    ROOT_CXXFLAGS=""
    ROOT_LDFLAGS="-L$ROOTPATH/lib -lroot -lm"
  fi

  ## Try to do better, using the root-config script
  AC_PATH_PROG(ROOTCONFIG, root-config, [], [$ROOTCONFIG_PATH])
  if test -x "$ROOTCONFIG"; then
    AC_MSG_NOTICE(Using $ROOTCONFIG to find ROOT flags)
    ROOTPATH=`$ROOTCONFIG --prefix`
    ROOTINCPATH="$ROOTPATH/include"
    ROOT_CPPFLAGS=`$ROOTCONFIG --cflags`
    ROOT_CXXFLAGS=`$ROOTCONFIG --cflags`
    ROOT_LDFLAGS=`$ROOTCONFIG --glibs`
  fi

  ## If it's worked, propagate the variables and execute success arg
  if test -e "$ROOTPATH"; then
    ## Otherwise  execute the fail arg
    AC_SUBST([ROOTPATH])
    AC_SUBST([ROOTINCPATH])
    AC_SUBST([ROOT_CPPFLAGS])
    AC_SUBST([ROOT_CXXFLAGS])
    AC_SUBST([ROOT_LDFLAGS])
    AM_CONDITIONAL([WITH_ROOT], true)
    AM_CONDITIONAL([WITH_ROOTLIB], true)
    AM_CONDITIONAL([WITH_ROOTINC], true)
    AM_CONDITIONAL([WITHOUT_ROOT], false)
    AM_CONDITIONAL([WITHOUT_ROOTLIB], false)
    AM_CONDITIONAL([WITHOUT_ROOTINC], false)
    AC_MSG_NOTICE([ROOT include path is $ROOTINCPATH])
    AC_MSG_NOTICE([ROOT CPPFLAGS is $ROOT_CPPFLAGS])
    AC_MSG_NOTICE([ROOT CXXFLAGS is $ROOT_CXXFLAGS])
    AC_MSG_NOTICE([ROOT LDFLAGS is $ROOT_LDFLAGS])
    $1
  else
    AM_CONDITIONAL([WITH_ROOT], false)
    AM_CONDITIONAL([WITH_ROOTLIB], false)
    AM_CONDITIONAL([WITH_ROOTINC], false)
    AM_CONDITIONAL([WITHOUT_ROOT], true)
    AM_CONDITIONAL([WITHOUT_ROOTLIB], true)
    AM_CONDITIONAL([WITHOUT_ROOTINC], true)
    AC_MSG_ERROR('ROOT not installed!')
    $2
  fi
])
