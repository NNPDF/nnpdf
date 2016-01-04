#AC_SEARCH_GSL(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_SQLITE3],[
  AC_ARG_WITH([sqlite3], AC_HELP_STRING(--with-sqlite3, [path to sqlite3 library and header files]))

  PKG_CHECK_MODULES(SQLITE3, sqlite3)

  ## Use a specified --with-sqlite3 arg to set basic paths, if provided
  SQLITE3CONFIG_PATH=$PATH
  SQLITE3INCPATH="$SQLITE3PATH/include"
  SQLITE3_CPPFLAGS=`pkg-config --cflags sqlite3`
  SQLITE3_CXXFLAGS=`pkg-config --cflags sqlite3`
  SQLITE3_LDFLAGS=`pkg-config --libs sqlite3`

  ## If it's worked, propagate the variables and execute success arg
  ## Otherwise  execute the fail arg
    AC_SUBST([SQLITE3PATH])
    AC_SUBST([SQLITE3INCPATH])
    AC_SUBST([SQLITE3_CPPFLAGS])
    AC_SUBST([SQLITE3_CXXFLAGS])
    AC_SUBST([SQLITE3_LDFLAGS])
    AM_CONDITIONAL([WITH_SQLITE3], true)
    AM_CONDITIONAL([WITH_SQLITE3LIB], true)
    AM_CONDITIONAL([WITH_SQLITE3INC], true)
    AM_CONDITIONAL([WITHOUT_SQLITE3], false)
    AM_CONDITIONAL([WITHOUT_SQLITE3LIB], false)
    AM_CONDITIONAL([WITHOUT_SQLITE3INC], false)
    AC_MSG_NOTICE([SQLITE3 include path is $SQLITE3INCPATH])
    AC_MSG_NOTICE([SQLITE3 CPPFLAGS is $SQLITE3_CPPFLAGS])
    AC_MSG_NOTICE([SQLITE3 CXXFLAGS is $SQLITE3_CXXFLAGS])
    AC_MSG_NOTICE([SQLITE3 LDFLAGS is $SQLITE3_LDFLAGS])
    $1
])
