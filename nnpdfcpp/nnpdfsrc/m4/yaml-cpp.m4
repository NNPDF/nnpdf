#AC_SEARCH_GSL(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_YAML],[
  AC_ARG_WITH([yaml-cpp], AC_HELP_STRING(--with-yaml-cpp, [path to yaml-cpp library and header files]))

  PKG_CHECK_MODULES(YAML, yaml-cpp)

  ## Use a specified --with-sqlite3 arg to set basic paths, if provided
  YAMLCONFIG_PATH=$PATH
  YAMLINCPATH="$YAMLPATH/include"
  YAML_CPPFLAGS=`pkg-config --cflags yaml-cpp`
  YAML_CXXFLAGS=`pkg-config --cflags yaml-cpp`
  YAML_LDFLAGS=`pkg-config --libs yaml-cpp`

  ## If it's worked, propagate the variables and execute success arg
  ## Otherwise  execute the fail arg
    AC_SUBST([YAMLPATH])
    AC_SUBST([YAMLINCPATH])
    AC_SUBST([YAML_CPPFLAGS])
    AC_SUBST([YAML_CXXFLAGS])
    AC_SUBST([YAML_LDFLAGS])
    AM_CONDITIONAL([WITH_YAML], true)
    AM_CONDITIONAL([WITH_YAMLLIB], true)
    AM_CONDITIONAL([WITH_YAMLINC], true)
    AM_CONDITIONAL([WITHOUT_YAML], false)
    AM_CONDITIONAL([WITHOUT_YAMLLIB], false)
    AM_CONDITIONAL([WITHOUT_YAMLINC], false)
    AC_MSG_NOTICE([YAML include path is $YAMLINCPATH])
    AC_MSG_NOTICE([YAML CPPFLAGS is $YAML_CPPFLAGS])
    AC_MSG_NOTICE([YAML CXXFLAGS is $YAML_CXXFLAGS])
    AC_MSG_NOTICE([YAML LDFLAGS is $YAML_LDFLAGS])
    $1
])
