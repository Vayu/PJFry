#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.65])

# Package version number
# Major.Minor goes into symbols version
m4_define([pjfry_major_version], [1])
m4_define([pjfry_minor_version], [0])
m4_define([pjfry_micro_version], [0-a2])

# Libtool interface number ( current : revision : age )
# - If the library source code has changed at all since the last update,
#   then increment revision
# - If any interfaces have been added, removed, or changed since the last update,
#   increment current, and set revision to 0.
# - If any interfaces have been added since the last public release,
#   then increment age.
# - If any interfaces have been removed or changed since the last public release,
#   then set age to 0.
m4_define([pjfry_interface_current],  [0])
m4_define([pjfry_interface_revision], [0])
m4_define([pjfry_interface_age],      [0])

m4_define([pjfry_version], m4_defn([pjfry_major_version]).m4_defn([pjfry_minor_version]).m4_defn([pjfry_micro_version]))

AC_INIT([pjfry], m4_defn([pjfry_version]))

PJFRY_MAJOR=pjfry_major_version
PJFRY_MINOR=pjfry_minor_version
PJFRY_MICRO=pjfry_micro_version

PJFRY_CURRENT=pjfry_interface_current
PJFRY_REVISION=pjfry_interface_revision
PJFRY_AGE=pjfry_interface_age

PJFRY_VERSION=$PJFRY_MAJOR.$PJFRY_MINOR.$PJFRY_MICRO
PJFRY_VERSION_INFO=$PJFRY_CURRENT:$PJFRY_REVISION:$PJFRY_AGE

AC_CONFIG_SRCDIR([src/common.h])

AC_CONFIG_AUX_DIR([auxfiles])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wall -Werror silent-rules foreign])

AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_F77
AC_PROG_FC
AC_F77_LIBRARY_LDFLAGS
AC_F77_WRAPPERS

# AC_PROG_RANLIB
AC_PROG_LIBTOOL

AC_PROG_LD
AC_PROG_CPP
AC_PROG_INSTALL
#AC_CHECK_TOOL(SED, sed, :)
AC_PROG_LN_S

# Set language
AC_LANG([C++])

# Checks for header files.
AC_CHECK_HEADERS([inttypes.h])
AC_HEADER_STDC

# Checks for typedefs, structures, and compiler characteristics.
#AC_C_INLINE
#AC_TYPE_UINT64_T

# Checks for library functions.
#AC_CHECK_FUNCS([sqrt])

AC_MSG_CHECKING([if libraries can be versioned])
# Special case for PE/COFF platforms: ld reports
# support for version-script, but doesn't actually
# do anything with it.
case $host in
*cygwin* | *mingw32* )
    have_ld_version_script=no
    AC_MSG_RESULT(no)
;;
* )
GLD=`$LD --help < /dev/null 2>/dev/null | grep version-script`
if test "$GLD"; then
    have_ld_version_script=yes
    AC_MSG_RESULT(yes)
else
    have_ld_version_script=no
    AC_MSG_RESULT(no)
    AC_MSG_WARN(*** You have not enabled versioned symbols.)
fi
;;
esac

AM_CONDITIONAL(HAVE_LD_VERSION_SCRIPT, test "$have_ld_version_script" = "yes")

# ------------------------- # ------------------------- # -------------------------
AC_ARG_WITH([integrals],
            [AS_HELP_STRING([--with-integrals],
            [library for scalar loop integrals @<:@default=auto@:>@.
@<:@qcdloop1@:>@ tries to use -lqcdloop1.
@<:@oneloop@:>@ tries to use -lavh_olo.
@<:@qcdloop@:>@ use libqcdloop.a and libff.a, CXXFLAGS and LDFLAGS need to be set manually.
])],
            [test "x$with_integrals" == "xyes" && with_integrals=auto],
            [with_integrals=auto])

case $with_integrals in
auto)
  AC_CHECK_LIB([qcdloop1], [qlinit_],[LIBS="$LIBS -lqcdloop1" use_qcdloop1=yes])
  AC_CHECK_HEADERS([qcdloop1.h])
  if test "x$use_qcdloop1" == "xyes" ; then
    AC_DEFINE(USE_QCDLOOP,"1",[use libqcdloop for scalar integrals])
  else
    AC_F77_FUNC(qlinit)
    AC_CHECK_LIB([qcdloop], [$qlinit],[LIBS="$LIBS -lqcdloop -lff" use_qcdloop=yes],[],[-L$FLIBS])
    if test "x$use_qcdloop" == "xyes" ; then
      AC_DEFINE(USE_QCDLOOP,"1",[use libqcdloop for scalar integrals])
    else
      AC_F77_FUNC(avh_olo_mu_set)
      AC_CHECK_LIB([avh_olo], [$avh_olo_mu_set],[LIBS="$LIBS -lavh_olo" use_oneloop=yes],[],[-L$FLIBS])
      if test "x$use_oneloop" == "xyes" ; then
        AC_DEFINE(USE_ONELOOP,"1",[use libavh_olo for scalar integrals])
      else
        AC_MSG_ERROR([can't find usable scalar integrals library])
      fi
    fi
  fi
;;
qcdloop1)
  AC_CHECK_LIB([qcdloop1], [qlinit_],
    [LIBS="$LIBS -lqcdloop1" use_qcdloop1=yes AC_DEFINE(USE_QCDLOOP,"1",[use libqcdloop for scalar integrals])],
    [AC_MSG_ERROR([libqcdloop1 library is missing])])
  AC_CHECK_HEADERS([qcdloop1.h])
;;
oneloop)
  AC_F77_FUNC(avh_olo_mu_set)
  AC_CHECK_LIB([avh_olo], [$avh_olo_mu_set],
    [LIBS="$LIBS -lavh_olo" use_oneloop=yes AC_DEFINE(USE_ONELOOP,"1",[use libavh_olo for scalar integrals])],
    [AC_MSG_ERROR([libavh_olo library is missing])],[-L$FLIBS])
;;
qcdloop)
  AC_F77_FUNC(qlinit)
  AC_CHECK_LIB([qcdloop], [$qlinit],
    [LIBS="$LIBS -lqcdloop -lff" use_qcdloop=yes AC_DEFINE(USE_QCDLOOP,"1",[use libqcdloop for scalar integrals])],
    [AC_MSG_ERROR([libqcdloop library is missing])],[-L$FLIBS])
;;
*)
  AC_MSG_ERROR([Bad value for --with-integrals: $with_integrals])
;;
esac
# ------------------------- # ------------------------- # -------------------------

# ------------------------- # ------------------------- # -------------------------
AC_ARG_WITH([mcc-path],
            [AS_HELP_STRING([--with-mcc-path],
            [build Mathematica MathLink interface @<:@default=auto@:>@. Set this to full path to mcc tool])],
            [test "x$with_mcc_path" == "xyes" && with_mcc_path=auto],
            [with_mcc_path=auto])
# ------------------------- # -------------------------
if test "x$with_mcc_path" == "xauto" ; then
AC_MSG_CHECKING([for mcc in system path])
AC_PATH_PROGS_FEATURE_CHECK([MCC], [mcc mcc7.0 mcc6.0],
[[
  cat - <<_ACEOF > conftest.cpp
    #include <mathlink.h>
    int main () { return 0; }
_ACEOF
  rm -f conftest$ac_exeext
  $ac_path_MCC conftest.cpp -o conftest$ac_exeext > conftest.err
  if test -s conftest$ac_exeext; then
    grep -v '^ *+' conftest.err >conftest.er1
    cat conftest.er1 >&5
    mv -f conftest.er1 conftest.err
  fi
  test ! -s conftest.err && test -s conftest$ac_exeext \
    && ac_cv_path_MCC=$ac_path_MCC ac_path_MCC_found=:
  rm -f conftest.cpp conftest$ac_exeext conftest.err
]], [with_mcc_path=])
if test "x$ac_cv_path_MCC" != "x" ; then
  AC_SUBST([MCC],[$ac_cv_path_MCC])
  AC_MSG_RESULT($ac_cv_path_MCC)
else
  AC_MSG_RESULT(no)
fi
if test `uname -s` == "Linux" ; then
  case `uname -m` in
    ia64)   SystemIDList="Linux-IA64";;
    i?86)   SystemIDList="Linux";;
    x86_64) SystemIDList="Linux-x86-64";;
    *)      SystemIDList="Unknown";;
  esac
  math_p1a=/opt/products/mathematica
  math_p1b=/usr/local/Wolfram/Mathematica
  math_p2=SystemFiles/Links/MathLink/DeveloperKit/$SystemIDList/CompilerAdditions
  math_ps=$math_p1a/7.0/$math_p2$PATH_SEPARATOR$math_p1b/7.0/$math_p2$PATH_SEPARATOR$math_p1a/6.0/$math_p2$PATH_SEPARATOR$math_p1b/6.0/$math_p2
  AC_PATH_PROGS(MPREP,[mprep],[],[$math_ps])
  if test "x$ac_cv_path_MCC" == "x" ; then
    AC_PATH_PROGS(MCC,[mcc mcc7.0 mcc6.0],[],[$math_ps])
  fi
fi
# -------------------------
else
  if test "x$with_mcc_path" != "xno" ; then
    if test -d $with_mcc_path ; then
      ac_cv_path_MCC=$with_mcc_path/mcc
      ac_cv_path_MPREP=$with_mcc_path/mprep
    else
      ac_cv_path_MCC=$with_mcc_path
    fi
    if test "x$ac_cv_path_MCC" != "x" && test -x $ac_cv_path_MCC ; then
      AC_SUBST([MCC],[$ac_cv_path_MCC])
    fi
    if test "x$ac_cv_path_MPREP" != "x" && test -x $ac_cv_path_MPREP ; then
      AC_SUBST([MPREP],[$ac_cv_path_MPREP])
    fi
  fi
fi
# ------------------------- # ------------------------- # -------------------------
AM_CONDITIONAL(HAVE_MCC,   test "x$with_mcc_path" != "xno" && test "x$ac_cv_path_MCC"   != "x")
AM_CONDITIONAL(HAVE_MPREP, test "x$with_mcc_path" != "xno" && test "x$ac_cv_path_MPREP" != "x")
AM_CONDITIONAL(HAVE_FLIBS, test "x$use_qcdloop" == "xyes" || test "x$use_oneloop" == "xyes")


AC_ARG_ENABLE(golem-mode,[AS_HELP_STRING([--enable-golem-mode],
            [build golem95-like interface @<:@default=no@:>@])],
            [test "x$enable_golem_mode" == "xyes" && AC_DEFINE(USE_GOLEM_MODE,"1","build golem95-like interface")],
            [enable_golem_mode=no])
AM_CONDITIONAL(EMU_GOLEM, test "x$enable_golem_mode" == "xyes")


AC_ARG_ENABLE(cache-high,[AS_HELP_STRING([--enable-cache-high],
            [enable caching of 5- and 6-point coefficients @<:@default=yes@:>@])],
            [],
            [enable_cache_high=yes])
if test "x$enable_cache_high" == "xyes" ; then
  AC_DEFINE(USE_CACHE_HIGH,"1","enable cache for 5 and 6-point coefficients")
fi

AC_ARG_ENABLE(cache-low,[AS_HELP_STRING([--enable-cache-low],
            [enable caching of 2-,3-,4-point coefficients @<:@default=no@:>@])],
            [test "x$enable_cache_low" == "xyes" && AC_DEFINE(USE_CACHE_LOW,"1","enable cache for 2-,3-,4-point coefficients")],
            [enable_cache_low=no])


# Substitutions for .in files
AC_SUBST(PJFRY_MAJOR)
AC_SUBST(PJFRY_MINOR)
AC_SUBST(PJFRY_MICRO)
AC_SUBST(PJFRY_CURRENT)
AC_SUBST(PJFRY_REVISION)
AC_SUBST(PJFRY_AGE)
AC_SUBST(PJFRY_VERSION_INFO)

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 mlink/Makefile
                 examples/Makefile])

AC_OUTPUT
