# Autoconf file for the Fixmath library.

#  Copyright (C) 2005-2011, Axis Communications AB, LUND, SWEDEN
#
#  This file is part of Fixmath.
#
#  Fixmath is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published
#  by the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  You can use the comments under either the terms of the GNU Lesser General
#  Public License version 3 as published by the Free Software Foundation,
#  either version 3 of the License or (at your option) any later version, or
#  the GNU Free Documentation License version 1.3 or any later version
#  published by the Free Software Foundation; with no Invariant Sections, no
#  Front-Cover Texts, and no Back-Cover Texts.
#  A copy of the license is included in the documentation section entitled
#  "GNU Free Documentation License".
#
#  Fixmath is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License and a copy of the GNU Free Documentation License along
#  with Fixmath. If not, see <http://www.gnu.org/licenses/>.

# Initialize
AC_INIT([fixmath], [1.4])
AC_CONFIG_SRCDIR([include/fixmath.h])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# -----------------------------------------------------------

# VERSION NUMBERS
# ===============
# Sets the compile-time and run-time major and minor version numbers.
# Update values accordingly when making a new release.

# Major version number
FX_VERSION=1

# Minor version number
FX_MINOR=4

# LIBTOOL VERSIONING
# ==================
# Controls the libtool version numbers CURRENT.REVISION.AGE.
# Do not try to synchronize these version numbers with the
# major and minor version numbers above - they serve two different
# purposes. For more information, consult the libtool manual
# http://www.gnu.org/software/libtool/manual.

# Current interface version
FX_LIB_CUR=1

# Interface revision
FX_LIB_REV=1

# Interface age
FX_LIB_AGE=1

# -----------------------------------------------------------

# Check version number consistency
if test "${FX_VERSION}.${FX_MINOR}" != ${VERSION}; then
    AC_MSG_ERROR([Conflicting version numbers set])
fi

# Save the user CFLAGS before AC_PROG_C89 modifies it
USER_CFLAGS=${CFLAGS}

# Check for programs
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_LIBTOOL
AC_PROG_CC_C89

# Restore CFLAGS
CFLAGS=${USER_CFLAGS}

# Check for headers
AC_HEADER_STDC
AC_CHECK_HEADER([stdint.h],,AC_MSG_ERROR([Required header stdint.h missing]))

# Check the C compiler vendor
AX_COMPILER_VENDOR

# Check for the CHECK unit test framework
m4_pattern_allow(PKG_CONFIG_LIBDIR)
PKG_CONFIG_LIBDIR_SAVE=$PKG_CONFIG_LIBDIR
unset PKG_CONFIG_LIBDIR
PKG_CHECK_MODULES([CHECK],
                  [check >= 0.9.4],
                  [have_check=yes],
                  [have_check=no])
AM_CONDITIONAL(HAVE_CHECK, test x"$have_check" = "xyes")

# Check for Doxygen
AC_CHECK_PROGS(DOXYGEN, doxygen)

# Check for doxygen option
AC_ARG_WITH(doxygen,
  AC_HELP_STRING([--with-doxygen], [Generate Doxygen documentation]),
  [case "${withval}" in
   yes)
       if test -n "${DOXYGEN}"; then
           OPT_DOXYGEN="yes"
       fi
       ;;
   no)
       OPT_DOXYGEN=
       ;;
   *)
       AC_MSG_ERROR(Bad value ${withval} for --with-doxygen option)
       ;;
   esac])
AM_CONDITIONAL(HAVE_DOXYGEN, test -n "${OPT_DOXYGEN}")

# Check for debug option
OPT_DEBUG=default
asserted_returns=0
AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([[--enable-debug@<:@=no@:>@]], [Enable debugging:])
   AS_HELP_STRING(,[no:               Disable debugging])
   AS_HELP_STRING(,[yes:              Enable debugging])
   AS_HELP_STRING(,[asserted-returns: For erroneous parameters,]
                   [ call abort() instead of returning])],
  [OPT_DEBUG=${enableval}
   case ${enableval} in
   yes|no) ;;
   asserted-returns)
       OPT_DEBUG=yes
       asserted_returns=1
       ;;
   *)
       AC_MSG_ERROR([[Bad value "${enableval}" for --enable-debug option]])
       ;;
   esac])
AC_DEFINE_UNQUOTED([RC_ASSERTED_RETURNS], [${asserted_returns}],
                   [Define to 1 to not return but call abort()]
                   [ for erroneous parameters])

# Add the BUILD flag
AM_CPPFLAGS="${AM_CPPFLAGS} -DFX_BUILD"

# Set the default compiler flags
fx_gflag=
case ${ax_cv_c_compiler_vendor} in
gnu)
    AM_CFLAGS="${AM_CFLAGS} -std=c89 -W -Wall"
    AX_CFLAGS_GCC_OPTION([-fvisibility=hidden], [AM_CFLAGS])
    AX_CFLAGS_GCC_OPTION([-g3], [fx_gflag])
    ;;
esac

# Set debug/release compiler flags
if test "${OPT_DEBUG}" = "yes"; then
    AM_CFLAGS="${AM_CFLAGS} -O0 -g"
else
    AM_CPPFLAGS="${AM_CPPFLAGS} -DNDEBUG=1"
    AM_CFLAGS="${AM_CFLAGS} -O2"
    case ${ax_cv_c_compiler_vendor} in
    gnu)
        AM_CFLAGS="${AM_CFLAGS} -funroll-loops -fomit-frame-pointer"
        ;;
    esac
fi

test "x${fx_gflag}" = x && test "$ac_cv_prog_cc_g" != no && fx_gflag=-g
test ${OPT_DEBUG} != no && AM_CFLAGS="${AM_CFLAGS} ${fx_gflag}"

# Substitute variables in output files
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([FX_VERSION])
AC_SUBST([FX_MINOR])
AC_SUBST([FX_LIB_CUR])
AC_SUBST([FX_LIB_REV])
AC_SUBST([FX_LIB_AGE])

# Output
AM_CONFIG_HEADER(config.h)
AC_CONFIG_FILES([Makefile
                 Doxyfile
                 fixmath.pc
                 include/fixmath_ver.h
                 src/Makefile
                 prof/Makefile
                 checktests/Makefile])
AC_OUTPUT
