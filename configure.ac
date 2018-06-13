dnl Copyright (C) 2009, 2010, 2011, 2012
dnl Computer Graphics Group, University of Siegen, Germany.
dnl Written by Martin Lambers <martin.lambers@uni-siegen.de>.
dnl See http://www.cg.informatik.uni-siegen.de/ for contact information.
dnl
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program.  If not, see <http://www.gnu.org/licenses/>.

dnl Autotools init stuff
AC_PREREQ([2.65])
AC_INIT([libecmdb], [0.7.1], [martin.lambers@uni-siegen.de], [libecmdb], [https://marlam.de/ecm])
AC_CONFIG_SRCDIR([src/ecmdb/ecmdb/ecm.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AC_CANONICAL_TARGET
AC_USE_SYSTEM_EXTENSIONS
AC_SYS_LARGEFILE
AM_INIT_AUTOMAKE([1.11.1 color-tests silent-rules tar-pax no-texinfo.tex -Wall -Wno-portability])
AM_SILENT_RULES([yes])
LT_PREREQ([2.2.6])
LT_INIT([win32-dll])
AC_PROG_CXX
AC_PROG_INSTALL
AC_LANG([C++])

dnl We want C++0x, and blindly add -std=gnu++0x. This only works for gcc (and possibly clang)!
CXXFLAGS="-std=gnu++0x $CXXFLAGS"

dnl Configure switches
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug], [Enable debugging checks. Disabled by default.])],
    [if test "$enableval" = "yes"; then ndebug=0; else ndebug=1; fi], [ndebug=1])
if test "$ndebug" = "1"; then
    AC_DEFINE_UNQUOTED([NDEBUG], [1], [Define to 1 to disable assertions and other debugging checks.])
fi
AC_ARG_ENABLE([la-file],
    [AS_HELP_STRING([--enable-la-file], [Install the .la file. Disabled by default when a shared library is built.])],
    [if test "$enableval" = "yes"; then install_la_file="yes"; else install_la_file="no"; fi],
    [if test "$enable_shared" = "yes"; then install_la_file="no"; else install_la_file="yes"; fi])
AM_CONDITIONAL([INSTALL_LA_FILE], [test "$install_la_file" = "yes"])
AC_ARG_ENABLE([reference],
    [AS_HELP_STRING([--enable-reference], [Build the reference. Enabled by default if doxygen is available.])],
    [if test "$enableval" = "yes"; then build_ref="yes"; else build_ref="no"; fi], [build_ref="yes"])
AC_ARG_VAR([DOXYGEN], [Doxygen command])
if test "$build_ref" = "yes"; then
    AC_CHECK_PROGS([DOXYGEN], [doxygen], [:])
    if test "$DOXYGEN" = ":"; then
        build_ref="no"
        AC_MSG_WARN([Disabled building the reference because doxygen was not found.])
    fi
fi
AM_CONDITIONAL([BUILD_REF], [test "$build_ref" = "yes"])

dnl Feature checks needed by the base modules
AC_LANG_PUSH([C])
dnl - dbg
AC_CHECK_FUNCS([backtrace sigaction])
dnl - fio
case "${target}" in *-*-mingw*) LIBS="$LIBS -lshlwapi" ;; esac
AC_CHECK_FUNCS([fdatasync fnmatch fseeko ftello getpwuid link mmap posix_fadvise readdir_r symlink])
dnl - str
AM_ICONV([])
AC_CHECK_FUNCS([nl_langinfo vasprintf wcswidth])
AC_LANG_POP([C])

dnl Feature checks needed by ecm.cpp
AC_CHECK_FUNCS([sincos])

dnl pkg-config (required to detect libraries)
PKG_PROG_PKG_CONFIG([])
if test -z "$PKG_CONFIG"; then
    AC_MSG_ERROR([pkg-config not found])
fi

dnl libgta
PKG_CHECK_MODULES([libgta], [gta >= 0.0], [HAVE_LIBGTA=1], [HAVE_LIBGTA=0])
if test "$HAVE_LIBGTA" != "1"; then
    AC_MSG_ERROR([libgta not found: $libgta_PKG_ERRORS])
fi

dnl libjpeg
AC_LIB_HAVE_LINKFLAGS([jpeg], [],
    [#include <stdio.h>
     #include <jpeglib.h>],
    [jpeg_create_compress(0);])
if test "$HAVE_LIBJPEG" != "yes"; then
    AC_MSG_ERROR([libjpeg not found])
fi

dnl libecmdb package version
AC_SUBST([ECMDB_VERSION], [$PACKAGE_VERSION])
AC_SUBST([ECMDB_VERSION_MAJOR], [`echo $ECMDB_VERSION | sed -e 's/\(.*\)\..*\..*/\1/'`])
AC_SUBST([ECMDB_VERSION_MINOR], [`echo $ECMDB_VERSION | sed -e 's/.*\.\(.*\)\..*/\1/'`])
AC_SUBST([ECMDB_VERSION_PATCH], [`echo $ECMDB_VERSION | sed -e 's/.*\..*\.\(.*\)/\1/' -e 's/[[^0-9]].*//'`])
AC_SUBST([ECMDB_VERSION_NUMBER], [`printf "0x%02x%02x%02x" $ECMDB_VERSION_MAJOR $ECMDB_VERSION_MINOR $ECMDB_VERSION_PATCH`])

dnl libecmdb shared library version
dnl See http://www.gnu.org/software/libtool/manual/libtool.html#Updating-version-info
dnl for a description of library versioning. In summary:
dnl Library code modified:            REVISION++
dnl Interfaces changed/added/removed: CURRENT++ REVISION=0
dnl Interfaces added:                 AGE++
dnl Interfaces removed:               AGE=0
AC_SUBST([LT_CURRENT], [1])
AC_SUBST([LT_REVISION], [0])
AC_SUBST([LT_AGE], [0])

dnl Output
AC_CONFIG_FILES([Makefile \
	doc/Makefile \
	doc/doxyfile \
	src/base/Makefile \
	src/glvm/Makefile \
	src/ecmdb/Makefile \
	src/ecmdb/ecmdb/version.h \
	src/ecmdb/ecmdb.pc \
	src/Makefile \
	tests/Makefile \
	])
AC_OUTPUT
