# Copyright (C) 2011, 2012
# Martin Lambers <martin.lambers@uni-siegen.de>
#
# Copying and distribution of this file, with or without modification, are
# permitted in any medium without royalty provided the copyright notice and this
# notice are preserved. This file is offered as-is, without any warranty.

check_PROGRAMS = \
	ecm

AM_DEFAULT_SOURCE_EXT = .cpp

TESTS = $(check_PROGRAMS)

AM_CPPFLAGS = \
	-I$(top_srcdir)/src/base \
	-I$(top_srcdir)/src/glvm \
	-I$(top_srcdir)/src/ecm/ecm

LDADD = ../src/ecm/libecmdb.la \
       ../src/base/libbase.la \
	$(LIBICONV)

# Prevent libtool from building annoying wrapper scripts,
# which would prevent us to check with valgrind.
AM_LDFLAGS = -no-install
