# Copyright (C) 2010, 2011, 2012
# Martin Lambers <martin.lambers@uni-siegen.de>
#
# Copying and distribution of this file, with or without modification, are
# permitted in any medium without royalty provided the copyright notice and this
# notice are preserved. This file is offered as-is, without any warranty.

EXTRA_DIST = doxyfile.in

if BUILD_REF

BUILT_SOURCES = reference/index.html
nobase_nodist_doc_DATA = $(wildcard reference/*)
CLEANFILES = $(wildcard reference/*)

reference/index.html: doxyfile $(top_srcdir)/src/ecmdb/ecmdb/ecm.h $(top_srcdir)/src/ecmdb/ecmdb/ecmdb.h
	$(DOXYGEN) doxyfile

endif
