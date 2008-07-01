# Last modified: 2008-07-01.09

# Copyright (c) 2008 Stefan Bienert <bienert@zbh.uni-hamburg.de>
# Copyright (c) 2008 Center for Bioinformatics, University of Hamburg 

# See COPYING file in the top level directory of this tree for licence.

# This Makefile creates the applications library. This is, the set of
# applications provided by CoRb. Therefore the library contains the
# main-function pendants of all the tools. See README for more information.


ACLOCAL_AMFLAGS = -I share/misc/m4

AM_CPPFLAGS = -I$(top_srcdir)/src

noinst_LIBRARIES = libcrbapps.a

libcrbapps_a_SOURCES = \
	brot.c         \
	brot_cmdline.c \
	fold.c         \
        fold_cmdline.c

noinst_HEADERS =      \
	brot_cmdline.h\
	brot.h \
	fold_cmdline.h\
        fold.h


# Do we need testing here?

# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:
