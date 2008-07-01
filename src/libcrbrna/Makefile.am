# Last modified: 2008-06-30.15

# Copyright (c) 2008 Stefan Bienert <bienert@zbh.uni-hamburg.de>
# Copyright (c) 2008 Center for Bioinformatics, University of Hamburg 

# See COPYING file in the top level directory of this tree for licence.

# This Makefile creates the core library of the CoRB project. This
# library holds standard functions, e.g. for file handling, memory management
# and error handling. Further, it defines some standard data types like arrays,
# stacks and hash tables.

ACLOCAL_AMFLAGS = -I share/misc/m4

AM_CPPFLAGS = -I$(top_srcdir)/src

noinst_LIBRARIES = libcrbrna.a

libcrbrna_a_SOURCES = \
	alphabet.c

noinst_HEADERS =      \
	alphabet.h    \
	crbrna.h

# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:
