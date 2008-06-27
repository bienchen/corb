# Last modified: 2008-06-26.09

# Copyright (c) 2008 Stefan Bienert <bienert@zbh.uni-hamburg.de>
# Copyright (c) 2008 Center for Bioinformatics, University of Hamburg 

# See COPYING file in the top level directory of this tree for licence.

# This Makefile creates the core library of the CoRB project. This
# library holds standard functions, e.g. for file handling, memory management
# and error handling. Further, it defines some standard data types like arrays,
# stacks and hash tables.

# AM_CFLAGS=@CFLAGS@


ACLOCAL_AMFLAGS = -I share/misc/m4

AM_CPPFLAGS = -I$(top_srcdir)/src

noinst_LIBRARIES = libcrbbasic.a # lib_LIBRARIES

# libcrbbasic_adir = $(includedir)/libcrbbasic

libcrbbasic_a_SOURCES =                    \
        memmgr.c                           \
        mprintf.c                          \
        errormsg.c 	                   \
        argvprsr.c	                   \
	str.c

# libcrbbasic_a
noinst_HEADERS =                           \
	memmgr.h                           \
	mprintf.h  	                   \
	mprintf_dummy.h                    \
	errormsg.h 	                   \
	inc_strg.h 	                   \
	inc_pthr.h 	                   \
	inc_bool.h                         \
	gcckywrd.h 	                   \
	argvprsr.h	                   \
	str.h                              \
        undef.h                            \
        crbbasic.h                         \
        genarray.h


# tests for the library
LDADD =                                    \
	libcrbbasic.a                      \
	../libcrbfallback/libcrbfallback.a

check_PROGRAMS =                           \
	test_mprintf                       \
	test_errormsg                      \
	test_memmgr                        \
	test_argvprsr                      \
	test_str                           \
        test_genarray

test_memmgr_SOURCES   = test_memmgr.c

test_mprintf_SOURCES  = test_mprintf.c

test_errormsg_SOURCES = test_errormsg.c

test_argvprsr_SOURCES = test_argvprsr.c

test_str_SOURCES      = test_str.c

test_genarray_SOURCES = test_genarray.c

# add test for correct failing of functions?

TESTS = $(check_PROGRAMS)

# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:
