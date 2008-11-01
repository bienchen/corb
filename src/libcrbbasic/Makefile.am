# Last modified: 2008-11-01.12

#dnl Copyright (C) 2008 Stefan Bienert
#dnl 
#dnl This file is part of CoRB.
#dnl
#dnl CoRB is free software: you can redistribute it and/or modify
#dnl it under the terms of the GNU General Public License as published by
#dnl the Free Software Foundation, either version 3 of the License, or
#dnl (at your option) any later version.
#dnl
#dnl CoRB is distributed in the hope that it will be useful,
#dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
#dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#dnl GNU General Public License for more details.
#dnl 
#dnl You should have received a copy of the GNU General Public License
#dnl along with CoRB.  If not, see <http://www.gnu.org/licenses/>.



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
        gfile.c  	                   \
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
        gfile.h  	                   \
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

#dnl Local variables:
#dnl eval: (add-hook 'write-file-hooks 'time-stamp)
#dnl time-stamp-start: "Last modified: "
#dnl time-stamp-format: "%:y-%02m-%02d.%02H"
#dnl time-stamp-end: "$"
#dnl End:
