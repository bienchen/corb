# Last modified: 2009-01-29.17

## Copyright (C) 2008 Stefan Bienert
##
## This file is part of CoRB.
##
## CoRB is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## CoRB is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with CoRB.  If not, see <http://www.gnu.org/licenses/>.


# This Makefile creates the core library of the CoRB project. This
# library holds standard functions, e.g. for file handling, memory management
# and error handling. Further, it defines some standard data types like arrays,
# stacks and hash tables.

@SET_MAKE@

ACLOCAL_AMFLAGS = -I share/misc/m4

AM_CPPFLAGS = -I$(top_srcdir)/src

noinst_LIBRARIES = libcrbrna.a

libcrbrna_a_SOURCES =                      \
	alphabet.c                         \
	nn_scores.c                        \
        secstruct.c                        \
	rna.c

noinst_HEADERS =                           \
	alphabet.h                         \
	nn_scores.h                        \
        secstruct.h                        \
	rna.h                              \
	crbrna.h

# tests for the library
LDADD =                                    \
	libcrbrna.a                        \
	../libcrbbasic/libcrbbasic.a       \
	../libcrbfallback/libcrbfallback.a

check_PROGRAMS =                           \
	test_alphabet                      \
	test_secstruct                     \
	test_nn_scores                     \
	test_rna

test_alphabet_SOURCES  = test_alphabet.c

test_nn_scores_SOURCES = test_nn_scores.c

test_secstruct_SOURCES = test_secstruct.c

test_rna_SOURCES       = test_rna.c

TESTS = $(check_PROGRAMS)

## Local variables:
## eval: (add-hook 'write-file-hooks 'time-stamp)
## time-stamp-start: "Last modified: "
## time-stamp-format: "%:y-%02m-%02d.%02H"
## time-stamp-end: "$"
## End:
