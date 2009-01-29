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

@SET_MAKE@

ACLOCAL_AMFLAGS = -I share/misc/m4

SUBDIRS = libcrbfallback libcrbbasic libcrbrna libcrbbrot libcrbapps

noinst_HEADERS = crb_cmdline.h

# is there any particular reason we used LDADD and not corb_LDADD in the first
# place?
LDADD =                                 \
	libcrbapps/libcrbapps.a         \
	libcrbbrot/libcrbbrot.a         \
	libcrbrna/libcrbrna.a           \
	libcrbbasic/libcrbbasic.a       \
	libcrbfallback/libcrbfallback.a

#corb_LDADD =                                 \
#	libcrbapps/libcrbapps.a         \
#	libcrbbasic/libcrbbasic.a       \
#	libcrbbrot/libcrbbrot.a         \
#	libcrbrna/libcrbrna.a           \
#	libcrbfallback/libcrbfallback.a

bin_PROGRAMS = corb

corb_SOURCES =       \
	crb_cmdline.c\
	corb.c


## Local variables:
## eval: (add-hook 'write-file-hooks 'time-stamp)
## time-stamp-start: "Last modified: "
## time-stamp-format: "%:y-%02m-%02d.%02H"
## time-stamp-end: "$"
## End:
