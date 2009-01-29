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




# This Makefile creates the applications library. This is, the set of
# applications provided by CoRb. Therefore the library contains the
# main-function pendants of all the tools. See README for more information.

@SET_MAKE@

ACLOCAL_AMFLAGS = -I share/misc/m4

AM_CPPFLAGS = -I$(top_srcdir)/src

noinst_LIBRARIES = libcrbapps.a

libcrbapps_a_SOURCES =  \
	brot.c          \
	brot_cmdline.c  \
	fold.c          \
        fold_cmdline.c  \
        er2de.c         \
        er2de_cmdline.c

noinst_HEADERS =        \
	brot_cmdline.h  \
	brot.h \
	fold_cmdline.h  \
        fold.h          \
        er2de.h         \
        er2de_cmdline.h


# Do we need testing here?

## Local variables:
## eval: (add-hook 'write-file-hooks 'time-stamp)
## time-stamp-start: "Last modified: "
## time-stamp-format: "%:y-%02m-%02d.%02H"
## time-stamp-end: "$"
## End:
