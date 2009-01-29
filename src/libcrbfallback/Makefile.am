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


# This Makefile creates a library of fallback functions. A "fallback function"
# is a function, which is not provided by the current system, but needed by
# this project. Functions covered by this library are declared via
# AC_REPLACE_FUNCS macros in configure.ac. If a function is not found by the
# configure script, an object is added to the library automatically. Of course
# an appropriate source file has to exist for that case. If all functions are
# available, we need a dummy object for the library, so the build system runs
# properly.

@SET_MAKE@

ACLOCAL_AMFLAGS = -I share/misc/m4

noinst_LIBRARIES = libcrbfallback.a

libcrbfallback_a_SOURCES = snprintf.c #dummy.c

noinst_HEADERS = snprintf.h

libcrbfallback_a_LIBADD  = @LIBOBJS@

## Local variables:
## eval: (add-hook 'write-file-hooks 'time-stamp)
## time-stamp-start: "Last modified: "
## time-stamp-format: "%:y-%02m-%02d.%02H"
## time-stamp-end: "$"
## End:
