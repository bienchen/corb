# Last modified: 2008-10-03.22

dnl Copyright (C) 2008 Stefan Bienert
dnl
dnl This file is part of CoRB.
dnl
dnl CoRB is free software: you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.
dnl
dnl CoRB is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with CoRB.  If not, see <http://www.gnu.org/licenses/>.


# This Makefile creates a library of fallback functions. A "fallback function"
# is a function, which is not provided by the current system, but needed by
# this project. Functions covered by this library are declared via
# AC_REPLACE_FUNCS macros in configure.ac. If a function is not found by the
# configure script, an object is added to the library automatically. Of course
# an appropriate source file has to exist for that case. If all functions are
# available, we need a dummy object for the library, so the build system runs
# properly.


ACLOCAL_AMFLAGS = -I share/misc/m4

noinst_LIBRARIES = libcrbfallback.a

libcrbfallback_a_SOURCES = snprintf.c #dummy.c

noinst_HEADERS = snprintf.h

libcrbfallback_a_LIBADD  = @LIBOBJS@

dnl Local variables:
dnl eval: (add-hook 'write-file-hooks 'time-stamp)
dnl time-stamp-start: "Last modified: "
dnl time-stamp-format: "%:y-%02m-%02d.%02H"
dnl time-stamp-end: "$"
dnl End:
