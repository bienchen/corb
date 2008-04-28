# Last modified: 2008-04-22.10

# Copyright (c) 2007 Stefan Bienert <bienert@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg 

# See COPYING file in the top level directory of this tree for licence.

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

# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:
