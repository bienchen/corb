# Last modified: 2008-10-06.14

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




# This Makefile creates the applications library. This is, the set of
# applications provided by CoRb. Therefore the library contains the
# main-function pendants of all the tools. See README for more information.


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

#dnl Local variables:
#dnl eval: (add-hook 'write-file-hooks 'time-stamp)
#dnl time-stamp-start: "Last modified: "
#dnl time-stamp-format: "%:y-%02m-%02d.%02H"
#dnl time-stamp-end: "$"
#dnl End:
