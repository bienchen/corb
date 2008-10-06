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



# This Makefile creates the core library of the CoRB project. This
# library holds standard functions, e.g. for file handling, memory management
# and error handling. Further, it defines some standard data types like arrays,
# stacks and hash tables.


ACLOCAL_AMFLAGS = -I share/misc/m4

AM_CPPFLAGS = -I$(top_srcdir)/src

noinst_LIBRARIES = libcrbbrot.a

libcrbbrot_a_SOURCES = \
	seqmatrix.c    \
	scmf_rna_opt.c

noinst_HEADERS =       \
	seqmatrix.h    \
        scmf_rna_opt.h \
	crbbrot.h

#dnl Local variables:
#dnl eval: (add-hook 'write-file-hooks 'time-stamp)
#dnl time-stamp-start: "Last modified: "
#dnl time-stamp-format: "%:y-%02m-%02d.%02H"
#dnl time-stamp-end: "$"
#dnl End:
