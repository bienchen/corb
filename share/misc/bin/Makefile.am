# Last modified: 2008-10-30.14

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



ACLOCAL_AMFLAGS = -I share/misc/m4

SUFFIXES=.pl

bin_SCRIPTS=
if CRB_REFORMAT
  bin_SCRIPTS += reformat
endif
if CRB_READ_VIENNAPARAMS
  bin_SCRIPTS += read_viennaparams
endif
if CRB_EX_CMP_ER2DE_RNAEVAL
  bin_SCRIPTS += ex_cmp_er2de_rnaeval
endif

# misc scripts, build on --with|enable-development

CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = reformat.pl read_viennaparams.pl

do_subst = sed \
  -e 's,[@]PERL[@],$(PERL),g' \
  -e 's,[@]SHELL[@],$(SHELL),g' \
  -e 's,[@]configure_input[@],Generated from $@.pl; do not edit by hand.,g'

.pl:
	rm -f $@ .$@.tpm
	$(do_subst) $@.pl >.$@.tmp
	chmod +x .$@.tmp
	chmod a-w .$@.tmp
	mv -f .$@.tmp $@


#dnl Local variables:
#dnl eval: (add-hook 'write-file-hooks 'time-stamp)
#dnl time-stamp-start: "Last modified: "
#dnl time-stamp-format: "%:y-%02m-%02d.%02H"
#dnl time-stamp-end: "$"
#dnl End:









