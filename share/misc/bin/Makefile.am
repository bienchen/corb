# Last modified: 2009-03-20.11

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

SUFFIXES=.pl

noinst_SCRIPTS=
if CRB_REFORMAT
  noinst_SCRIPTS += reformat
endif
if CRB_READ_VIENNAPARAMS
  noinst_SCRIPTS += read_viennaparams
endif
if CRB_EX_CMP_ER2DE_RNAEVAL
  noinst_SCRIPTS += ex_cmp_er2de_rnaeval
endif
if CRB_OPT_BROT_DESIGN_TERMS
  noinst_SCRIPTS += opt_brot_design_terms
  noinst_SCRIPTS += simplex_brot_params
endif

pl_mods = PBar.pm CorbIO.pm RNA.pm Simplex.pm

# misc scripts, build on --with|enable-development

CLEANFILES = $(noinst_SCRIPTS)

EXTRA_DIST = reformat.pl read_viennaparams.pl ex_cmp_er2de_rnaeval.pl \
	     simplex_brot_params.pl opt_brot_design_terms.pl $(pl_mods)

do_subst = sed \
  -e 's,[@]PERL[@],$(PERL),g' \
  -e 's,[@]SHELL[@],$(SHELL),g' \
  -e 's,[@]configure_input[@],Generated from $@.pl; do not edit by hand.,g' \
  -e 's,[@]CORB_PERL5LIB[@],$(abs_srcdir),g'

.pl:$(pl_mods)
	rm -f $@ .$@.tmp
	$(do_subst) $(srcdir)/$@.pl >.$@.tmp
	chmod +x .$@.tmp
	chmod a-w .$@.tmp
	mv -f .$@.tmp $@


## Local variables:
## eval: (add-hook 'write-file-hooks 'time-stamp)
## time-stamp-start: "Last modified: "
## time-stamp-format: "%:y-%02m-%02d.%02H"
## time-stamp-end: "$"
## End:
