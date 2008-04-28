# Last modified: 2008-04-24.11

# Copyright (c) 2007 Stefan Bienert <bienert@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg 

# See COPYING file in the top level directory of this tree for licence.


ACLOCAL_AMFLAGS = -I share/misc/m4

SUFFIXES=.pl

if CRB_REFORMAT
  bin_SCRIPTS = reformat
endif
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = reformat.pl

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


# Local variables:
# eval: (add-hook 'write-file-hooks 'time-stamp)
# time-stamp-start: "Last modified: "
# time-stamp-format: "%:y-%02m-%02d.%02H"
# time-stamp-end: "$"
# End:
