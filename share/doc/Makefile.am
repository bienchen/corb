# Last modified: 2008-11-29.13

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

AM_MAKEINFOHTMLFLAGS = $(AM_MAKEINFOFLAGS) --no-split

# please note: The version*.texi files are automatically dependencies
info_TEXINFOS = corb-developer-handbook.txi corb-user-manual.txi

corb_developer_handbook_TEXINFOS = corb.txi        \
                                   impatient.txi   \
                                   preface_dvl.txi \
                                   onomastics.txi  \
                                   coding.txi      \
                                   indeces.txi     \
                                   fdl_wrapper.txi

corb_user_manual_TEXINFOS = corb.txi

# create pdf and html docs along with `make all`
all-local: pdf html

# install pdf and html docs along with `make install`
install-data-hook:
	@$(NORMAL_INSTALL)
	@$(MAKE) $(AM_MAKEFLAGS) install-pdf install-html

# distribute html and pdf docs along with this package
EXTRA_DIST = corb-developer-handbook.pdf corb-developer-handbook.html


#dnl Local variables:
#dnl eval: (add-hook 'write-file-hooks 'time-stamp)
#dnl time-stamp-start: "Last modified: "
#dnl time-stamp-format: "%:y-%02m-%02d.%02H"
#dnl time-stamp-end: "$"
#dnl End:
