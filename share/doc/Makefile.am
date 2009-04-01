# Last modified: 2009-03-31.00

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

AM_MAKEINFOHTMLFLAGS = $(AM_MAKEINFOFLAGS) --no-split

# please note: The version*.texi files are automatically dependencies
info_TEXINFOS = corb-developer-handbook.txi corb-user-manual.txi

corb_developer_handbook_TEXINFOS = corb.txi         \
                                   corb-copying.txi \
                                   impatient.txi    \
                                   preface_dvl.txi  \
                                   onomastics.txi   \
                                   coding.txi       \
                                   licensing.txi    \
                                   documenting.txi  \
                                   buildsystem.txi  \
                                   indeces.txi      \
                                   fdl_wrapper.txi  \
                                   fdl-1.3.txi

corb_user_manual_TEXINFOS = corb.txi

# create pdf and html docs along with `make all`
all-local: pdf html

# install pdf and html docs along with `make install`
install-data-hook:
	@$(NORMAL_INSTALL)
	@$(MAKE) $(AM_MAKEFLAGS) install-pdf install-html

# distribute html and pdf docs along with this package
EXTRA_DIST = corb-developer-handbook.pdf corb-developer-handbook.html


## Local variables:
## eval: (add-hook 'write-file-hooks 'time-stamp)
## time-stamp-start: "Last modified: "
## time-stamp-format: "%:y-%02m-%02d.%02H"
## time-stamp-end: "$"
## End:
