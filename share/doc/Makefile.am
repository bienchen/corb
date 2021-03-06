# Last modified: 2010-01-28.20


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

AUTOMAKE_OPTIONS = -Wno-override

ACLOCAL_AMFLAGS = -I share/misc/m4

AM_MAKEINFOHTMLFLAGS = $(AM_MAKEINFOFLAGS) --no-split

# please note: The version*.texi files are automatically dependencies
info_TEXINFOS = corb-user-manual.txi corb-developer-handbook.txi

corb_developer_handbook_TEXINFOS = corb.txi                      \
                                   corb-cpyng_wrapper_dvl.txi    \
                                   corb-copying.txi              \
                                   impatient.txi                 \
                                   preface_dvl.txi               \
                                   onomastics.txi                \
                                   abbreviations_wrapper_dvl.txi \
                                   abbreviations.txi             \
                                   using_git.txi                 \
                                   coding.txi                    \
                                   licensing.txi                 \
                                   documenting.txi               \
                                   buildsystem.txi               \
                                   indices-cp-dvl.txi            \
                                   indices-op-dvl.txi            \
                                   fdl_wrapper_dvl.txi           \
                                   fdl-1.3.txi

corb_user_manual_TEXINFOS = corb.txi                      \
                            corb-cpyng_wrapper_use.txi    \
                            corb-copying.txi              \
                            preface_use.txi               \
                            abbreviations_wrapper_use.txi \
                            abbreviations.txi             \
                            installation.txi              \
                            corb_tools.txi                \
                            brot-mission-pdf.pdf          \
                            brot-mission.png              \
                            and_everything.txi            \
                            vienna_notation-pdf.pdf       \
                            vienna_notation.png           \
                            indices-cp-use.txi            \
                            indices-pg-use.txi            \
                            indices-op-use.txi            \
                            indices-fig-use.txi           \
                            fdl_wrapper_use.txi           \
                            fdl-1.3.txi

# disable the dvi rule, since we do not support it (would need all pictures of
# the manuals in eps format)
dvi:

# create pdf and html docs along with `make all`
all-local: pdf html

# install pdf and html docs along with `make install`
install-data-hook:
	@$(NORMAL_INSTALL)
	@$(MAKE) $(AM_MAKEFLAGS) install-pdf install-html

# distribute html and pdf docs along with this package
EXTRA_DIST = corb-developer-handbook.pdf corb-developer-handbook.html \
             corb-user-manual.pdf corb-user-manual.html


## Local variables:
## eval: (add-hook 'write-file-hooks 'time-stamp)
## time-stamp-start: "Last modified: "
## time-stamp-format: "%:y-%02m-%02d.%02H"
## time-stamp-end: "$"
## End:
