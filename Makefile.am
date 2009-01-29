# Last modified: 2009-01-29.17

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

SUBDIRS=share src 

#  One target to rule them all, 
#    one target to check them,
#  one target to verify them all,
# and before the commit find them.
# Or in other words: A check to be performed before commiting anyting to git
check-commit:
	@echo "[ ============= Resetting to a virginal state ============= ]"
	@make -s maintainer-clean
	@echo "[ ===================== Bootstrapping ===================== ]"
	@$(top_srcdir)/bootstrap
	@echo "[ ============== Configuring for Development ============== ]"
	@$(top_srcdir)/configure --with-development --enable-development
	@echo "[ ================= Running  \`make check\` ================= ]"
	@make -s check
	@echo "[ =============== Running  \`make distcheck\` =============== ]"
	@make -s distcheck
	@echo "[ ====================== Cleaning up ====================== ]"
	@make -s clean
	@echo "[ ============ Configuring  Productive Version ============ ]"
	@$(top_srcdir)/configure
	@echo "[ ================= Running  \`make check\` ================= ]"
	@make -s check
	@echo "[ =============== Running  \`make distcheck\` =============== ]"
	@make -s distcheck


## Local variables:
## eval: (add-hook 'write-file-hooks 'time-stamp)
## time-stamp-start: "Last modified: "
## time-stamp-format: "%:y-%02m-%02d.%02H"
## time-stamp-end: "$"
## End:
