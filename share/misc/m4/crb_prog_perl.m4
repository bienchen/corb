# Last modified: 2008-10-03.22

dnl Copyright (C) 2007 Stefan Bienert
dnl 
dnl This file is part of CoRB.
dnl 
dnl CoRB is free software: you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.
dnl 
dnl CoRB is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl 
dnl You should have received a copy of the GNU General Public License
dnl along with CoRB.  If not, see <http://www.gnu.org/licenses/>.

# CRB_PROG_PERL
# -------------
# CRB_PROG_PERL ()
# Tries to find the path to the Perl executable and stores it with the binary
# name in a variable PERL. If Perl is not found, an error is raised.
AC_DEFUN([CRB_PROG_PERL],
[dnl# macro-body
 AC_ARG_VAR(PERL, [Path and name of the Perl executable. Used for ]
                  [substitution in Perl scripts.])
 AC_PATH_PROG(PERL, perl)
]dnl# macro-body
        )


dnl# Local variables:
dnl# eval: (add-hook 'write-file-hooks 'time-stamp)
dnl# time-stamp-start: "Last modified: "
dnl# time-stamp-format: "%:y-%02m-%02d.%02H"
dnl# time-stamp-end: "$"
dnl# End:
