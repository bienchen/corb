# Last modified: 2009-10-29.21

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

# CRB_PERL_PROG
# -------------
# CRB_PERL_PROG ()
# Tries to find the path to the Perl executable and stores it with the binary
# name in a variable PERL.
AC_DEFUN([CRB_PERL_PROG],
[dnl# macro-body
 AC_ARG_VAR(PERL, [Path and name of the Perl executable. Used for ]
                  [substitution in Perl scripts.])
 AC_PATH_PROG(PERL, [perl])
]dnl# macro-body
        )

# CRB_CHECK_MATCHING_PERL
# -----------------------
# CRB_CHECK_MATCHING_PERL ()
# Check that Perl is available in the version needed by CoRB's Perl scripts.
AC_DEFUN([CRB_CHECK_MATCHINGPERL],
[dnl# macro-body
  AC_REQUIRE([CRB_PERL_PROG])
  AS_IF([test -z $PERL],
        AC_MSG_ERROR([perl not found])dnl# run-if-true
       )dnl# AS_IF

  AS_IF([$PERL -e 'require 5.006;'],
        [],dnl# run-if-true
        [AC_MSG_ERROR([Perl v5.6 or better is required; Perl v5.8.2 or better is recommended. If you have several perl versions installed, select the one Automake should use using ./configure PERL=/path/to/perl])]dnl# run-if-false
       )dnl# AS_IF

]dnl# macro-body 
        )

# CRB_PERL_CHECK_MODULE
# ---------------------
# CRB_PERL_CHECK_MODULE (module[, run-if-true][, run-if-false])
# Check wheter Perl can find a certain module.
# Please note: We do not test for the Perl binary, we just use $PERL.
AC_DEFUN([CRB_PERL_CHECK_MODULE],
[dnl# macro-body
 AC_MSG_CHECKING([for Perl module $1])
 $PERL -e 'require $1;' 2>/dev/null
 AS_IF([test $? -eq 0],
       AC_MSG_RESULT(yes)
       m4_ifvaln([$2],[$2]),dnl# run-if-true
       AC_MSG_RESULT(no)
       m4_ifvaln([$3],[$3])dnl # run-if-false
      )dnl# AS_IF

]dnl# macro-body 
        )


dnl# Local variables:
dnl# eval: (add-hook 'write-file-hooks 'time-stamp)
dnl# time-stamp-start: "Last modified: "
dnl# time-stamp-format: "%:y-%02m-%02d.%02H"
dnl# time-stamp-end: "$"
dnl# End:
