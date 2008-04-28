#
# Copyright (C) 2008 Stefan Bienert
# Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
#
# See COPYING file in the top level directory of this tree for licence.
#
# Last modified: 2008-04-23.21
#


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
