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

# CRB_PROG
# ________
# CRB_PROG (program-to-find)
# Tries to find the path to "program-to-find" and stores it in a variable
# "CRB_program-to-find" with all-upper case letters as well as '-' and '.'
# replaced by '_'. If "program-to-find" is not found, AC_MSG_ERROR is raised.
AC_DEFUN([CRB_PROG],
[dnl# macro-body
 AC_PATH_PROG(CRB_[]m4_toupper(m4_translit([$1],[-.],[__])), [$1], no)
 AS_IF([test x$CRB_[]m4_toupper(m4_translit([$1],[-.],[__])) = xno],
       [AC_MSG_ERROR([Unable to find the $1 command])]
      )dnl# AS_IF

]dnl# macro-body
        )

# CRB_PROG_HEAD
# -------------
# CRB_PROG_HEAD ()
# Tries to find the path to the head executable and stores it with the binary
# name in a variable CRB_HEAD. If head is not found, an error is raised.
AC_DEFUN([CRB_PROG_HEAD],
[dnl# macro-body
 CRB_PROG([head])
]dnl# macro-body
        )

# CRB_PROG_CUT
# -------------
# CRB_PROG_CUT ()
# Tries to find the path to the cut executable and stores it with the binary
# name in a variable CRB_CUT. If cut is not found, an error is raised.
AC_DEFUN([CRB_PROG_CUT],
[dnl# macro-body
 CRB_PROG([cut])
]dnl# macro-body
        )

# CRB_PROG_UNIQ
# -------------
# CRB_PROG_UNIQ ()
# Tries to find the path to the uniq executable and stores it with the binary
# name in a variable CRB_UNIQ. If cut is not found, an error is raised.
AC_DEFUN([CRB_PROG_UNIQ],
[dnl# macro-body
 CRB_PROG([uniq])
]dnl# macro-body
        )

# CRB_PROG_PRINTF
# -------------
# CRB_PROG_PRINTF ()
# Tries to find the path to the printf executable and stores it with the binary
# name in a variable CRB_PRINTF. If printf is not found, an error is raised.
AC_DEFUN([CRB_PROG_PRINTF],
[dnl# macro-body
 CRB_PROG([printf])
]dnl# macro-body
        )

dnl# Local variables:
dnl# eval: (add-hook 'write-file-hooks 'time-stamp)
dnl# time-stamp-start: "Last modified: "
dnl# time-stamp-format: "%:y-%02m-%02d.%02H"
dnl# time-stamp-end: "$"
dnl# End:
