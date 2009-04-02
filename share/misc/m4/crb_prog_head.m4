# Last modified: 2009-03-17.15

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

# CRB_PROG_HEAD
# -------------
# CRB_PROG_HEAD ()
# Tries to find the path to the head executable and stores it with the binary
# name in a variable head. If head is not found, an error is raised.
AC_DEFUN([CRB_PROG_HEAD],
[dnl# macro-body
 AC_PATH_PROG(CRB_HEAD, [head], no)
 AS_IF([test x$CRB_HEAD = xno],
       [AC_MSG_ERROR([Unable to find the head command])]
      )dnl# AS_IF
]dnl# macro-body
        )

dnl# Local variables:
dnl# eval: (add-hook 'write-file-hooks 'time-stamp)
dnl# time-stamp-start: "Last modified: "
dnl# time-stamp-format: "%:y-%02m-%02d.%02H"
dnl# time-stamp-end: "$"
dnl# End:
