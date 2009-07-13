# Last modified: 2009-07-13.15

dnl Copyright (C) 2009 Stefan Bienert
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


# CRB_EVAL_VAR (variable, variable-to-expand)
# ________
# CRB_EVAL_VAR (variable, variable-to-expand)
# Takes a "variable-to-expand" and expands all variables in it iteratively. The
# result will be stored in "variable".
# Please note: Variable expansion will at max. be called 10 times.
# Please note, that "variable" must not contain ...
AC_DEFUN([CRB_EVAL_VAR],
[dnl# macro-body
 eval "$1=`eval "echo $2"`"
 dnl# extend if you find a variable which needs more evaluations
 for ____crb_null in 0 1 2 3 4 5 6 7 8 9 10
 do
   AS_CASE([@S|@$1],
           [*\@S|@*], [eval "$1=`eval "echo @S|@$1"`"],dnl#pattern1, if-matched1
           [break])dnl#default
 done
]dnl# macro-body
        )


dnl# Local variables:
dnl# eval: (add-hook 'write-file-hooks 'time-stamp)
dnl# time-stamp-start: "Last modified: "
dnl# time-stamp-format: "%:y-%02m-%02d.%02H"
dnl# time-stamp-end: "$"
dnl# End:
