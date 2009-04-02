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

# CRB_GET_VERSION (program, variable)
# -------------------------
# Fetch the version string of an application. For retrieving the version of a
# program, the following options are tested: "--version", "-v" and "-V". If the
# option of your desires is not in, just extend the list (and notify the author
# to adopt it). The detection mechanism just fetches the first value which
# "looks" like a version number (Major.Minro.Whatever). If a version is found,
# stores it in variable.
AC_DEFUN([CRB_GET_VERSION],
[dnl# macro-body
 _crb_version_options="--version -v -V"
 AC_REQUIRE([AC_PROG_GREP])
 AC_REQUIRE([AC_PROG_SED])
 AC_REQUIRE([CRB_PROG_HEAD])
 _crb_prog=`type -p $1 2>/dev/null | tail -n 1 | awk '{print $NF}'`
 AS_IF([test -n "$_crb_prog"],
       [dnl# run-if-true
        for _crb_opt in $_crb_version_options
          do
            _crb_ver=[`"$_crb_prog" $_crb_opt 2>&1                    \
                      | $GREP '\(^\| \)[0-9][0-9]*\.[0-9]'            \
                      | $CRB_HEAD -n 1                                \
                      | tr ' ' '\n'                                   \
                      | $SED 's/\([0-9][0-9]*\.[0-9][0-9.]*\).*/\1/g' \
                      | $GREP '\(^\| \)[0-9][0-9]*\.[0-9]'            \
                      | $CRB_HEAD -n 1`]
            AS_IF([test -n "$_crb_ver"],
                  [dnl# run-if-true
                   AS_TR_SH($2)=$_crb_ver
                   break
                  ]dnl# run-if-true
                 )dnl# AS_IF   
          done
       ],dnl# run-if-true
       [dnl# run-if-not-true
        AC_MSG_ERROR([program "$1" not found in PATH])
       ]dnl# run-if-not-true
      )dnl# AS_IF
 AS_IF([test -z "$AS_TR_SH($2)"],
       [dnl# run-if-true
        AC_MSG_ERROR([no version number found for program "$1"])        
       ]dnl# run-if-true
      )dnl# AS_IF
]dnl# macro-body
        )
 
dnl# Local variables:
dnl# eval: (add-hook 'write-file-hooks 'time-stamp)
dnl# time-stamp-start: "Last modified: "
dnl# time-stamp-format: "%:y-%02m-%02d.%02H"
dnl# time-stamp-end: "$"
dnl# End:
