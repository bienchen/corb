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

# CRB_CHECK_LIBCAIRO
# ------------------
# CRB_CHECK_LIBCAIRO (version, special-function, function1 function2 ...,
# [other-libraries])
# Checks for the presence of libcairo with extensive checking for all required
# functions. version should be a version number for libcairo which works with
# all requirements. special-function should be a function which is assumed to
# be in all versions of libcairo so that we can detect if the library is
# installed but of wrong version. function1 function2 ... is a whitespace
# separated list of all required functions. [other-libraries] are the
# libraries required to correctly link against libcairo.
AC_DEFUN([CRB_CHECK_LIBCAIRO],
[dnl# macro body
 # STEFAN
 AS_VAR_PUSHDEF([ac_Lib], [ac_cv_lib_cairo])
 # STEFAN
 AC_CACHE_CHECK([for -lcairo], ac_Lib,
                [dnl# check functions
                 ac_check_lib_save_LIBS=$LIBS
                 LIBS="-lcairo $4 $LIBS"
                 AC_LINK_IFELSE([AC_LANG_CALL([], [$2])],
                                              [AS_VAR_SET(ac_Lib, yes)],
 	                                      [AS_VAR_SET(ac_Lib, no)])
                 crb_cairo_fnc=
                 m4_foreach_w([__crb_cairo_fnc], [$3],
                            [dnl# expression
                             AS_IF([test AS_VAR_GET(ac_Lib) = yes],
                             [dnl# try to compile
                             AC_LINK_IFELSE([AC_LANG_CALL([], __crb_cairo_fnc)],
                                            [AS_VAR_SET(ac_Lib, yes)],
 	                                    [AS_VAR_SET(ac_Lib, no)])
                             AS_IF([test AS_VAR_GET(ac_Lib) = no],
                                   [dnl# set name of function not found
                                    crb_cairo_fnc=__crb_cairo_fnc
                                   ]
                                  )dnl# AS_IF
                             ]
                                  )dnl# AS_IF
                            ]
                           )dnl# m4_foreach
                 LIBS=$ac_check_lib_save_LIBS
                ]
               )
 AS_IF([test AS_VAR_GET(ac_Lib) = yes],
       [dnl# set -lcairo
        LIBS="-lcairo $LIBS"
       ],
       [dnl# action if false
        AS_IF([test "x$crb_cairo_fnc" != x],
              [dnl# action-if-true
               echo "         installed version of libcairo is missing function \`${crb_cairo_fnc}\`, recommended version: $1" >&AS_MESSAGE_FD
              ]
             )dnl# AS_IF
       ]
      )dnl# AS_IF
 AM_CONDITIONAL([CRB_WITH_CAIRO], [test AS_VAR_GET(ac_Lib) = yes])
 AS_VAR_POPDEF([ac_Lib])
]dnl# macro body
        )


dnl# Local variables:
dnl# eval: (add-hook 'write-file-hooks 'time-stamp)
dnl# time-stamp-start: "Last modified: "
dnl# time-stamp-format: "%:y-%02m-%02d.%02H"
dnl# time-stamp-end: "$"
dnl# End:
