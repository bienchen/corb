#
# Copyright (C) 2008 Stefan Bienert
# Copyright (C) 2008 Dominic Fabian
# Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
#
# See COPYING file in the top level directory of this tree for licence.
#
# Last modified: 2008-04-12.21
#


# CRB_CHECK_LIBCAIRO
# ------------------
# CRB_CHECK_LIBCAIRO ([other-libraries])
# $2: function
# $3: action-if-found
# $4: action-if-not-foun
# $1: other-libraries
AC_DEFUN([CRB_CHECK_LIBCAIRO],
[dnl# macro body
 AS_VAR_PUSHDEF([ac_Lib], [ac_cv_lib_cairo])
 AC_CACHE_CHECK([for -lcairo], ac_Lib,
                [dnl# check functions
                 ac_check_lib_save_LIBS=$LIBS
                 LIBS="-lcairo $1 $LIBS"
                 AS_VAR_SET(ac_Lib, yes)
                 m4_foreach([__crb_cairo_fnc], [dnl# list
                                                [function1],
                                                [function2]
                                               ],
                            [dnl# expression
                             AS_IF([test AS_VAR_GET(ac_Lib) = yes],
                             [dnl# try to compile
                             AC_LINK_IFELSE([AC_LANG_CALL([], __crb_cairo_fnc)],
                                            [AS_VAR_SET(ac_Lib, yes)],
 	                                    [AS_VAR_SET(ac_Lib, no)])
                             # throw msg in case test did not work
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
