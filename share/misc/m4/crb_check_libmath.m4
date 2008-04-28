#
# Copyright (C) 2008 Stefan Bienert
# Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
#
# See COPYING file in the top level directory of this tree for licence.
#
# Last modified: 2008-04-12.21
#


# CRB_CHECK_LIBMATH
# -----------------
# check for math library on different operating systems (will be extended with
# new OS as we learn them). If library is found, LIB is modified.
AC_DEFUN([CRB_CHECK_LIBMATH],
[dnl# macro body
 AC_REQUIRE([AC_CANONICAL_HOST])dnl
 __crb_LIBMATH=
 AS_CASE($host,
         [*-*-beos*],      [],
         [*-*-cygwin*],    [],
         [*-*-pw32*],      [],
         [*-*-darwin*],    [],
         [*-ncr-sysv4.3*], [dnl
                            AC_CHECK_LIB(mw, _mwvalidcheckl,
                                         __crb_LIBMATH="-lmw")
                            AC_CHECK_LIB(m, cos,
                                        __crb_LIBMATH="$__crb_LIBMATH -lm")  
                           ],
         [dnl# default
          AC_CHECK_LIB(m, log10, __crb_LIBMATH="-lm")
         ]
        )dnl# AS_CASE
 LIBS="$LIBS $__crb_LIBMATH"
]dnl# macro body
         )
 
dnl# Local variables:
dnl# eval: (add-hook 'write-file-hooks 'time-stamp)
dnl# time-stamp-start: "Last modified: "
dnl# time-stamp-format: "%:y-%02m-%02d.%02H"
dnl# time-stamp-end: "$"
dnl# End:
