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
