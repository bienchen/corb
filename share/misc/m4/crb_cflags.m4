# Last modified: 2009-10-29.23

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

################################################################################
#######     Internal macros: Macros only to be used within this file     #######
################################################################################

# _CRB_CFLAGS_SET_WARNINGS_GCC(compiler, variable)
# ---------------------
# Set warning-flags for gcc. compiler is the gcc with or without its path.
# parameter 2 stores the version of gcc.
AC_DEFUN([_CRB_CFLAGS_SET_WARNINGS_GCC],
[dnl# macro-body
 CRB_GET_VERSION([$1], [$2])
 dnl# set default flags, assumed to work with everything lower than gcc-4.1.2,
 dnl# tested with gcc-3.3. If you encounter a gcc < 3.3, set up a comparison
 dnl# for gcc-3.3.
 CFLAGS="${CFLAGS} -Wall -Werror -pedantic -std=c99"
 CFLAGS="${CFLAGS} -Wstrict-prototypes -Wundef -Wshadow"
 CFLAGS="${CFLAGS} -Wcast-align -Wsign-compare"
 CFLAGS="${CFLAGS} -Wnested-externs -Winline"
 CFLAGS="${CFLAGS} -Wformat -Wformat-security -Wformat-y2k"
 CFLAGS="${CFLAGS} -Wbad-function-cast -Wcast-qual"
 AX_COMPARE_VERSION($AS_TR_SH($2), [ge], [4.1.2],
                    [dnl# action-if-true
                     CFLAGS="${CFLAGS} -Wextra -Winit-self"
                     CFLAGS="${CFLAGS} -Wdeclaration-after-statement"
                     CFLAGS="${CFLAGS} -Wmissing-include-dirs"
# what about -Wpointer-arith -Wc++-compat -Wredundant-decls?
# at least in 4.3 we have -Wno-variadic-macros. With this we could use -ansi
                    ]dnl# action-if-true
                   )dnl# AX_COMPARE_VERSION
]dnl# macro-body
        )


# _CRB_CFLAGS_SET_WARNINGS_ICC(compiler, variable)
# ---------------------
# Set warning-flags for icc. compiler is the icc with or without its path.
# parameter 2 stores the version of icc.
AC_DEFUN([_CRB_CFLAGS_SET_WARNINGS_ICC],
[dnl# macro-body
 CRB_GET_VERSION([$1], [$2])
 CFLAGS="${CFLAGS} -std=c99 -Wall -Wcheck -diag-enable port-win -fstack-security-check -Wdeprecated -Wextra-tokens -Wformat -Winline -Wmissing-declarations -Wmissing-prototypes -Wreturn-type -Wshadow -Wstrict-prototypes -Wuninitialized -Wunknown-pragmas -Wunused-function -Wwrite-strings -Werror -no-gcc -wd981"
# do we need '-debug full' for debugging? -Werror -no-gcc -Wunused-variable
]dnl# macro-body
        )

# _CRB_CFLAGS_SET_WARNINGS_SUNCC(compiler, variable)
# ---------------------
# Set warning-flags for gcc. compiler is the gcc with or without its path.
# parameter 2 stores the version of gcc.
AC_DEFUN([_CRB_CFLAGS_SET_WARNINGS_SUNCC],
[dnl# macro-body
 CRB_GET_VERSION([$1], [$2])
 CFLAGS="${CFLAGS} -xO2 -Xc -errfmt=error -errwarn"
]dnl# macro-body
        )

# _CRB_CFLAGS_INVOKE_POSIX
# ---------------------------
# Set the the macro _XOPEN_SOURCE to a value of 600. This invokes all POSIX
# compliant functions (ISO C) plus the X/Open System Interfaces extensions
# (XSI). The macro has to be defined at the compiler level since it is
# required to be set before any header is included. See
# http://www.opengroup.org/onlinepubs/009695399/functions/xsh_chap02_02.html
# for further informations.
AC_DEFUN([_CRB_CFLAGS_INVOKE_POSIX],
[dnl# macro-body
 CFLAGS="${CFLAGS} -D_XOPEN_SOURCE=600"
]dnl# macro-body
        )

################################################################################
#########        Public macros: Macros to set certain CFLAGS           #########
################################################################################

# CRB_CFLAGS_SET_WARNINGS
# -------------------
# Set warning-flags for the compiler. The options to request or suppress
# warnings go immedeately into the CFLAGS variable. At the moment, gcc, icc and
# sun cc are supported. For unknown compilers nothing will be set. Should be
# called only after all AC_REPLACE.. macros or other checks which include
# compiling. Most of such tests will not work with these compiler settings.
AC_DEFUN([CRB_CFLAGS_SET_WARNINGS],
[dnl# macro-body
 _crb_cc_basename=["`expr "//${CC#distcc}" : '.*/[[:space:]]*\([^/]*\)'`"]
 AC_MSG_CHECKING([version of $CC])
 AS_CASE([$_crb_cc_basename],
         [gcc], [_CRB_CFLAGS_SET_WARNINGS_GCC([$CC], [_crb_ver])],
         [icc], [_CRB_CFLAGS_SET_WARNINGS_ICC([$CC], [_crb_ver])],
         [cc],  [_CRB_CFLAGS_SET_WARNINGS_SUNCC([$CC], [_crb_ver])],
         [_crb_ver="${CC} unknown, no automatic development flag setting"]
        )dnl# AS_CASE

 dnl# _XOPEN_SOURCE is not a compiler option but a macro to be defined for all
 dnl# compilers
 _CRB_CFLAGS_INVOKE_POSIX
 AC_MSG_RESULT([$_crb_ver])
]dnl# macro-body
        )

# CRB_CFLAGS_CREATE_LIST (variable, offset)
# ---------------
# Convert the CFLAGS into a list ready for output. Takes an offset as argument
# to create aligned output. This offset has to be a string. The list is stored
# in variable. If CFLAGS is empty,
# variable is set to "none"
AC_DEFUN([CRB_CFLAGS_CREATE_LIST],
[dnl# macro-body
 for _crb_cflag in $CFLAGS
   do
AS_TR_SH($1)="$_crb_cflag
$2$[]AS_TR_SH($1)"
   done
 AS_IF([test -z "$AS_TR_SH($1)"],
       [dnl# run-if-true
        AS_TR_SH($1)=none
       ]dnl# run-if-true
      )dnl# AS_IF
 # then handle empty list
]dnl# macro-body
        )

dnl# Local variables:
dnl# eval: (add-hook 'write-file-hooks 'time-stamp)
dnl# time-stamp-start: "Last modified: "
dnl# time-stamp-format: "%:y-%02m-%02d.%02H"
dnl# time-stamp-end: "$"
dnl# End:
