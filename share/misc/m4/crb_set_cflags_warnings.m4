#
# Copyright (C) 2008 Stefan Bienert
# Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
#
# See COPYING file in the top level directory of this tree for licence.
#
# Last modified: 2008-04-12.21
#

# _CRB_SET_CFLAGS_WARNINGS_GCC(compiler, variable)
# ---------------------
# Set warning-flags for gcc. compiler is the gcc with or without its path.
# parameter 2 stores the version of gcc.
AC_DEFUN([_CRB_SET_CFLAGS_WARNINGS_GCC],
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
                    ]dnl# action-if-true
                   )dnl# AX_COMPARE_VERSION
]dnl# macro-body
        )


# _CRB_SET_CFLAGS_WARNINGS_ICC(compiler, variable)
# ---------------------
# Set warning-flags for icc. compiler is the icc with or without its path.
# parameter 2 stores the version of icc.
AC_DEFUN([_CRB_SET_CFLAGS_WARNINGS_ICC],
[dnl# macro-body
 CRB_GET_VERSION([$1], [$2])
 CFLAGS="${CFLAGS} -std=c99 -Wall -Wcheck -diag-enable port-win -fstack-security-check -Wdeprecated -Wextra-tokens -Wformat -Winline -Wmissing-declarations -Wmissing-prototypes -Wreturn-type -Wshadow -Wstrict-prototypes -Wuninitialized -Wunknown-pragmas -Wunused-function -Wwrite-strings -Werror -no-gcc -wd981"
# do we need '-debug full' for debugging? -Werror -no-gcc -Wunused-variable
]dnl# macro-body
        )

# _CRB_SET_CFLAGS_WARNINGS_SUNCC(compiler, variable)
# ---------------------
# Set warning-flags for gcc. compiler is the gcc with or without its path.
# parameter 2 stores the version of gcc.
AC_DEFUN([_CRB_SET_CFLAGS_WARNINGS_SUNCC],
[dnl# macro-body
 CRB_GET_VERSION([$1], [$2])
 CFLAGS="${CFLAGS} -xO2 -Xc -errfmt=error -errwarn"
]dnl# macro-body
        )

# CRB_SET_CFLAGS_WARNINGS
# -------------------
# Set warning-flags for the compiler. The options to request or suppress
# warnings go immedeately into the CFLAGS variable. At the moment, gcc, icc and
# sun cc are supported. For unknown compilers nothing will be set. Should be
# called only after all AC_REPLACE.. macros or other checks which include
# compiling. Most of such tests will not work with these compiler settings.
AC_DEFUN([CRB_SET_CFLAGS_WARNINGS],
[dnl# macro-body
 _crb_cc_basename=["`expr "//$CC" : '.*/\([^/]*\)'`"]
 AC_MSG_CHECKING([version of $CC])
 AS_CASE([$_crb_cc_basename],
         [gcc], [_CRB_SET_CFLAGS_WARNINGS_GCC([$CC], [_crb_ver])],
         [icc], [_CRB_SET_CFLAGS_WARNINGS_ICC([$CC], [_crb_ver])],
         [cc],  [_CRB_SET_CFLAGS_WARNINGS_SUNCC([$CC], [_crb_ver])],
         [_crb_ver="${CC} unknown, no automatic development flag setting"]
        )dnl# AS_CASE
 AC_MSG_RESULT([$_crb_ver])
]dnl# macro-body
        )
 
dnl# Local variables:
dnl# eval: (add-hook 'write-file-hooks 'time-stamp)
dnl# time-stamp-start: "Last modified: "
dnl# time-stamp-format: "%:y-%02m-%02d.%02H"
dnl# time-stamp-end: "$"
dnl# End:
