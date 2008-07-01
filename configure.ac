#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
# Last modified: 2008-06-30.15
#

AC_PREREQ(2.60)
AC_INIT([corb], [0.1], [Stefan Bienert <bienert@zbh.uni-hamburg.de>])

AC_CONFIG_AUX_DIR([share/misc/config])
AC_CONFIG_MACRO_DIR([share/misc/m4])


AM_INIT_AUTOMAKE([filename-length-max=99 -Wall -Werror])

AM_PATH_LISPDIR

AC_CONFIG_SRCDIR([src/corb.c])
AC_CONFIG_HEADERS(src/config.h:src/config.h.in)

AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# Checks for programs.
AC_PROG_SED
AC_PROG_GREP

# development options
# development features
CRB_ARG_ENABLE_CUM([development], [cflags-warnings assert memcheck])

CRB_ARG_ENABLE_CFLAGS_WARNINGS_OPT
CRB_ARG_ENABLE_MEMCHECK
CRB_ARG_ENABLE_ASSERT

# development tools
CRB_ARG_WITH_CUM([development], [pc-elisp reformat])

CRB_ARG_WITH_PC_ELISP
CRB_ARG_WITH_REFORMAT

# Checks for header files
CRB_CHECK_HEADERS([assert.h float.h limits.h stddef.h])
AC_HEADER_STDBOOL
gl_SIZE_MAX() # check for SIZE_MAX macro

# Checks for libraries.
CRB_ARG_ENABLE_CHECK_PTHREADS
CRB_CHECK_LIBMATH
#CRB_CHECK_LIBCAIRO([1.5.1], [cairo_create], [function1 function2])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC  # AC_REPLACE_FUNCS(malloc) with check for GNU libc compatibility
AC_FUNC_REALLOC # like AC_FUNC_MALLOC
#AC_CHECK_FUNCS([snprintf], , [AC_LIBOBJ(snprintf)
#                              AC_REPLACE_FUNCS([memchr])
#                              AC_REPLACE_FUNCS([strchr])])
AC_CHECK_FUNCS([snprintf],[AC_DEFINE([HAVE_SNPRINTF],
                          [1],
                          [Define to 1 if you have the `snprintf' function.])])
AC_DEFINE([PREFER_PORTABLE_SNPRINTF],
          [1],
          [Define to 1 to get a portable snprintf.])
AC_REPLACE_FUNCS([memset])
AC_REPLACE_FUNCS([vprintf])
AC_REPLACE_FUNCS([vfprintf])
AC_REPLACE_FUNCS([strerror])


CRB_ARG_ENABLE_CFLAGS_WARNINGS_EVAL


# enhance config.h
AH_TOP([
#ifdef __cplusplus
extern "C" {
#endif

#ifndef CRB_CONFIG_H
#define CRB_CONFIG_H])
AH_BOTTOM([
#endif  /* CRB_CONFIG_H */

#ifdef __cplusplus
}
#endif])

AC_CONFIG_FILES([
        Makefile
        share/Makefile
        share/misc/Makefile
        share/misc/bin/Makefile
        share/misc/emacs/Makefile
        src/libcrbfallback/Makefile
        src/libcrbbasic/Makefile
        src/libcrbrna/Makefile
        src/libcrbbrot/Makefile
        src/libcrbapps/Makefile
        src/Makefile
                ])

AC_OUTPUT


eval __crb_bindir=`eval "echo $bindir"`

AC_MSG_NOTICE(
[Finished successful

           Configuration Report
           ====================
           Build:
              ${PACKAGE_NAME} applications:

              development tools:
                       reformat:               $with_reformat
                       emacs lisp byte code:   $with_pc_elisp
           
           Build settings:
              Enable paranoid CFLAGS:          $enable_cflags_warnings
              Enable memory checking:          $enable_memcheck
              Enable assertions:               $enable_assert
              Enable usage of POSIX threads:   $enable_pthreads

           Installation directories:
              Executables:   $__crb_bindir]
            )dnl# AC_MSG_NOTICE

 
dnl# Local variables:
dnl# eval: (add-hook 'write-file-hooks 'time-stamp)
dnl# time-stamp-start: "Last modified: "
dnl# time-stamp-format: "%:y-%02m-%02d.%02H"
dnl# time-stamp-end: "$"
dnl# End:
