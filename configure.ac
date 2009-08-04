# Last modified: 2009-08-04.15
#
dnl Copyright (C) 2008 Stefan Bienert
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
CRB_ARG_WITH_CUM([development],
                 [pc_elisp reformat read_viennaparams ex_cmp_er2de_rnaeval
                  opt_brot_design_terms brot_iterative eval_structure_file])

CRB_ARG_WITH_PC_ELISP
CRB_ARG_WITH_PERL_PROG([reformat], [build reformat tool for source code])
CRB_ARG_WITH_PERL_PROG([read_viennaparams],
                       [build conversion tool for Vienna energy tables])
CRB_ARG_WITH_PERL_PROG([ex_cmp_er2de_rnaeval],
                 [compare Nearest Neighbour model parameters to Vienna package])
CRB_ARG_WITH_PERL_PROG([opt_brot_design_terms],
                       [Optimiser for parameters in brot])
CRB_ARG_WITH_PERL_PROG([brot_iterative], [build tool to run brot iteratively])
CRB_ARG_WITH_PERL_PROG([eval_structure_file],
                       [build tool to test brot on a set of structures])

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
AC_REPLACE_FUNCS([strtoul])

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
        share/doc/Makefile
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


CRB_EVAL_VAR([__crb_bindir],  [$bindir])
CRB_EVAL_VAR([__crb_infodir], [$infodir])
CRB_EVAL_VAR([__crb_pdfdir],  [$pdfdir])
CRB_EVAL_VAR([__crb_htmldir], [$htmldir])

CRB_CFLAGS_CREATE_LIST([__crb_fmtd_cflags],
                       [                                               ])

AC_MSG_NOTICE(
[Finished successful

           Configuration Report
           ====================
           Build:
              ${PACKAGE_NAME} applications:

              development tools:
                       reformat:               $with_reformat
                       read_viennaparams:      $with_read_viennaparams
                       ex_cmp_er2de_rnaeval:   $with_ex_cmp_er2de_rnaeval
                       brot_iterative:         $with_brot_iterative
                       eval_structure_file:    $with_eval_structure_file
                       emacs lisp byte code:   $with_pc_elisp
           
           Build settings:
dnl              Perl:                            $PERL
              Enable memory checking:          $enable_memcheck
              Enable assertions:               $enable_assert
              Enable usage of POSIX threads:   $enable_pthreads
              CC:                              $CC
              CFLAGS:                          $__crb_fmtd_cflags
           Installation directories:
              Executables:     $__crb_bindir
              Manual (pdf):    $__crb_pdfdir
              Manual (info):   $__crb_infodir
              Manual (html):   $__crb_htmldir]
            )dnl# AC_MSG_NOTICE
 
dnl# Local variables:
dnl# eval: (add-hook 'write-file-hooks 'time-stamp)
dnl# time-stamp-start: "Last modified: "
dnl# time-stamp-format: "%:y-%02m-%02d.%02H"
dnl# time-stamp-end: "$"
dnl# End:
