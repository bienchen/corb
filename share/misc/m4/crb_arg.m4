#    a
# Copyright (C) 2008 Stefan Bienert
# Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
#
# See COPYING file in the top level directory of this tree for licence.
#
# Last modified: 2008-04-26.23
#


################################################################################
#######     Internal macros: Macros only to be used within this file     #######
################################################################################

# _CRB_ARG
# --------
# _CRB_ARG (type, prefix, option, message, help-text)
# Internal macro. Simple wrapper for AC_ARG_"type". Instead of the
# "action-if-given" and "action-if-not-given" cases, a variable "type"_"option"
# is set to 'yes' or 'no'. As argument the resulting options only take these
# values. As help string, only the help text should be given. The string is
# invoked using AS_HELP_STRING([--"prefix"-"option"], ["help-text"]). In the
# end, AC_MSG_RESULT is called with either 'yes' or 'no'. Before AC_ARG_"type",
# AC_MSG_CHECKING is set with "message".
AC_DEFUN([_CRB_ARG],
[dnl# macro-body
 AC_MSG_CHECKING($4)
 AC_ARG_$1($3,
               AS_HELP_STRING([--$2-$3], $5),
               [dnl# action-if-given
                AS_CASE($$2val,
                      [yes], [m4_tolower($2)_[]m4_translit([$3],[-.],[__])=yes],
                      [no],  [m4_tolower($2)_[]m4_translit([$3],[-.],[__])=no],
                        [dnl# default
        AC_MSG_ERROR([bad value "${m4_tolower($2)val}" for --$2-$3])
                        ]dnl
                       )dnl# AS_CASE
               ],
               [dnl# action-if-not-given
                m4_tolower($2)_[]m4_translit([$3],[-.],[__])=no
               ]
              )dnl# AC_ARG_$1
 AC_MSG_RESULT([$m4_tolower($2)_[]m4_translit([$3],[-.],[__])])
]dnl# macro-body
        )

# _CRB_ARG_WWO
# ------------
# _CRB_ARG_WWO (application, message, help-string, type)
# Wrapper for AC_ARG_WITH. Instead of the "action-if-given" and
# "action-if-not-given" cases, a variable with_"application" is set to 'yes' or
# 'no'. As argument the resulting options only take these values. As help
# string, only the help text should be given. The string is invoked using
# AS_HELP_STRING([--"type"-"applicaiton"], ["help-text"]). In the end,
# AC_MSG_RESULT is called with either 'yes' or 'no'. Before AC_ARG_WITH,
# AC_MSG_CHECKING is set with "message".
AC_DEFUN([_CRB_ARG_WWO],
[dnl# macro-body
 _CRB_ARG([WITH], $4, $1, $2, $3)
]dnl# macro-body
        )

# _CRB_ARG_CHECKED_LOOP_SET
# -------------------------
# _CRB_ARG_CHECKED_LOOP_SET (type, feature1 feature2 ..., [value])
# Internal macro. Loops over a whitespace separated list, creates a variable
# "type_<list itme>" for each member of the list and sets it to "value".
AC_DEFUN([_CRB_ARG_CHECKED_LOOP_SET],
[dnl# macro-body
 m4_foreach_w([crb_item],
              m4_translit([$2],[-.],[__]),
              [dnl# expression
               AS_IF([test "${$1_[]crb_item+set}" = set],
                     [],dnl# run-if-true
                     [dnl# run-if-false
                      $1_[]crb_item=$3
                     ]
                    )dnl# AS_IF
              ]dnl# expression
             )dnl# m4_foreach_w                
]dnl# macro-body
        )


################################################################################
######     General macros: Comfortable wrappers for AC_ARG_... macros     ######
################################################################################

# CRB_ARG_WITH
# ------------
# CRB_ARG_WITH (application, message, help-string)
# Simple wrapper for AC_ARG_WITH. Instead of the "action-if-given" and
# "action-if-not-given" cases, a variable with_"application" is set to 'yes' or
# 'no'. As argument the resulting options only take these values. As help
# string, only the help text should be given. The string is invoked using
# AS_HELP_STRING([--with-"applicaiton"], ["help-text"]). In the end,
# AC_MSG_RESULT is called with either 'yes' or 'no'. Before AC_ARG_WITH,
# AC_MSG_CHECKING is set with "message".
AC_DEFUN([CRB_ARG_WITH],
[dnl# macro-body
  _CRB_ARG_WWO($1, $2, $3, [with])
]dnl# macro-body
        )


# CRB_ARG_WITHOUT
# ---------------
# CRB_ARG_WITHOUT (application, message, help-string)
# Simple wrapper for AC_ARG_WITH. Instead of the "action-if-given" and
# "action-if-not-given" cases, a variable with_"application" is set to 'yes' or
# 'no'. As argument the resulting options only take these values. As help
#  string, only the help text should be given. The string is invoked using
# AS_HELP_STRING([--without-"applicaiton"], ["help-text"]). In the end,
# AC_MSG_RESULT is called with either 'yes' or 'no'. Before AC_ARG_WITH,
# AC_MSG_CHECKING is set with "message".
AC_DEFUN([CRB_ARG_WITHOUT],
[dnl# macro-body
  _CRB_ARG_WWO($1, $2, $3, [without])
]dnl# macro-body
        )

# CRB_ARG_ENABLE
# --------------
# CRB_ARG_ENABLE (feature, message, help-string)
# Simple wrapper for AC_ARG_ENABLE. Instead of the "action-if-given" and
# "action-if-not-given" cases, a variable enable_"feature" is set to 'yes' or
# 'no'. As argument the resulting options only take these values. As help
# string, only the help text should be given. The string is invoked using
# AS_HELP_STRING([--enable-"feature"], ["help-text"]). In the end,
# AC_MSG_RESULT is called with either 'yes' or 'no'. Before AC_ARG_ENABLE,
# AC_MSG_CHECKING is set with "message".
AC_DEFUN([CRB_ARG_ENABLE],
[dnl# macro-body
 _CRB_ARG([ENABLE], [enable], $1, $2, $3)
]dnl# macro-body
        )

# CRB_ARG_ENABLE_CUM
# ------------------
# CRB_ARG_ENABLE_CUM ([cumulative-feature], [feature1 ...])
# Creates a feature (--enable-"cumulative-feature") to enable a whole bunch of
# features at once. That is, all variables $enable_"feature1" ... are set
# explicitly to 'yes'. Invoking "cumulative-feature" as a negociation will set
# all variables to 'no'. Does not override members which are set for themself.
# CRB_ARG_ENABLE_CUM has to be place BEFORE any setting of features in
# configure.ac.
AC_DEFUN([CRB_ARG_ENABLE_CUM],
[dnl# macro-body
 CRB_ARG_ENABLE($1,
                [whether to enable all $1 features],
                [invoke]
                m4_foreach_w([crb_feature],
                             [$2],
                             ["--enable-[]crb_feature" ])
               )dnl# CRB_ARG_ENABLE
 dnl# after CRB_ARG_ENABLE, $enable_$1 can only by 'yes' or 'no'
 _CRB_ARG_CHECKED_LOOP_SET([enable],
                           $2,
                           [$enable_[]m4_translit([$1],[-.],[__])]
                          )dnl# _CRB_ARG_CHECKED_LOOP_SET
]dnl# macro-body
        )

# CRB_ARG_WITH_CUM
# ----------------
# CRB_ARG_WITH ([cumulative-package], [package1 ...])
# Creates a option (--with-"cumulative-package") to enable a whole bunch of
# packages to be build at once. That is, all variables $with_"package1" ... are
# set explicitly to 'yes'. Invoking "cumulative-package" as a negociation will
# set all variables to 'no'. Does not override members which are set for
# themself. CRB_ARG_WITH_CUM has to be place BEFORE any setting of packages in
# configure.ac.
AC_DEFUN([CRB_ARG_WITH_CUM],
[dnl# macro-body
 CRB_ARG_WITH($1,
                [whether to enable all $1 packages],
                [invoke]
                m4_foreach_w([crb_feature],
                             [$2],
                             ["--with-[]crb_feature" ])
               )dnl# CRB_ARG_WITH
 dnl# after CRB_ARG_WITH, $with_$1 can only by 'yes' or 'no'
 _CRB_ARG_CHECKED_LOOP_SET([with],
                           $2,
                           [$with_[]m4_translit([$1],[-.],[__])]
                          )dnl# _CRB_ARG_CHECKED_LOOP_SET
]dnl# macro-body
        )
 

################################################################################
#########        Explicit macros: Macros with complete function        #########
################################################################################

# CRB_ARG_ENABLE_ASSERT
# ---------------------
# CRB_ARG_ENABLE_ASSERT ()
# Check whether to enable assertions.
# Creates option "--enable-assert" for the configure script to check whether
# to enable assertions or not. If set, in fact nothing happens. If not set, C
# macro NDEBUG is defined in the configurational header (config.h). This
# signals the compiler to turn off assertions. The result of the checks'
# evalation is stored as "yes" or "no" in a variable enable_assert.
# Don't forget to check for the header by yourself!
AC_DEFUN([CRB_ARG_ENABLE_ASSERT],
[dnl# macro-body
 CRB_ARG_ENABLE([assert],
                [whether to enable assertions],
                [turn on assertions]
               )dnl# CRB_ARG_ENABLE
 AS_IF([test "x$enable_assert" = xyes],
       [],dnl# do nothing
       [dnl# run-if-true
        AC_DEFINE([NDEBUG],
                  [1],
                  [Define to 1 if assertions should be disabled.])
       ]
      )dnl# AS_IF         
]dnl# macro-body
        )

# CRB_ARG_ENABLE_CFLAGS_WARNINGS
# ------------------------------
# CRB_ARG_ENABLE_CFLAGS_WARNINGS ()
# The idea is to enable a set of very pedantic compiler flags for checking C
# code. This is considered to be usefull during development. Since the CFLAGS
# should only be changed after all autoconf checks for functions, libraries,
# etc. took place, enabling the pedantic flags is splitted into the
# option-recording macro and the evaluation of the option.
# At this moment, gcc, intels' icc and suns' cc are supported by specialised
# settings.
AC_DEFUN([CRB_ARG_ENABLE_CFLAGS_WARNINGS_OPT],
[dnl# macro-body
 CRB_ARG_ENABLE([cflags-warnings],
                [whether to enable pedantic C code checking],
                [turn on pedantic C code verification]
               )dnl# CRB_ARG_ENABLE
]dnl# macro-body
        )

# CRB_ARG_ENABLE_CFLAGS_WARNINGS_EVAL
# -----------------------------------
# CRB_ARG_ENABLE_CFLAGS_WARNINGS_EVAL ()
# Evaluates whether "--enable-cflags-warnings" was set and acts on it.
AC_DEFUN([CRB_ARG_ENABLE_CFLAGS_WARNINGS_EVAL],
[dnl# macro-body
 AS_IF([test "x$enable_cflags_warnings" = xyes],
       [CRB_SET_CFLAGS_WARNINGS]dnl# run-if-true
      )dnl# AS_IF

]dnl# macro-body
        )

# CRB_ARG_ENABLE_CHECK_PTHREADS
# -----------------------------
# CRB_ARG_ENABLE_CHECK_PTHREADS ()
# Create a option "--enable-pthreads". If set, ACX_PTHREAD is called to setup
# compilation with POSIX threads. If pthreads are available, proper CFLAGS are
# set, LIBS is extended and a C macro HAVE_PTHREAD is defined. The result of
# the checks' evalation is stored as 'yes' or 'no' in a variable
# enable_pthreads.
AC_DEFUN([CRB_ARG_ENABLE_CHECK_PTHREADS],
[dnl# macro-body
 CRB_ARG_ENABLE([pthreads],
                [whether to check for POSIX threads],
                [check and invoke POSIX threads in C programs]
               )dnl# CRB_ARG_ENABLE
 AS_IF([test "x$enable_pthreads" == xyes],
       [dnl# action-if-true
        ACX_PTHREAD([LIBS="$PTHREAD_LIBS $LIBS"
                     CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
                     CC="$PTHREAD_CC"
                     AC_DEFINE([HAVE_PTHREAD],
                               [1],
                  [Define if you have POSIX threads libraries and header files.]
                              )
                    ]
                   )
       ]
      )dnl# AS_IF        
]dnl# macro-body
        )

# CRB_ARG_ENABLE_MEMCHECK
# -----------------------
# CRB_ARG_ENABLE_MEMCHECK ()
# Create a feature option (--enable-memcheck) to indicate whether to
# enable memory checking (leaks & double free's). Defines a C macro MEMCHECK in
# the configurational header (config.h) if enabled. The result of
# the checks' evalation is stored as 'yes' or 'no' in a variable
# enable_memcheck.
AC_DEFUN([CRB_ARG_ENABLE_MEMCHECK],
[dnl# macro-body1
 CRB_ARG_ENABLE([memcheck],
                [whether to enable memory checking],
                [enable check for memory corruption in C programs]
               )dnl# CRB_ARG_ENABLE1
 AS_IF([test "x$enable_memcheck" = xyes],
       [dnl# run-if-true1
        AC_DEFINE([MEMCHECK],
                  [1],
                  [Define to 1 to enable memory checking.])
       ]
      )dnl# AS_IF1  
]dnl# macro-body1
        )

# CRB_ARG_WITH_PC_ELISP
# -----------------------
# CRB_ARG_WITH_PC_ELISP ()
# Create an option (--with-pc-elisp) to indicate whether to
# precompile elisp code. Sets an automake conditional (AM_CONDITIONAL())
# CRB_PC_ELISP. The result of the checks' evalation is stored as 'yes' or 'no'
# in a variable with_pc_elisp.
AC_DEFUN([CRB_ARG_WITH_PC_ELISP],
[dnl# macro-body
 CRB_ARG_WITH([pc-elisp],
                [whether to enable elisp precompilation],
                [enable elisp precompilation]
               )dnl# CRB_ARG_WITH
 #crb_pc_elisp=$with_pc_elisp
 AM_CONDITIONAL([CRB_PC_ELISP], [test x$with_pc_elisp = xyes])
]dnl# macro-body
        )

# CRB_ARG_WITH_REFORMAT
# ---------------------
# CRB_ARG_WITH_REFORMAT ()
# Creates an option "--with-reformat" to control whether the reformat tool for
# source code is built or not. Additionally the option is evaluated here as
# well. If option is set, an automake conditional (AM_CONDITIONAL())
# CRB_REFORMAT is set. The result of the checks' evalation is stored as 'yes'
# or 'no' in a variable with_reformat.
AC_DEFUN([CRB_ARG_WITH_REFORMAT],
[dnl# macro-body
 AC_REQUIRE([CRB_PROG_PERL])
 AC_REQUIRE([AC_PROG_SED])
 CRB_ARG_WITH([reformat],
              [whether to build reformat tool],
              [build reformat tool for source code]
             )dnl# CRB_ARG_WITH
 AS_IF([test x$with_reformat = xyes],
       [dnl# run-if-true
        AS_IF([test -z $PERL],
              AC_MSG_ERROR([perl not found])dnl# run-if-true
             )dnl# AS_IF
        AS_IF([$PERL -e 'require 5.006;'],
              [],dnl# run-if-true
              [AC_MSG_ERROR([Perl v5.6 or better is required; Perl v5.8.2 or better is recommended. If you have several perl versions installed, select the one Automake should use using ./configure PERL=/path/to/perl])]dnl# run-if-false
             )dnl# AS_IF
       ]dnl# run-if-true
      )dnl# AS_IF   
 AM_CONDITIONAL([CRB_REFORMAT], [test x$with_reformat = xyes])
]dnl# macro-body
        )


dnl# Local variables:
dnl# eval: (add-hook 'write-file-hooks 'time-stamp)
dnl# time-stamp-start: "Last modified: "
dnl# time-stamp-format: "%:y-%02m-%02d.%02H"
dnl# time-stamp-end: "$"
dnl# End:
