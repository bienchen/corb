#
# Copyright (C) 2008 Stefan Bienert
# Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
#
# See COPYING file in the top level directory of this tree for licence.
#
# Last modified: 2008-04-22.18
#


# CRB_CHECK_HEADERS(HEADER-FILE...)
# -------------------------
# Check the presence of a set of headers. If a header is not found, an error
# message is thrown and the confgure script is forced to exit imedeately. For
# each header found, a macro will be created in the confgurational header. he
# name of such a macro is HAVE_HEADER.
# This macro is just a rough rewrite of the original AC_CHECK_HEADERS macro,
# which does not provide the name of a header in an action-if-found or
# action-if-not-found statement.
AC_DEFUN([CRB_CHECK_HEADERS],
[dnl# macro-body
 AH_CHECK_HEADERS([$1])dnl
 for __crb_header in $1
 do
 AC_CHECK_HEADER($__crb_header,
	 	 [AC_DEFINE_UNQUOTED(AS_TR_CPP(HAVE_$__crb_header))],
		 [AC_MSG_ERROR([Mandatory header not found: $__crb_header])],
		 [AC_INCLUDES_DEFAULT]
                )dnl# AC_CHECK_HEADER
done
]dnl# macro-body
        )


dnl# Local variables:
dnl# eval: (add-hook 'write-file-hooks 'time-stamp)
dnl# time-stamp-start: "Last modified: "
dnl# time-stamp-format: "%:y-%02m-%02d.%02H"
dnl# time-stamp-end: "$"
dnl# End:
