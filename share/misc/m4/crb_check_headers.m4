# Last modified: 2009-10-29.20

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

# CRB_CHECK_HEADERS(HEADER-FILE...)
# -------------------------
# Check the presence of a set of headers. If a header is not found, an error
# message is thrown and the confgure script is forced to exit imedeately. For
# each header found, a macro will be created in the confgurational header. The
# name of such a macro is HAVE_HEADER.
# This macro is just a rough rewrite of the original AC_CHECK_HEADERS macro,
# which does not provide the name of a header in an action-if-found or
# action-if-not-found statement.
#  formerly started with AH_CHECK_HEADERS([$1]) # not needed in autoconf 2.64 anymore
AC_DEFUN([CRB_CHECK_HEADERS],
[dnl# macro-body
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
