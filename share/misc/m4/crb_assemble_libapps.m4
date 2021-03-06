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

AC_DEFUN([CRB_LIBAPPSSOURCE], [])


AC_DEFUN([_CRB_LIBAPPSOBJ],
[AS_LITERAL_IF([$1],
               [CRB_LIBAPPSSOURCE([$1.c])],
               [$2])dnl
case " $LIBCRBAPPS@&t@OBJS " in
  *" $1.$ac_objext "* ) ;;
  *) AC_SUBST([LIBCRBAPPS@&t@OBJS], ["$LIBCRBAPPS@&t@OBJS $1.$ac_objext"]) ;;
esac
])

AC_DEFUN([CRB_LIBAPPSOBJ],            
[_CRB_LIBAPPSOBJ([$1],
            [AC_DIAGNOSE(syntax,                                                
                         [$0($1): you should use literals])])dnl
])

 
dnl# Local variables:
dnl# eval: (add-hook 'write-file-hooks 'time-stamp)
dnl# time-stamp-start: "Last modified: "
dnl# time-stamp-format: "%:y-%02m-%02d.%02H"
dnl# time-stamp-end: "$"
dnl# End:
