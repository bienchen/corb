# Copyright (C) 2008 Stefan Bienert
# Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
#
# See COPYING file in the top level directory of this tree for licence.
#
# Last modified: 2008-04-20.14
#


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
