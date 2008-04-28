/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbasic/inc_bool.h
 *
 *  @brief Wrapper to include stdbool.h properly.
 *
 *  Module:  *** Module name ***
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author 
 *
 *  @date 2008-02-22
 *
 *
 *  Revision History:
 *         - 2008Feb22 bienert: created
 *
 */



#ifdef __cplusplus
extern "C" {
#endif

#ifndef INC_BOOL_H
#define INC_BOOL_H

#include <config.h>

#ifdef HAVE_STDBOOL_H
# include <stdbool.h>
#else
# ifndef HAVE__BOOL
#  ifdef __cplusplus
   typedef bool _Bool;
#  else
#   define _Bool signed char
#  endif
# endif
# define bool _Bool
# define false 0
# define true 1
# define __bool_true_false_are_defined 1
#endif


#endif /* INC_BOOL_H */

#ifdef __cplusplus
}
#endif
