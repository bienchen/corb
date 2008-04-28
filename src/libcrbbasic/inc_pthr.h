/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbasic/inc_pthr.h
 *
 *  @brief Wrapper to include pthread.h properly.
 *
 *  Module:  *** Module name ***
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author 
 *
 *  @date 2008-01-31
 *
 *
 *  Revision History:
 *         - 2008Jan31 bienert: created
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef INC_PTHR_H
#define INC_PTHR_H

#include <config.h>

#ifdef HAVE_PTHREAD
#include <sys/types.h>
#include <pthread.h>
#endif


#endif /* INC_PTHR_H */

#ifdef __cplusplus
}
#endif
