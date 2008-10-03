/*
 * Copyright (C) 2008 Stefan Bienert
 *
 * This file is part of CoRB.
 *
 * CoRB is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CoRB is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CoRB.  If not, see <http://www.gnu.org/licenses/>.
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
