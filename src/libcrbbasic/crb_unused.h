/*
 * Copyright (C) 2010 Stefan Bienert
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
 *  @file libcrbbasic/crb_unused.h
 *
 *  @brief CRB_UNUSED_ARG macro for functions overwritten as empty.
 *
 *  Module: none
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2010-07-14
 *
 *
 *  Revision History:
 *         - 2010Jul14 bienert: created
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef CRBUNUSED_H
#define CRBUNUSED_H

#include <config.h>

#ifdef CRB_DEF_UNUSED
#define CRB_UNUSED(X) (void) X
#else
#define CRB_UNUSED(X)
#endif
  
   
#endif /* CRBUNUSED_H */

#ifdef __cplusplus
}
#endif
