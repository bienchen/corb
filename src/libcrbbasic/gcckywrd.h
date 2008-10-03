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
 *  @file libcrbbasic/gcckywrd.h
 *
 *  @brief Undefining GCC keywords if not compiling with gcc
 *
 *  Module:  none
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-01-25
 *
 *
 *  Revision History:
 *         - 2008Jan25 bienert: created
 *
 */


#ifdef __cplusplus
extern "C" {
#endif

#ifndef GCCKYWRD_H
#define GCCKYWRD_H

#ifndef __GNUC__
#define __inline__ 
#endif
   
#endif /* GCCKYWRD_H */

#ifdef __cplusplus
}
#endif
