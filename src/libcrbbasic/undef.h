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
 *  @file libcrbbasic/undef.h
 *
 *  @brief Undefiend values for basic datatypes
 *
 *  Module: none
 *
 *  Library: libcrbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-03-04
 *
 *
 *  Revision History:
 *         - 2008Mar04 bienert: created
 *
 */


#include <limits.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef UNDEF_H
#define UNDEF_H

#define ULONG_UNDEF ULONG_MAX
#define CHAR_UNDEF  CHAR_MAX
#define UCHAR_UNDEF UCHAR_MAX


#endif /* UNDEF_H */

#ifdef __cplusplus
}
#endif

 
