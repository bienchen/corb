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
 *  @file libcrbbasic/mprintf_dummy.h
 *
 *  @brief Dummy header file for mprintf source.
 *
 *  Module:  *** Module name ***
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAnalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-04-08
 *
 *
 *  Revision History:
 *         - 2008Apr08 bienert: created
 *
 */

/*
  This header just provides declarations of the external functions of mprintf.
  That means only mprintf.c includes this file to work with some very pedantic
  compilers. These compilers require to see an external declaration for a
  function which is not set to be static. But mprintf can not include the
  usual mprintf header since it disables all the standard printing functions
  at its the end...
*/

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MPRINTF_DUMMY_H
#define MPRINTF_DUMMY_H

extern int
mvfprintf(FILE*, const char*, va_list)
#ifdef __GNUC__
__attribute__ ((format (printf, 2, 0)))
#endif /* __GNUC__ */
;

extern int
mfprintf(FILE*, const char*, ...)
#ifdef __GNUC__
__attribute__ ((format (printf, 2, 3)))
#endif /* __GNUC__ */
;

extern int
mvprintf (const char*, va_list)
#ifdef __GNUC__
__attribute__ ((format (printf, 1, 0)))
#endif /* __GNUC__ */
;

extern int
mprintf (const char*, ...)
#ifdef __GNUC__
__attribute__ ((format (printf, 1, 2)))
#endif /* __GNUC__ */
;

extern int
mvsprintf (char*, const char*, va_list)
#ifdef __GNUC__
__attribute__ ((format (printf, 2, 0)))
#endif /* __GNUC__ */
;

extern int
msprintf (char*, const char*, ...)
#ifdef __GNUC__
__attribute__ ((format (printf, 2, 3)))
#endif /* __GNUC__ */
;

extern int
mvsnprintf (char*, const size_t, const char*, va_list)
#ifdef __GNUC__
__attribute__ ((format (printf, 3, 0)))
#endif /* __GNUC__ */
;

extern int
msnprintf (char*, const size_t, const char*, ...)
#ifdef __GNUC__
__attribute__ ((format (printf, 3, 4)))
#endif /* __GNUC__ */
;
   
   
#endif /* MPRINTF_DUMMY_H */

#ifdef __cplusplus
}
#endif
