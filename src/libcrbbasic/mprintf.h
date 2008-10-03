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
 *  @file libcrbbasic/mprintf.h
 *
 *  @brief Printf wrapper to allow us to control output
 *
 *  Module:  mprintf
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-01-03
 *
 *
 *  Revision History:
 *         - 2008Jan03 bienert: created
 *
 */


#ifdef __cplusplus
extern "C" {
#endif

#ifndef MPRINTF_H
#define MPRINTF_H

#include <config.h>
#include <stdio.h>
#include <stdarg.h>

int
mvfprintf(FILE *stream, const char *format, va_list ap)
#ifdef __GNUC__
__attribute__ ((format (printf, 2, 0)))
#endif /* __GNUC__ */
;

int
mfprintf(FILE *stream, const char *format, ...)
#ifdef __GNUC__
__attribute__ ((format (printf, 2, 3)))
#endif /* __GNUC__ */
;

int
mvprintf (const char *format, va_list ap)
#ifdef __GNUC__
__attribute__ ((format (printf, 1, 0)))
#endif /* __GNUC__ */
;

int
mprintf (const char *format, ...)
#ifdef __GNUC__
__attribute__ ((format (printf, 1, 2)))
#endif /* __GNUC__ */
;

int
mvsprintf (char *str, const char *format, va_list ap)
#ifdef __GNUC__
__attribute__ ((format (printf, 2, 0)))
#endif /* __GNUC__ */
;

int
msprintf (char *str, const char *format, ...)
#ifdef __GNUC__
__attribute__ ((format (printf, 2, 3)))
#endif /* __GNUC__ */
;

int
mvsnprintf (char *str, const size_t size, const char *format, va_list ap)
#ifdef __GNUC__
__attribute__ ((format (printf, 3, 0)))
#endif /* __GNUC__ */
;

int
msnprintf (char *str, const size_t size, const char *format, ...)
#ifdef __GNUC__
__attribute__ ((format (printf, 3, 4)))
#endif /* __GNUC__ */
;

/* Now, the scarey part.  We do not want routines to do their own i/o.
 * We use the following #defines to prevent this.
 */

#define vfprintf(a)  ERROR_DONT_USE_VFPRINTF!!_Use_mvfprintf.
#define fprintf(a)   ERROR_DONT_USE_FPRINTF!!_Use_mfprintf.
#define vprintf(a)   ERROR_DONT_USE_VPRINTF!!_Use_mvprintf.
#define printf(a)    ERROR_DONT_USE_PRINTF!!_Use_mprintf.
#define vsprintf(a)  ERROR_DONT_USE_VSPRINTF!!_Use_mvsprintf.
#define sprintf(a)   ERROR_DONT_USE_SPRINTF!!_Use_msprintf.
#define snprintf(a)  ERROR_DONT_USE_SNPRINTF!!_Use_msnfprintf.
#define vsnprintf(a) ERROR_DONT_USE_VSNPRINTF!!_Use_mvsnfprintf.

#endif  /* MPRINTF_H */

#ifdef __cplusplus
}
#endif
