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
 *  @file libcrbbasic/errorh.h
 *
 *  @brief Error/ Exception handler and related functions
 *
 *  Module: errormsg
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-01-08
 *
 *
 *  Revision History:
 *         - 2008Jan08 bienert: created
 *
 */


#ifdef __cplusplus
extern "C" {
#endif

#ifndef ERRORH_H
#define ERRORH_H

#include <config.h>
#include <stdarg.h>
#include <stddef.h>

enum errormsg_retvals{
   ERR_PNAME_SET = 1,           /* program name not NULL */
   ERR_NAME_EMPTY,              /* string argument is NULL */
   ERR_PNAME_EMPTY,             /* program name is NULL */
   ERR_PS_TO_LONG               /* program name + string to large to store */
};


int
set_progname (const char*);

char*
get_progname (void);

int
add_2_progname (const char*);

int
add_name_2_progname (const char*);

size_t
get_progname_len (void);

void
free_progname (void);

int
(*get_error_msg_func(void))(const char*, int, const char*, va_list);

void
set_error_msg_func (int (*error_msgr)(const char*, int, const char*, va_list));

int call_error_msgr (const char*, int, const char*, ...)
#ifdef __GNUC__
__attribute__ ((format (printf, 3, 4)))
#endif /* __GNUC__ */
;

int
(*get_warn_msg_func(void))(const char*, int, const char*, va_list);

void
set_warn_msg_func (int (*error_msgr)(const char*, int, const char*, va_list));

int call_warn_msgr (const char*, int, const char*, ...)
#ifdef __GNUC__
__attribute__ ((format (printf, 3, 4)))
#endif /* __GNUC__ */
;

/* macros to use call_error_msgr and call_warn_msgr */
#ifdef NDEBUG
#define THROW_ERROR_MSG(...)                \
   call_error_msgr(NULL, 0, __VA_ARGS__)
#define THROW_WARN_MSG(...)                  \
   call_warn_msgr(NULL, 0, __VA_ARGS__)
#else
#define THROW_ERROR_MSG(...)                          \
   call_error_msgr(__FILE__, __LINE__, __VA_ARGS__)
#define THROW_WARN_MSG(...)                          \
   call_warn_msgr(__FILE__, __LINE__, __VA_ARGS__)
#endif

#endif /* ERRORH_H */

#ifdef __cplusplus
}
#endif
