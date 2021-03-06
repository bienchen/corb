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
 *  @file libcrbbasic/gfile.h
 *
 *  @brief Generic file handling
 *
 *  Module: gfile
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-11-01
 *
 *
 *  Revision History:
 *         - 2008Nov01 bienert: created
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef GFILE_H
#define GFILE_H

#include <stdio.h>
#include "str.h"

typedef enum {
   GFILE_VOID = 0,
   GFILE_UNCOMPRESSED
} GFileType;

enum {
   GFILE_UNKNOWN_TYPE = 1,
   GFILE_READ_ERROR,
   GFILE_MEM_ERROR,
   GFILE_REWIND_ERROR
};

enum {
   GFILE_TR_FROM = 0,
   GFILE_TR_TO,
   GFILE_TR_N
};

/* line feed character */
#define CRB_LF '\n'
#define CRB_TAB '\t'
#define CRB_COM '#'

typedef struct GFile GFile;


GFileType
gfile_get_type (const char*, unsigned long);

const char*
gfile_get_type_str (const char*, unsigned long);

unsigned long
gfile_ext_from_list (const char*, const unsigned long, const char**,
                     const unsigned long);

Str*
gfile_get_path (const GFile*);

GFile*
gfile_open (const char*, const unsigned long, const GFileType, const char*,
            const char*, const int);

#define GFILE_OPEN(PATH, LEN, TYPE, MODE) \
   gfile_open (PATH, LEN, TYPE, MODE, __FILE__, __LINE__)

int
gfile_close (GFile*);

int
gfile_rewind (GFile*);

unsigned long
gfile_getline_verbatim (int*, char**, size_t*, GFile*);

unsigned long
gfile_getline_tab (int*, char**, size_t*, GFile*);

unsigned long
gfile_getline (int*, char**, size_t*, const GFile*);

int
gfile_printf (GFile*, const char*, ...)
#ifdef __GNUC__
__attribute__ ((format (printf, 2, 3)))
#endif /* __GNUC__ */
;

#endif /* GFILE_H */

#ifdef __cplusplus
}
#endif
