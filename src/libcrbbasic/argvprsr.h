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
 *  @file libcrbbasic/argvprsr.h
 *
 *  @brief Parsing the command line (argv)
 *
 *  Module: argvprsr
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-02-04
 *
 *
 *  Revision History:
 *         - 2008Feb04 bienert: created
 *
 */


#ifdef __cplusplus
extern "C" {
#endif

#ifndef ARGVPRSR_H
#define ARGVPRSR_H


/* error numbers */
enum argvprsr_retvals{
   ERR_AP_ALLOC = 1 /* Memory (re)allocation failure */
};


typedef struct ArgvParser ArgvParser;


extern ArgvParser*
argvparser_new (const char*, const int);

#define ARGVPARSER_NEW argvparser_new(__FILE__, __LINE__)

extern void
argvparser_delete (ArgvParser*);

extern char*
argvparser_get_maintainername (ArgvParser*);

extern int
argvparser_set_maintainername (char*, ArgvParser*);

#endif /* ARGVPRSR_H */

#ifdef __cplusplus
}
#endif
