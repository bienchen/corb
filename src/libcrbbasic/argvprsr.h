/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
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
