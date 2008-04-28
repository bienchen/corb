/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbasic/test_argvprsr.c
 *
 *  @brief Test program for the argvprsr module.
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


#include <stdio.h>
#include <stdlib.h>
#include "errormsg.h"
#include "mprintf.h"
#include "memmgr.h"
#include "argvprsr.h"

int main(int argc,char *argv[])
{
   ArgvParser* cmdline;

   set_progname ("test_argvprsr");

   cmdline = ARGVPARSER_NEW;
   if (cmdline == NULL)
   {
      return EXIT_FAILURE;
   }

   mprintf ("%d, %s\n", argc, argv[0]);

   argvparser_delete (cmdline);

   free_progname();

   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}
