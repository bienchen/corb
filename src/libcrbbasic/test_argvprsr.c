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
