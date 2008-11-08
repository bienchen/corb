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
 *  @file libcrbbasic/test_gfile.c
 *
 *  @brief Testing the gfile module
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


#include <config.h>
#include <stdlib.h>
#include <assert.h>
#include "inc_strg.h"
#include "memmgr.h"
#include "errormsg.h"
/*#include "mprintf.h"*/
#include "gfile.h"
 

int main(int argc __attribute__((unused)),char *argv[] __attribute__((unused)))
{
   GFile* file;
   char c_file[] = "test_gfile.c";
   int ret_val = 0;

   /* open file */
   file = GFILE_OPEN(c_file, strlen (c_file), GFILE_VOID, "r");
   if (file == NULL)
   {
      THROW_ERROR_MSG ("Could not open \"%s\".", c_file);
      return EXIT_FAILURE;
   }   

   ret_val = gfile_close (file);
   if (ret_val == EOF)
   {
      THROW_ERROR_MSG ("Could not close \"%s\".", c_file);      
   }

   /* create compressed files */

   /* read compr. files */

   /* delete files */

   /* try to open again */


   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}
