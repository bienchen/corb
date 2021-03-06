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
#include "mprintf.h"
#include "gfile.h"


int main(int argc __attribute__((unused)),char *argv[] __attribute__((unused)))
{
   GFile* file;
   char c_file[] = "Makefile";
   int ret_val = 0;
   int error;
   size_t sizeof_buf = 0;
   char* test = NULL;

   /* open file */
   THROW_WARN_MSG ("Trying to OPEN   file \"%s\".", c_file);
   file = GFILE_OPEN(c_file, strlen (c_file), GFILE_VOID, "r");
   if (file == NULL)
   {
      THROW_ERROR_MSG ("Could not open \"%s\".", c_file);
      return EXIT_FAILURE;
   }

   THROW_WARN_MSG ("Trying to READ   file \"%s\".", c_file);
   while (gfile_getline_verbatim (&error, &test, &sizeof_buf, file) > 0)
   {
      mprintf ("%s\n", test);
   }

   THROW_WARN_MSG ("Trying to REWIND file \"%s\".", c_file);
   ret_val = gfile_rewind (file);
   if (ret_val)
   {
      THROW_ERROR_MSG ("Could not rewind \"%s\".", c_file);      
   }

   THROW_WARN_MSG ("Trying to READ   file \"%s\" with tab translation.",
                   c_file);
   while (gfile_getline_tab (&error, &test, &sizeof_buf, file))
   {
      mprintf ("%s\n", test);
   }

   THROW_WARN_MSG ("Trying to REWIND file \"%s\".", c_file);
   ret_val = gfile_rewind (file);
   if (ret_val)
   {
      THROW_ERROR_MSG ("Could not rewind \"%s\".", c_file);      
   }

   THROW_WARN_MSG ("Trying to READ   file \"%s\" with checks for comments.",
                   c_file);
   while (gfile_getline (&error, &test, &sizeof_buf, file))
   {
      mprintf ("%s\n", test);
   }

   THROW_WARN_MSG ("Trying to CLOSE file \"%s\".", c_file);
   ret_val = gfile_close (file);
   if (ret_val == EOF)
   {
      THROW_ERROR_MSG ("Could not close \"%s\".", c_file);      
   }

   /* create compressed files */

   /* read compr. files */

   /* delete files */

   /* try to open again */
   THROW_WARN_MSG ("Trying to OPEN   file \"%s\" for writing.", "test.test");
   file = GFILE_OPEN("test.test", strlen ("test.test"), GFILE_VOID, "w");
   if (file == NULL)
   {
      THROW_ERROR_MSG ("Could not open \"%s\".", "test.test");
      return EXIT_FAILURE;
   }

   if (gfile_printf (file, "Hallo") < 0)
   {
      return EXIT_FAILURE;
   }
   if (gfile_printf (file, " Welt") < 0)
   {
      return EXIT_FAILURE;
   }
   if (gfile_printf (file, "!\n") < 0)
   {
      return EXIT_FAILURE;
   }


   THROW_WARN_MSG ("Trying to CLOSE file \"%s\".", "test.test");
   ret_val = gfile_close (file);
   if (ret_val == EOF)
   {
      THROW_ERROR_MSG ("Could not close \"%s\".", "test.test");      
   }

   XFREE (test);
   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}
