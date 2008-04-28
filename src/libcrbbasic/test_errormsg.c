/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbasic/test_errorc.c
 *
 *  @brief Test program for the errormsg module.
 *
 *  Module: errormsg
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-01-10
 *
 *
 *  Revision History:
 *         - 2008Jan10 bienert: created
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "inc_strg.h"
#include "mprintf.h"
#include "memmgr.h"
#include "errormsg.h"

static int
alt_msgr(const char* file, int line, const char *format, va_list ap)
{
   if (file != NULL)
   {
      mfprintf (stderr, "I see a file...%s and a line, it's %d! ", file, line);
   }

   if (format != NULL)
   {
      mfprintf (stderr, "Uh, and there's a message for you: ");
      mvfprintf (stderr, format, ap);
   }

   mfprintf (stderr, "\n");

   return 1;
}

int main(int argc,char *argv[])
{
   int retval;
   unsigned long len;

   if (argc > 1)
   {
      mfprintf (stderr, "%s called with an argument.\n", argv[0]);      
      return EXIT_FAILURE;
   }

   /* test program name */
   retval = set_progname (argv[0]);
   if (retval != 0)
   {
      mfprintf (stderr, "Failed to init messenger: %s.\n", argv[0]);      
      return EXIT_FAILURE;      
   }

   retval = add_name_2_progname (argv[0]);
   if (retval != 0)
   {
      mfprintf (stderr, "Failed to add name to program name: %s.\n", argv[0]);
      return EXIT_FAILURE;      
   }

   len = mprintf ("%s\n", get_progname());
   if (len != ((strlen (argv[0]) * 2) + 2))
   {
      mfprintf (stderr, "Failed to get program name: %s.\n", argv[0]);
      return EXIT_FAILURE;      
   }

   len = get_progname_len();
   if (len != ((strlen (argv[0]) * 2) + 1))
   {
      mfprintf (stderr, "Failed to get program name length: %s.\n", argv[0]);
      return EXIT_FAILURE;      
   }

   if (THROW_ERROR_MSG ("%s is %d test. Don%ct worry!", "This", 1, '\'') < 0)
   {
      mfprintf (stderr, "Failed to print error message: %s\n", argv[0]);
      return EXIT_FAILURE;
   }

   set_error_msg_func (alt_msgr);
   THROW_ERROR_MSG ("Hello chap!");

   if (THROW_WARN_MSG ("Now testing %d %s.", 1, "Warning") < 0)
   {
      mfprintf (stderr, "Failed to print error message: %s\n", argv[0]);
      return EXIT_FAILURE;
   }


   set_warn_msg_func (alt_msgr);
   THROW_WARN_MSG ("Hello bloke!");

   free_progname();

   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}

