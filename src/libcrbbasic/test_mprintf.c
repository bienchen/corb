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
 *  @file libcrbbasic/test_mprintf.c
 *
 *  @brief Test program for the mprintf module
 *
 *  Module: mprintf
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


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "inc_strg.h"
#include "mprintf.h"

int main(int argc, char *argv[])
{
   size_t written_chars;
   size_t test_msg_len;
   char* buf;
   const char* test_msg = "Hello Dolly!\n";

   if (argc > 1)
   {
      mfprintf (stderr, "%s called with an argument.\n", argv[0]);      
      return EXIT_FAILURE;
   }

   test_msg_len = strlen (test_msg);

   /* test mfprintf */
   written_chars = mfprintf(stderr, "%s", test_msg);
   if (written_chars != test_msg_len)
   {
      mfprintf (stderr, "Testing \"mfprintf\" failed, test message length: "
                        "%zd, written: %zd\n", test_msg_len, written_chars);
      return EXIT_FAILURE;
   }

   /* test mprintf */
   written_chars = mprintf("%s", test_msg);
   if (written_chars != test_msg_len)
   {
      mfprintf (stderr, "Testing \"mprintf\" failed, test message length: "
                        "%zd, written: %zd\n", test_msg_len, written_chars);
      return EXIT_FAILURE;
   }

   /* allocate buffer for sprintf's */
   buf = calloc ((test_msg_len + 1), sizeof(char));

   /* test msprintf */
   written_chars = msprintf (buf, "%s", test_msg);
   if (written_chars != test_msg_len)
   {
      mfprintf (stderr, "Testing \"msprintf\" failed, test message length: "
                        "%zd, written: %zd\n", test_msg_len, written_chars);
      return EXIT_FAILURE;
   }   

   /* test snprintf if necessary */
/* #ifndef HAVE_SNPRINTF */
   written_chars = msnprintf (buf, (test_msg_len + 1), "%s", test_msg);
   if (written_chars != test_msg_len)
   {
      mfprintf (stderr, "Testing \"msnprintf\" failed, test message length: "
                        "%zd, written: %zd\n", test_msg_len, written_chars);
      return EXIT_FAILURE;
   }
/* #endif */

   free (buf);

   return EXIT_SUCCESS;
}
