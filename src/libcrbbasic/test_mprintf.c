/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
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
