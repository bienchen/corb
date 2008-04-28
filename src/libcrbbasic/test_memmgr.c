/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbasic/test_memmgr.c
 *
 *  @brief Test program for the memmgr module.
 *
 *  Module: memmgr
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-01-18
 *
 *
 *  Revision History:
 *         - 2008Jan18 bienert: created
 *
 */


#include <stdlib.h>
#include "mprintf.h"
#include "memmgr.h"

int main(void)
{
   char** string;
   unsigned long i;
   unsigned long sample = 100000;

   /* test allocating */
   string = XMALLOC (sample * sizeof (char*));

   for (i=0; i < sample; i++)
   {
      string[i] = XMALLOC (10 * sizeof (char));
      if (string == NULL)
      {
         mfprintf (stderr, "Failed to use xmalloc()\n");
         return EXIT_FAILURE;
      }
   }

   for (i=0; i < sample; i++)
   {
      XFREE (string[i]);
   }

   XFREE (string);

   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}
