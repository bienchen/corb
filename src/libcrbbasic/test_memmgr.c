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
