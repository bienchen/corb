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
   unsigned long width = 10;
   size_t dim[] = {100000, 10};
   char test_text[] = "Hello World";

   /* test allocating */
   string = XMALLOC (sample * sizeof (*string));
   if (string == NULL)
   {
      mfprintf (stderr, "Failed to use xmalloc()\n");
      return EXIT_FAILURE;
   }

   for (i=0; i < sample; i++)
   {
      string[i] = XMALLOC (width * sizeof (**string));
      if (string[i] == NULL)
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

   /* check 2D allocation */
   string = (char**) XMALLOC_2D(sample, width, sizeof (**string));
   if (string == NULL)
   {
      mfprintf (stderr, "Failed to use xmalloc_2d()\n");
      return EXIT_FAILURE;
   }

   for (i = 0; i < sample; i++)
   {
      string[i][0] = 'H';
   }

   XFREE_2D ((void**) string);

   /* doing 2D allocation with rnd function */
   string = (char**) XMALLOC_RND (sizeof (**string),
                                  sizeof (dim)/sizeof (*dim),
                                  dim);
   if (string == NULL)
   {
      mfprintf (stderr, "Failed to use xmalloc_rnd()\n");
      return EXIT_FAILURE;
   }

   for (i = 0; i < sample; i++)
   {
      string[i][0] = 'H';
   }

   XFREE_ND (sizeof (dim)/sizeof (*dim), (void**) string);

   /* doing 2D allocation using nd function */
   string = (char**) XMALLOC_ND (sizeof (**string), 2, dim[0], dim[1]);
   if (string == NULL)
   {
      mfprintf (stderr, "Failed to use xmalloc_nd()\n");
      return EXIT_FAILURE;
   }

   for (i = 0; i < sample; i++)
   {
      string[i][0] = 'H';
   }

   XFREE_ND (2, (void**) string);

   /* trying to allocate an 2D array to store text (char*) */
   string = (char**) XMALLOC_ND (sizeof (*string), 1, dim[0]);
   if (string == NULL)
   {
      mfprintf (stderr, "Failed to use xmalloc_nd()\n");
      return EXIT_FAILURE;
   }

   for (i = 0; i < sample; i++)
   {
      string[i] = test_text;
   }

   /* we may only free the pointer we allocated! */
   XFREE ((void*) string);

   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}
