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
 *  @file libcrbbasic/test_genarray.c
 *
 *  @brief Test genarray implementation
 *
 *  Module: none
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-06-26
 *
 *
 *  Revision History:
 *         - 2008Jun26 bienert: created
 *
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "memmgr.h"
#include "errormsg.h"
#include "genarray.h"

typedef struct { int n; } tint;

ARRAY_CREATE_CLASS(tint);

int main(int argc __attribute__((unused)), char *argv[] __attribute__((unused)))
{
   Arraytint test_array;
   int def_value = 1;
   unsigned long init_size = 2;
   tint u = { def_value };
   unsigned long current = 0;
   unsigned long i;

   set_progname (argv[0]);

   ARRAY_INIT (test_array, init_size, tint);

   if (ARRAY_IS_NULL(test_array))
   {
      THROW_ERROR_MSG ("Could not create test array");
      return EXIT_FAILURE;
   }

   ARRAY_PUSH (test_array, u, tint,
               {
                  THROW_ERROR_MSG ("Could not push to test array");
                  return EXIT_FAILURE;                  
               });

   if (ARRAY_SIZE (test_array) != init_size)
   {
      THROW_ERROR_MSG ("Unexpected change in size on push. Expected: %lu, "
                       "found: %lu", init_size, ARRAY_SIZE (test_array));
      return EXIT_FAILURE;      
   }

   current = ARRAY_CURRENT (test_array);
   if (current == 0)
   {
      THROW_ERROR_MSG ("Indexing after push failed");
      return EXIT_FAILURE;
   }

   current--;

   ARRAY_GET (u, test_array, current);
   if (u.n != def_value)
   {
      THROW_ERROR_MSG ("Wrong value stored on push, expected: %d, found: %d",
                       def_value, u.n);
      return EXIT_FAILURE;    
   }

   u.n = 2 * def_value;
   ARRAY_SET (test_array, u, current);
   ARRAY_POP (u, test_array);
   if (u.n != (2 * def_value))
   {
      THROW_ERROR_MSG ("Wrong value stored on pop, expected: %d, found: %d",
                       (2 * def_value) , u.n);
      return EXIT_FAILURE;      
   }   

   if (ARRAY_CURRENT (test_array) != current)
   {
      THROW_ERROR_MSG ("Index change after pop failed, expected: %lu, found: "
                       "%lu", current, ARRAY_CURRENT (test_array));
      return EXIT_FAILURE;      
   }

   init_size = ARRAY_SIZE (test_array);
   for (i = 0; i <= init_size; i++)
   {
      ARRAY_PUSH (test_array, u, tint,
                  {
                     THROW_ERROR_MSG ("Could not reallocate test array");
                     return EXIT_FAILURE;                     
                  });
   }

   if (ARRAY_SIZE (test_array) == init_size)
   {
      THROW_ERROR_MSG ("Index change after multiple push failed, expected "
                       "value > %lu, found: %lu", init_size,
                       ARRAY_SIZE (test_array));
      return EXIT_FAILURE;      
   }

   ARRAY_DELETE (test_array);

   if (ARRAY_SIZE (test_array) != 0)
   {
      THROW_ERROR_MSG ("Array does not have size 0 after deletion, found "
                       "size: %lu", ARRAY_SIZE (test_array));
      return EXIT_FAILURE;       
   }

   free_progname();

   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}
