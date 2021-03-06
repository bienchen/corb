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
 *  @file libcrbbasic/test_str.c
 *
 *  @brief Test program for the str module
 *
 *  Module: str
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-02-05
 *
 *
 *  Revision History:
 *         - 2008Feb05 bienert: created
 *
 */


#include <config.h>
#include <stdlib.h>
#include <assert.h>
#include "inc_strg.h"
#include "memmgr.h"
#include "errormsg.h"
#include "mprintf.h"
#include "str.h"


int main(int argc __attribute__((unused)),char *argv[] __attribute__((unused)))
{
   Str* test_str;
   Str* test_str2;

   /* test creation */
   test_str = STR_NEW_CHAR ('.', 1000);
   if (test_str == NULL)
   {
      THROW_ERROR_MSG ("Could not create first string:");
      return EXIT_FAILURE;
   }   

   str_delete (test_str);
   test_str = STR_NEW_CSTR ("HELLO WORLD");
   if (test_str == NULL)
   {
      THROW_ERROR_MSG ("Could not create first string:");
      return EXIT_FAILURE;
   }
   
   test_str2 = STR_NEW_STR (test_str);
   if (test_str2 == NULL)
   {
      THROW_ERROR_MSG ("Could not create second string:");
      return EXIT_FAILURE;
   }

   if (str_set (test_str2, "Hello") != 0)
   {
      THROW_ERROR_MSG ("Could not set strings.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE;
   }

   if (str_cpy (test_str, test_str2) != 0)
   {
      THROW_ERROR_MSG ("Could not copy strings.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE;
   }

   if (str_at (test_str, 2, 'A') != 0)
   {
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE;     
   }

   if (str_set (test_str, "Hello") != 0)
   {
      THROW_ERROR_MSG ("Could not set string.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE;     
   }

   if (str_append_cstr (test_str, " Macao") != 0)
   {
      THROW_ERROR_MSG ("Could not append strings.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE;     
   }

   mprintf ("%s length: %lu size: %zu\n", str_get (test_str),
            str_length (test_str), str_capacity (test_str));

   str_resize (test_str, 6, '\0');

   if (str_empty (test_str) == true)
   {
      THROW_ERROR_MSG ("String not supposed to be empty.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE; 
   }

   str_clear (test_str);
   if (str_empty (test_str) == false)
   {
      THROW_ERROR_MSG ("Clearing string did not work.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE; 
   }

   str_set (test_str, "GCATGACATAGAGAGGAGAGAGTAGACGCTACG");
   str_set (test_str2, "AGAGAGT");
   if (str_find_str (test_str, test_str2) != 17)
   {
      THROW_ERROR_MSG ("String search did not work.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE; 
   }

   str_set (test_str,  "GCATCGCAGATGAGAGAGGAGAGATACAGTACG");
   str_set (test_str2, "TGAGAGA");
   if (str_rfind_str (test_str, test_str2) != 11)
   {
      THROW_ERROR_MSG ("String search did not work.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE; 
   }


   str_set (test_str2, "AUGC");
   if (str_find_first_of_str (test_str, test_str2, 2) != 2)
   {
      THROW_ERROR_MSG ("Element search did not work.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE; 
   }

   if (str_find_last_of_str (test_str, test_str2, 10) != 10)
   {
      THROW_ERROR_MSG ("Element search did not work.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE; 
   }
   /* find last not of? */
   if (str_find_first_not_of_str (test_str, test_str2, 2) != 4)
   {
      THROW_ERROR_MSG ("Element search did not work.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE; 
   }

   /* compare strings */
   str_set (test_str, "AAAA");
   str_set (test_str2, "AAAA");
   if (str_compare_str (test_str, test_str2) != 0)
   {
      THROW_ERROR_MSG ("Comparsion of equal strings did not work.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE;     
   }

   str_set (test_str2, "BBBB");
   if (str_compare_str (test_str, test_str2) != -1)
   {
      THROW_ERROR_MSG ("Comparsion of strings did not work.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE;     
   }

   if (str_compare_str (test_str2, test_str) != 1)
   {
      THROW_ERROR_MSG ("Comparsion of strings did not work.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE;     
   }

   /* compare substrings */
   if (str_compare_substr (test_str, 0, 2, test_str2, 0, 2) != -1)
   {
      THROW_ERROR_MSG ("Comparsion of substrings did not work.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE;     
   }   

   str_set (test_str2, "BAAB");
   if (str_compare_substr (test_str, 2, 2, test_str2, 2, 2) != 0)
   {
      THROW_ERROR_MSG ("Comparsion of substrings did not work.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE;     
   } 

   /* out of bounds test */
   str_set (test_str2, "BAAA");
   if (str_compare_substr (test_str, 2, 20, test_str2, 2, 20) != 0)
   {
      THROW_ERROR_MSG ("Comparsion of substrings did not work.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE;     
   }    

   if (str_compare_substr (test_str, 20, 2, test_str2, 2, 20) != -1)
   {
      THROW_ERROR_MSG ("Comparsion of substrings did not work.");
      str_delete (test_str);
      str_delete (test_str2);
      return EXIT_FAILURE;     
   } 

   /* delete strings */
   str_delete (test_str);
   str_delete (test_str2);

   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}
