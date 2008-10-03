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
 *  @file libcrbrna/test_alphabet.c
 *
 *  @brief Test program for the alphabet module
 *
 *  Module: alphabet
 *
 *  Library: crbrna
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-07-03
 *
 *
 *  Revision History:
 *         - 2008Jul03 bienert: created
 *
 */


#include <config.h>
#include <stdlib.h>
#include <libcrbbasic/crbbasic.h>
#include "alphabet.h"

int main(int argc __attribute__((unused)),char *argv[] __attribute__((unused)))
{
   Alphabet* test_sigma1;

   /* test creation */
   test_sigma1 = ALPHABET_NEW_PAIR ("AUGC", "augc", 4);

   if (test_sigma1 == NULL)
   {
      THROW_ERROR_MSG ("Could not create standard alphabet");
      return EXIT_FAILURE;
   }

   if (! alphabet_is_standard_rna (test_sigma1))
   {
      return EXIT_FAILURE;      
   }
   
   alphabet_delete (test_sigma1);

   test_sigma1 = ALPHABET_NEW_PAIR ("AUGCT", "augct", 5);   

   if (alphabet_is_standard_rna (test_sigma1))
   {
      THROW_ERROR_MSG ("Non-standard alphabet was identified as standard "
                       "(size check)");
      return EXIT_FAILURE;      
   }

   alphabet_delete (test_sigma1);

   test_sigma1 = ALPHABET_NEW_PAIR ("ATGC", "atgc", 4);

   if (alphabet_is_standard_rna (test_sigma1))
   {
      THROW_ERROR_MSG ("Non-standard alphabet was identified as standard "
                       "(symbol check)");
      return EXIT_FAILURE;
   }

   /* delete alphabet */
   alphabet_delete (test_sigma1);

   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}
