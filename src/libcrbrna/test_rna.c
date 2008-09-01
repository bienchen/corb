/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbrna/test_rna.c
 *
 *  @brief Test program for the rna module
 *
 *  Module: rna
 *
 *  Library: crbrna
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-08-19
 *
 *
 *  Revision History:
 *         - 2008Aug19 bienert: created
 *
 */


#include <config.h>
#include <stdlib.h>
#include <libcrbbasic/crbbasic.h>
#include "rna.h"

int main(int argc __attribute__((unused)),char *argv[] __attribute__((unused)))
{
   char* test_string = "STRUCTURE";
   Alphabet* sigma;
   int error = 0;
   unsigned long l, r;
   Rna* test_obj = RNA_NEW;

   if (test_obj == NULL)
   {
      return EXIT_FAILURE;
   }

   /* test with malformed structure */
   THROW_WARN_MSG ("Testing recognition of malformed structure. You should "
                   "see an error message in the following.");
   error = RNA_INIT_PAIRLIST_VIENNA(test_string, strlen (test_string),
                                    test_obj);
   if (error != ERR_RNA_VIENNA_FORMAT)
   {
      THROW_ERROR_MSG ("Producing error by using \'%s\' as a structure "
                       "string failed!", test_string);
      return EXIT_FAILURE;
   }

   /* test with unbalanced structure 1 */
   THROW_WARN_MSG ("Testing recognition of unbalanced structure. You "
                   "should see an error message in the following.");
   test_string = "........)";
   error = RNA_INIT_PAIRLIST_VIENNA(test_string, strlen (test_string),
                                    test_obj);
   if (error != ERR_RNA_VIENNA_MMC)
   {
      THROW_ERROR_MSG ("Producing error by using \'%s\' as a structure "
                       "string failed!", test_string);
      return EXIT_FAILURE;
   }

   /* test with unbalanced structure 2 */
   THROW_WARN_MSG ("Testing recognition of unbalanced structure. You should "
                   "see an error message in the following.");
   test_string = "(........";
   error = RNA_INIT_PAIRLIST_VIENNA(test_string, strlen (test_string),
                                    test_obj);
   if (error != ERR_RNA_VIENNA_MMO)
   {
      THROW_ERROR_MSG ("Producing error by using \'%s\' as a structure "
                       "string failed!", test_string);
      return EXIT_FAILURE;
   }

   /* set correct structure */
   test_string = "((.....))";
   error = RNA_INIT_PAIRLIST_VIENNA(test_string, strlen (test_string),
                                    test_obj);
   if (error)
   {
      THROW_ERROR_MSG ("Unintentional error while testing with structure "
                       "\'%s\'!", test_string);
      return EXIT_FAILURE;      
   }

   /* testing accessing pairs */
   l = rna_base_pairs_with (0, test_obj);
   if (l != strlen (test_string))
   {
      THROW_ERROR_MSG ("Wrong index for pairing partner found. Checked pos. "
                       "%d, expected partner index %lu, got %lu. Structure "
                       "was: \'%s\'!",
                       0, (unsigned long) strlen (test_string), l, test_string);
      return EXIT_FAILURE;       
   }
   r = rna_base_pairs_with (strlen (test_string) - 1, test_obj);
   if (r != 1)
   {
      THROW_ERROR_MSG ("Wrong index for pairing partner found. Checked pos. "
                       "%d, expected partner index %d, got %lu. Structure "
                       "was: \'%s\'!",
                       (strlen (test_string) - 1), 1, r, test_string);
      return EXIT_FAILURE;       
   }

   
   /* testing sequence: Initialisation */
   test_string = "AUGCAUGCA";
   error = RNA_INIT_SEQUENCE(strlen (test_string), test_obj);
   if (error)
   {
      THROW_ERROR_MSG ("Unintentional error while trying to init the sequence "
                       "component of an Rna object!");
      return EXIT_FAILURE;      
   }

   /* set bases in sequence */
   for (l = 0; l < strlen (test_string); l++)
   {
      rna_set_sequence_base (test_string[l], l, test_obj);
   }

   /* check sequence stored */
   if (strcmp (rna_get_sequence (test_obj), test_string) != 0)
   {
      THROW_ERROR_MSG ("Unintentional error: Stored sequence and sequence to "
                       "be stored differ: src = \'%s\' dest = \'%s\'",
                       test_string, rna_get_sequence (test_obj));
      return EXIT_FAILURE;
   }

   /* try failed transformation */
   sigma = ALPHABET_NEW_SINGLE (RNA_ALPHABET, strlen(RNA_ALPHABET)/2);
   
   THROW_WARN_MSG ("Testing recognition of invalid sequences. You should "
                   "see an error message in the following.");
   l = 0;
   rna_set_sequence_base ('\n', l, test_obj);
   error = rna_transform_sequence_2_no (sigma, test_obj);
   if (error != ERR_RNA_NO_BASE)
   {
      THROW_ERROR_MSG ("Invalid base in sequence not recognised, was \'\\n\' "
                       "at position %lu.", l);
      return EXIT_FAILURE;
   }

   rna_set_sequence_base (test_string[l], l, test_obj);   
   if (strcmp (rna_get_sequence (test_obj), test_string) != 0)
   {
      THROW_ERROR_MSG ("Unintentional error: Sequence is changed after failed "
                       "transformation: src = \'%s\' dest = \'%s\'",
                       test_string, rna_get_sequence (test_obj));
      return EXIT_FAILURE;
   }

   /* test true transformation */
   error = rna_transform_sequence_2_no (sigma, test_obj);
   if (error)
   {
      THROW_ERROR_MSG ("Unintentional error on sequence transformation "
                       "(char -> no).");
      return EXIT_FAILURE;
   }
   
   for (l = 0; l < rna_get_size (test_obj); l++)
   {
      if ( rna_get_sequence_base (l, test_obj) 
           != alphabet_base_2_no (test_string[l], sigma))
      {
         THROW_ERROR_MSG ("Unintentional error on sequence transformation, "
                          "base type not preserved for pos. %lu: src = %d "
                          "dest = %d", l,
                          alphabet_base_2_no (test_string[l], sigma),
                          rna_get_sequence_base (l, test_obj));
         return EXIT_FAILURE;         
      }
   }

   /* test back transformation */
   error = rna_transform_sequence_2_bases (sigma, test_obj);
   if (error)
   {
      THROW_ERROR_MSG ("Unintentional error on sequence transformation "
                       "(no -> char).");
      return EXIT_FAILURE;
   }
   if (strcmp (rna_get_sequence (test_obj), test_string) != 0)
   {
      THROW_ERROR_MSG ("Unintentional error: Sequence is changed after failed "
                       "transformation: src = \'%s\' dest = \'%s\'",
                       test_string, rna_get_sequence (test_obj));
      return EXIT_FAILURE;
   }

   rna_delete (test_obj);
   alphabet_delete (sigma);

   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}
