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
   char* structure = "STRUCTURE";
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
   error = RNA_INIT_PAIRLIST_VIENNA(structure, strlen (structure), test_obj);
   if (error != ERR_RNA_VIENNA_FORMAT)
   {
      THROW_ERROR_MSG ("Producing error by using \'%s\' as a structure "
                       "string failed!", structure);
      return EXIT_FAILURE;
   }

   /* test with unbalanced structure 1 */
   THROW_WARN_MSG ("Testing recognition of unbalanced structure. You should "
                   "see an error message in the following.");
   structure = "........)";
   error = RNA_INIT_PAIRLIST_VIENNA(structure, strlen (structure), test_obj);
   if (error != ERR_RNA_VIENNA_MMC)
   {
      THROW_ERROR_MSG ("Producing error by using \'%s\' as a structure "
                       "string failed!", structure);
      return EXIT_FAILURE;
   }

   /* test with unbalanced structure 2 */
   THROW_WARN_MSG ("Testing recognition of unbalanced structure. You should "
                   "see an error message in the following.");
   structure = "(........";
   error = RNA_INIT_PAIRLIST_VIENNA(structure, strlen (structure), test_obj);
   if (error != ERR_RNA_VIENNA_MMO)
   {
      THROW_ERROR_MSG ("Producing error by using \'%s\' as a structure "
                       "string failed!", structure);
      return EXIT_FAILURE;
   }

   /* set correct structure */
   structure = "((.....))";
   error = RNA_INIT_PAIRLIST_VIENNA(structure, strlen (structure), test_obj);
   if (error)
   {
      THROW_ERROR_MSG ("Unintentional error while testing with structure "
                       "\'%s\'!", structure);
      return EXIT_FAILURE;      
   }

   /* testing accessing pairs */
   l = rna_base_pairs_with (0, test_obj);
   if (l != strlen (structure))
   {
      THROW_ERROR_MSG ("Wrong index for pairing partner found. Checked pos. "
                       "%d, expected partner index %lu, got %lu. Structure "
                       "was: \'%s\'!",
                       0, (unsigned long) strlen (structure), l, structure);
      return EXIT_FAILURE;       
   }
   r = rna_base_pairs_with (strlen (structure) - 1, test_obj);
   if (r != 1)
   {
      THROW_ERROR_MSG ("Wrong index for pairing partner found. Checked pos. "
                       "%d, expected partner index %d, got %lu. Structure "
                       "was: \'%s\'!",
                       (strlen (structure) - 1), 1, r, structure);
      return EXIT_FAILURE;       
   }

   rna_delete (test_obj);

   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}
