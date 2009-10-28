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
 *  @file libcrbrna/test_nn_scores.c
 *
 *  @brief Test program for the nn_scores module
 *
 *  Module: nn_scores
 *
 *  Library: crbrna
 *
 *  Project: CoRB - COllection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-07-03
 *
 *
 *  Revision History:
 *         - 2008Jul03 bienert: created
 *         - 2009Oct08 bienert: Added test for hashed tetra loops
 *
 */


#include <config.h>
#include <stdlib.h>
#include <libcrbbasic/crbbasic.h>
#include "alphabet.h"
#include "nn_scores.h"

int main(int argc __attribute__((unused)),char *argv[] __attribute__((unused)))
{

   Alphabet* sigma;
   NN_scores* scores;
   unsigned long bp, i;
   unsigned long no_of_tloops;
   char tloop[7] = {'\0', '\0', '\0', '\0', '\0', '\0', '\0'};
   bool is_tloop;
   const char* ref_t_loop;
   char asize;

   sigma = ALPHABET_NEW_PAIR ("ACGU", "acgu", 4);
   if (sigma == NULL)
   {
      THROW_ERROR_MSG ("Could not create alphabet");
      return EXIT_FAILURE;
   }
   
   scores = NN_SCORES_NEW_INIT (0, sigma);
   if (scores == NULL)
   {
      THROW_ERROR_MSG ("Could not setup scoring scheme");
      alphabet_delete (sigma);
      return EXIT_FAILURE;
   }

   mprintf ("Allowed base pairs:\n");
   nn_scores_fprintf_bp_allowed (stdout, scores, sigma);
   mprintf ("Base pair indices:\n");
   nn_scores_fprintf_bp_idx (stdout, scores, sigma);
   mprintf ("Stacking energies:\n");
   nn_scores_fprintf_G_stack (stdout, scores, sigma);
   mprintf ("Mismatch stacking energies:\n");
   nn_scores_fprintf_mm_G_stack (stdout, scores, sigma);
   mprintf ("Hairpin loop energies:\n");
   mprintf ("Size: Score\n");
   nn_scores_fprintf_G_hairpin_loop (stdout, scores);
   mprintf ("Mismatch hairpin energies:\n");
   nn_scores_fprintf_G_mismatch_hairpin (stdout, scores, sigma);
   mprintf ("Penalties for non-GC closing base pairs:\n");
   nn_scores_fprintf_non_gc_penalty_for_bp(stdout, scores, sigma);
   mprintf ("Bulge loop energies:\n");
   mprintf ("Size: Score\n");
   nn_scores_fprintf_G_bulge_loop (stdout, scores);
   mprintf ("1x1 internal loop energies:\n");
   nn_scores_fprintf_G_int11(stdout, scores, sigma);
   mprintf ("2x1 internal loop energies:\n");
   nn_scores_fprintf_G_int21(stdout, scores, sigma);
   mprintf ("2x2 internal loop energies:\n");
   nn_scores_fprintf_G_int22(stdout, scores, sigma);
   mprintf ("Generic internal loop energies:\n");
   mprintf ("Size: Score\n");
   nn_scores_fprintf_G_internal_loop (stdout, scores);
   mprintf ("Mismatch interior energies:\n");
   nn_scores_fprintf_G_mismatch_interior (stdout, scores, sigma);
   mprintf ("5' dangling end energies:\n");
   nn_scores_fprintf_G_dangle5(stdout, scores, sigma);
   mprintf ("3' dangling end energies:\n");
   nn_scores_fprintf_G_dangle3(stdout, scores, sigma);
   mprintf ("Non-unitable nucleotides (nun) penalties:\n");
   nn_scores_fprintf_nun_penalties (stdout, scores, sigma);
   mprintf ("Tetra loop bonus energies:\n");
   nn_scores_fprintf_tetra_loop(stdout, scores, sigma);
   mprintf ("Tetra loop hash function:\n");
   nn_scores_fprintf_tetra_loop_hashfunction(stdout, scores, sigma);

   /* SB 09-10-08 START 
    * Testing the tetra loop hashfunction. If the function stops working,
    * modify nn_scores.c according to the instruction written therein and
    * comment out this test.
    */
   THROW_WARN_MSG ("Checking tetra loop scores.");
   asize = (signed) alphabet_size (sigma);
   no_of_tloops = nn_scores_get_no_of_tetra_loops (scores);
   /* loop over all possible closing base pairs */
   for (bp = 0; bp < nn_scores_no_allowed_basepairs (scores); bp++)
   {
      nn_scores_get_allowed_basepair (bp, &(tloop[0]), &(tloop[5]), scores);

      /* loop over all possible unpaired bases */
      for (tloop[1] = 0; tloop[1] < asize; tloop[1]++)
      {
         for (tloop[2] = 0; tloop[2] < asize; tloop[2]++)
         {
            for (tloop[3] = 0; tloop[3] < asize; tloop[3]++)
            {
               for (tloop[4] = 0; tloop[4] < asize; tloop[4]++)
               {
                  /* check each loop to be one of the 30 */
                  is_tloop = false;
                  for (i = 0; i < no_of_tloops; i++)
                  {
                     ref_t_loop = nn_scores_get_tetra_loop (i, scores);
                     if (  (tloop[0] == ref_t_loop[0])
                         &&(tloop[1] == ref_t_loop[1])
                         &&(tloop[2] == ref_t_loop[2])
                         &&(tloop[3] == ref_t_loop[3])
                         &&(tloop[4] == ref_t_loop[4])
                         &&(tloop[5] == ref_t_loop[5]))
                     {
                                          mprintf ("%c-%c%c%c%c-%c: %.2f\n",
                           alphabet_no_2_base (tloop[0], sigma),
                           alphabet_no_2_base (tloop[1], sigma),
                           alphabet_no_2_base (tloop[2], sigma),
                           alphabet_no_2_base (tloop[3], sigma),
                           alphabet_no_2_base (tloop[4], sigma),
                           alphabet_no_2_base (tloop[5], sigma),
                           nn_scores_get_G_tetra_loop (tloop, 0, scores));
                                          is_tloop = true;
                                          i = no_of_tloops + 1;
                     }
                  }

                  if (! is_tloop)
                  {
                     if ((nn_scores_get_G_tetra_loop (tloop, 0, scores) < 0.0f)
                      || (nn_scores_get_G_tetra_loop (tloop, 0, scores) > 0.0f))
                     {
                        THROW_ERROR_MSG ("Tetra loop \"%c-%c%c%c%c-%c\" is "
                                         "supposed to have a bonus score of "
                                         "0, but gets %.2f",
                                         alphabet_no_2_base (tloop[0], sigma),
                                         alphabet_no_2_base (tloop[1], sigma),
                                         alphabet_no_2_base (tloop[2], sigma),
                                         alphabet_no_2_base (tloop[3], sigma),
                                         alphabet_no_2_base (tloop[4], sigma),
                                         alphabet_no_2_base (tloop[5], sigma),
                                   nn_scores_get_G_tetra_loop (tloop, 0, scores)
                           );
                        alphabet_delete (sigma);
                        nn_scores_delete (scores);
                        
                        FREE_MEMORY_MANAGER;
                        
                        return EXIT_FAILURE;
                     }
                  }
               }
            }
         }
      }
   }
   /* SB 09-10-08 END */

   alphabet_delete (sigma);
   nn_scores_delete (scores);

   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}
