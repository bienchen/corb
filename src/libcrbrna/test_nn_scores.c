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
   mprintf ("Base pair indeces:\n");
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
   mprintf ("Bulge loop energies:\n");
   mprintf ("Size: Score\n");
   nn_scores_fprintf_G_bulge_loop (stdout, scores);
   mprintf ("Penalties for non-GC closing base pairs:\n");
   nn_scores_fprintf_non_gc_penalty_for_bp(stdout, scores, sigma);
   mprintf ("Tetra loop bonus energies:\n");
   nn_scores_fprintf_tetra_loop(stdout, scores, sigma);
   mprintf ("5' dangling end energies:\n");
   nn_scores_fprintf_G_dangle5(stdout, scores, sigma);
   mprintf ("3' dangling end energies:\n");
   nn_scores_fprintf_G_dangle3(stdout, scores, sigma);
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

   alphabet_delete (sigma);
   nn_scores_delete (scores);

   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}
