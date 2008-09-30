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

   sigma = ALPHABET_NEW_PAIR ("AUGC", "augc", 4);
   if (sigma == NULL)
   {
      THROW_ERROR_MSG ("Could not create alphabet");
      return EXIT_FAILURE;
   }
   
   scores = NN_SCORES_NEW_INIT (sigma);
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

   alphabet_delete (sigma);
   nn_scores_delete (scores);

   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}
