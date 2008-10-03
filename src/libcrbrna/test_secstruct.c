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
 *  @file libcrbrna/test_secstruct.c
 *
 *  @brief Test program for the secstruct module
 *
 *  Module: secstruct
 *
 *  Library: crbrna
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-09-16
 *
 *
 *  Revision History:
 *         - 2008Sep16 bienert: created
 *
 */


#include <config.h>
#include <stdlib.h>
#include <libcrbbasic/crbbasic.h>
#include "rna.h"
#include "secstruct.h"

int main(int argc __attribute__((unused)),char *argv[] __attribute__((unused)))
{
   Rna* rna;
   SecStruct* structure;
   /*char test_string[] = "(((...((((((.....))).....))).....)))";*/
   char test_string[] = 
     "..(((...(((...)))...(((...(((...)))...(((...((((((.....))).....))).....)))...(((...)))...)))...(((...)))...)))";

   unsigned long stack[][2] = { { 2,109}, { 3,108}, { 8, 16}, { 9,15}, {20,91},
                                {21, 90}, {26, 34}, {27, 33}, {38,73}, {39,72},
                                {44, 65}, {45, 64}, {47, 57}, {48,56}, {77,85},
                                {78, 84}, {95,103}, {96,102} };
   unsigned long noof_st = sizeof (stack) / sizeof (stack[0]);

   /* store: 5', 3', size  */
   unsigned long hairpin[][3] = { {10,14,3}, {28,32,3}, {49,55,5}, {79,83,3},
                                  {97,101,3} };
   unsigned long noof_hp = sizeof (hairpin) / sizeof (hairpin[0]);

   /* multiloop 2 */
   unsigned long ml2_stems[][2] = { {26,34}, {38,73}, {77,85}, {22,89} };
   unsigned long noof_ml2_st = sizeof (ml2_stems) / sizeof (ml2_stems[0]);

   /* dangles */
   unsigned long ml2_dangles[][4] = { {26, 34, 25, 35}, {38, 73, 37, 74},
                                      {77, 85, 76, 86}, {89, 22, 88, 23} };
   unsigned long noof_ml2_dl = sizeof (ml2_dangles) / sizeof (ml2_dangles[0]);

   /* multiloop 1 */
   unsigned long ml1_stems[][2] = { { 8,16}, {20,91}, {95,103}, {4,107} };
   unsigned long noof_ml1_st = sizeof (ml1_stems) / sizeof (ml1_stems[0]);

   /* dangles */
   unsigned long ml1_dangles[][4] = { { 8,  16, 7, 17}, {20, 91, 19, 92},
                                      {95, 103, 94, 104}, {107, 4, 106, 5} };
   unsigned long noof_ml1_dl = sizeof (ml1_dangles) / sizeof (ml1_dangles[0]);


   unsigned long i;
   int error;

   rna = RNA_NEW;
   if (rna == NULL)
   {
      return EXIT_FAILURE;
   }

   error = RNA_INIT_PAIRLIST_VIENNA(test_string, strlen (test_string), rna);

   if (error)
   {
      rna_delete (rna);
      return EXIT_FAILURE;
   }

   structure = SECSTRUCT_NEW;

   if (structure == NULL)
   {
      return EXIT_FAILURE;
   }

   /* try to fetch a structure */
   THROW_WARN_MSG ("Trying to store structure: \"%s\"", test_string);
   error = secstruct_find_interactions (rna_get_pairlist (rna),
                                        rna_get_size (rna),
                                        structure);

   if (error)
   {
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }

   /* check hairpins */
   if (secstruct_get_noof_hairpins (structure) != noof_hp)
   {
      THROW_ERROR_MSG ("Failed find right number of hairpins. Expected: 3, "
                       "found: %lu", secstruct_get_noof_hairpins (structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }

   for (i = 0; i < noof_hp; i++)
   {
      if (  (secstruct_get_i_start_hairpin (i, structure) != hairpin[i][0])
          ||(secstruct_get_i_end_hairpin (i, structure) != hairpin[i][1])
          ||(secstruct_get_i_size_hairpin (i, structure) != hairpin[i][2]))
      {
         THROW_ERROR_MSG ("Hairpin was stored with wrong geometry, supposed: "
                          "start %lu, end %lu, size %lu, is: %lu, %lu, %lu",
                          hairpin[i][0], hairpin[i][1], hairpin[i][2],
                          secstruct_get_i_start_hairpin (i, structure),
                          secstruct_get_i_end_hairpin (i, structure),
                          secstruct_get_i_size_hairpin (i, structure));
         secstruct_delete (structure);
         return EXIT_FAILURE;
      }
   }
   
   if (secstruct_get_noof_bulges (structure) != 1)
   {
      THROW_ERROR_MSG ("Failed find right number of bulge loops. Expected: 1, "
                       "found: %lu", secstruct_get_noof_bulges (structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }

   if (   (secstruct_get_i_start_bulge (0, structure) != 46)
       || (secstruct_get_i_end_bulge (0, structure) != 63)
       || (secstruct_get_i_size_bulge (0, structure) != 5))
   {
      THROW_ERROR_MSG ("Bulge loop was stored with wrong geometry, supposed: "
                       "start 46, end 63, size 5, is: %lu, %lu, %lu",
                       secstruct_get_i_start_bulge (0, structure),
                       secstruct_get_i_end_bulge (0, structure),
                       secstruct_get_i_size_bulge (0, structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;      
   }

   if (secstruct_get_noof_internals (structure) != 1)
   {
      THROW_ERROR_MSG ("Failed find right number of internal loops. Expected: "
                       "1, found: %lu",
                       secstruct_get_noof_internals (structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }

   if (   (secstruct_get_i_start_internal (0, structure) != 40)
       || (secstruct_get_i_end_internal (0, structure) != 71)
       || (secstruct_get_i_size_internal (0, structure) != 3))
   {
      THROW_ERROR_MSG ("Internal loop was stored with wrong geometry, "
                       "expected: start 40, end 71, size 3, is: %lu, %lu, %lu",
                       secstruct_get_i_start_internal (0, structure),
                       secstruct_get_i_end_internal (0, structure),
                       secstruct_get_i_size_internal (0, structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;      
   }

   /* check stacks */
   if (secstruct_get_noof_stacks (structure) != noof_st)
   {
      THROW_ERROR_MSG ("Failed to find right number of stacks. Expected: %lu, "
                       "found: %lu", noof_st,
                       secstruct_get_noof_stacks (structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;      
   }

   for (i = 0; i < noof_st; i++)
   {
      if (   (secstruct_get_i_5p_stack (i, structure) != stack[i][0]) 
          || (secstruct_get_i_3p_stack (i, structure) != stack[i][1]))
      {
         THROW_ERROR_MSG ("Stack %lu was stored with wrong geometry, expected: "
                          "5'- %lu, 3'- %lu, found: %lu, %lu", i,
                          stack[i][0], stack[i][1],
                          secstruct_get_i_5p_stack (i, structure),
                          secstruct_get_i_3p_stack (i, structure));
         secstruct_delete (structure);
         return EXIT_FAILURE;         
      }
   }

   /* check multiloops */
   if (secstruct_get_noof_multiloops (structure) != 2)
   {
      THROW_ERROR_MSG ("Failed to find right number of multiloops. Expected: "
                       "1, found: %lu",
                       secstruct_get_noof_multiloops (structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;      
   }

   /* check unpaired ml1 */
   if (secstruct_get_i_noof_unpaired_multiloop (0, structure) != 12)
   {
      THROW_ERROR_MSG ("Failed to find right number of unpaired bases for "
                       "multiloop 0. Expected: 12, found: %lu",
                       secstruct_get_i_noof_unpaired_multiloop (0, structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }

   /* check stems ml1 */
   if (secstruct_get_i_noof_stems_multiloop (0, structure) != noof_ml1_st)
   {
      THROW_ERROR_MSG ("Failed to find right number of base pairs for "
                       "multiloop 0. Expected: %lu, found: %lu",
                       noof_ml1_st,
                       secstruct_get_i_noof_stems_multiloop (0, structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }

   for (i = 0; i < noof_ml1_st; i++)
   {
      if (  (secstruct_get_i_5p_stem_multiloop (i, 0, structure)
             != ml1_stems[i][0])
          ||(secstruct_get_i_3p_stem_multiloop (i, 0, structure)
             != ml1_stems[i][1]))
      {
         THROW_ERROR_MSG ("Stem %lu of multiloop 0 was stored with wrong "
                          "geometry, expected: 5'- %lu, 3'- %lu, found: %lu, "
                          "%lu", i, ml1_stems[i][0], ml1_stems[i][1],
                          secstruct_get_i_5p_stem_multiloop (i, 0, structure),
                          secstruct_get_i_3p_stem_multiloop (i, 0, structure));
         secstruct_delete (structure);
         return EXIT_FAILURE;    
      }
   }

   /* check dangles ml1 */
   /* 5' */
   if (secstruct_get_i_noof_5pdangles_multiloop (1, structure) != noof_ml1_dl)
   {
      THROW_ERROR_MSG ("Failed to find right number of 5' dangles for "
                       "multiloop 0. Expected: %lu, found: %lu",
                       noof_ml1_dl,
                       secstruct_get_i_noof_5pdangles_multiloop (0, structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }   

   /* 3' */
   if (secstruct_get_i_noof_3pdangles_multiloop (0, structure) != noof_ml1_dl)
   {
      THROW_ERROR_MSG ("Failed to find right number of 3' dangles for "
                       "multiloop 0. Expected: %lu, found: %lu",
                       noof_ml1_dl,
                       secstruct_get_i_noof_3pdangles_multiloop (0, structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }

   for (i = 0; i < noof_ml1_dl; i++)
   {
      /* 5' */
      if (   (secstruct_get_i_5p_5pdangle_multiloop (i, 0, structure)
              != ml1_dangles[i][0])
          || (secstruct_get_i_3p_5pdangle_multiloop (i, 0, structure)
              != ml1_dangles[i][1]) 
          || (secstruct_get_i_dangle_5pdangle_multiloop (i, 0, structure)
              != ml1_dangles[i][2]))
      {
         THROW_ERROR_MSG ("5' dangle %lu of multiloop 0 was stored with wrong "
                          "geometry, expected: 5'- %lu, 3'- %lu, dangling base "
                          "- %lu, found: %lu, %lu, %lu", i,
                          ml1_dangles[i][0],
                          ml1_dangles[i][1],
                          ml1_dangles[i][2],
                          secstruct_get_i_5p_5pdangle_multiloop (i, 0,
                                                                 structure),
                          secstruct_get_i_3p_5pdangle_multiloop (i, 0,
                                                                 structure),
                          secstruct_get_i_dangle_5pdangle_multiloop (i, 0,
                                                                    structure));
         secstruct_delete (structure);
         return EXIT_FAILURE;
      }

      /* 3' */
      if (  (secstruct_get_i_5p_3pdangle_multiloop (i, 0, structure)
              != ml1_dangles[i][0])
          || (secstruct_get_i_3p_3pdangle_multiloop (i, 0, structure)
              != ml1_dangles[i][1])
          || (secstruct_get_i_dangle_3pdangle_multiloop (i, 0, structure)
              != ml1_dangles[i][3]))
      {
         THROW_ERROR_MSG ("3' dangle %lu of multiloop 0 was stored with wrong "
                          "geometry, expected: 5'- %lu, 3'- %lu, dangling base "
                          "- %lu, found: %lu, %lu, %lu", i,
                          ml1_dangles[i][0],
                          ml1_dangles[i][1],
                          ml1_dangles[i][3],
                          secstruct_get_i_5p_3pdangle_multiloop (i, 0,
                                                                 structure),
                          secstruct_get_i_3p_3pdangle_multiloop (i, 0,
                                                                 structure),
                          secstruct_get_i_dangle_3pdangle_multiloop (i, 0,
                                                                    structure));
         secstruct_delete (structure);
         return EXIT_FAILURE;
      }
   }

   /* check unpaired ml2 */
   if (secstruct_get_i_noof_unpaired_multiloop (1, structure) != 12)
   {
      THROW_ERROR_MSG ("Failed to find right number of unpaired bases for "
                       "multiloop 1. Expected: 12, found: %lu",
                       secstruct_get_i_noof_unpaired_multiloop (1, structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }

   /* check stems ml2 */
   if (secstruct_get_i_noof_stems_multiloop (0, structure) != noof_ml2_st)
   {
      THROW_ERROR_MSG ("Failed to find right number of base pairs for "
                       "multiloop 1. Expected: %lu, found: %lu",
                       noof_ml2_st,
                       secstruct_get_i_noof_stems_multiloop (0, structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }

   for (i = 0; i < noof_ml2_st; i++)
   {
      if (  (secstruct_get_i_5p_stem_multiloop (i, 1, structure)
             != ml2_stems[i][0])
          ||(secstruct_get_i_3p_stem_multiloop (i, 1, structure)
             != ml2_stems[i][1]))
      {
         THROW_ERROR_MSG ("Stem %lu of multiloop 1 was stored with wrong "
                          "geometry, expected: 5'- %lu, 3'- %lu, found: %lu, "
                          "%lu", i, ml2_stems[i][0], ml2_stems[i][1],
                          secstruct_get_i_5p_stem_multiloop (i, 1, structure),
                          secstruct_get_i_3p_stem_multiloop (i, 1, structure));
         secstruct_delete (structure);
         return EXIT_FAILURE;
      }
   }

   /* check dangles ml2 */
   /* 5' */
   if (secstruct_get_i_noof_5pdangles_multiloop (1, structure) != noof_ml2_dl)
   {
      THROW_ERROR_MSG ("Failed to find right number of 5' dangles for "
                       "multiloop 1. Expected: %lu, found: %lu",
                       noof_ml2_dl,
                       secstruct_get_i_noof_5pdangles_multiloop (0, structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }

   /* 3' */
   if (secstruct_get_i_noof_3pdangles_multiloop (1, structure) != noof_ml2_dl)
   {
      THROW_ERROR_MSG ("Failed to find right number of 3' dangles for "
                       "multiloop 1. Expected: %lu, found: %lu",
                       noof_ml2_dl,
                       secstruct_get_i_noof_3pdangles_multiloop (0, structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }

   for (i = 0; i < noof_ml2_dl; i++)
   {
      /* 5' */
      if (   (secstruct_get_i_5p_5pdangle_multiloop (i, 1, structure)
              != ml2_dangles[i][0])
          || (secstruct_get_i_3p_5pdangle_multiloop (i, 1, structure)
              != ml2_dangles[i][1]) 
          || (secstruct_get_i_dangle_5pdangle_multiloop (i, 1, structure)
              != ml2_dangles[i][2]))
      {
         THROW_ERROR_MSG ("5' dangle %lu of multiloop 1 was stored with wrong "
                          "geometry, expected: 5'- %lu, 3'- %lu, dangling base "
                          "- %lu, found: %lu, %lu, %lu", i,
                          ml2_dangles[i][0],
                          ml2_dangles[i][1],
                          ml2_dangles[i][2],
                          secstruct_get_i_5p_5pdangle_multiloop (i, 1,
                                                                 structure),
                          secstruct_get_i_3p_5pdangle_multiloop (i, 1,
                                                                 structure),
                          secstruct_get_i_dangle_5pdangle_multiloop (i, 1,
                                                                    structure));
         secstruct_delete (structure);
         return EXIT_FAILURE;
      }

      /* 3' */
      if (  (secstruct_get_i_5p_3pdangle_multiloop (i, 1, structure)
              != ml2_dangles[i][0])
          || (secstruct_get_i_3p_3pdangle_multiloop (i, 1, structure)
              != ml2_dangles[i][1])
          || (secstruct_get_i_dangle_3pdangle_multiloop (i, 1, structure)
              != ml2_dangles[i][3]))
      {
         THROW_ERROR_MSG ("3' dangle %lu of multiloop 1 was stored with wrong "
                          "geometry, expected: 5'- %lu, 3'- %lu, dangling base "
                          "- %lu, found: %lu, %lu, %lu", i,
                          ml2_dangles[i][0],
                          ml2_dangles[i][1],
                          ml2_dangles[i][3],
                          secstruct_get_i_5p_3pdangle_multiloop (i, 1,
                                                                 structure),
                          secstruct_get_i_3p_3pdangle_multiloop (i, 1,
                                                                 structure),
                          secstruct_get_i_dangle_3pdangle_multiloop (i, 1,
                                                                    structure));
         secstruct_delete (structure);
         return EXIT_FAILURE;
      }
   }

   /* check for external loop */
   if (secstruct_get_i_noof_unpaired_extloop (structure) != 2)
   {
      THROW_ERROR_MSG ("Wrong no. of unpaired bases found for the external "
                       "loop. Expected: 2, found: %lu",
                       secstruct_get_i_noof_unpaired_extloop (structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }

   if (secstruct_get_i_noof_stems_extloop (structure) != 1)
   {
      THROW_ERROR_MSG ("Failed to find right number of base pairs for the "
                       "external loop. Expected: 1, found: %lu",
                       secstruct_get_i_noof_stems_extloop (structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }

   if (   (secstruct_get_i_5p_stem_extloop (0, structure) != 2)
       || (secstruct_get_i_3p_stem_extloop (0, structure) != 109))
   {
         THROW_ERROR_MSG ("Stem 0 of external loop was stored with wrong "
                          "geometry, expected: 5'- 2, 3'- 109, found: %lu, "
                          "%lu",
                          secstruct_get_i_5p_stem_extloop (0, structure),
                          secstruct_get_i_3p_stem_extloop (0, structure));
         secstruct_delete (structure);
         return EXIT_FAILURE;
   }

   /* check dangles */
   if (secstruct_get_noof_5pdangles_extloop (structure) != 1)
   {
      THROW_ERROR_MSG ("Wrong no. of 5' dangling ends found for external "
                       "loop. Expected: 1, found: %lu",
                       secstruct_get_noof_5pdangles_extloop (structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }

   if (secstruct_get_noof_3pdangles_extloop (structure) != 0)
   {
      THROW_ERROR_MSG ("Wrong no. of 3' dangling ends found for external "
                       "loop. Expected: 0, found: %lu",
                       secstruct_get_noof_3pdangles_extloop (structure));
      secstruct_delete (structure);
      return EXIT_FAILURE;
   }  

   if (   (secstruct_get_i_5p_5pdangle_extloop (0, structure) != 2) 
       || (secstruct_get_i_3p_5pdangle_extloop (0, structure) != 109))
   {
         THROW_ERROR_MSG ("5' dangle 0 of external loop was stored with wrong "
                          "geometry, expected: 5'- 2, 3'- 109, found: %lu, "
                          "%lu",
                          secstruct_get_i_5p_5pdangle_extloop (0, structure),
                          secstruct_get_i_3p_5pdangle_extloop (0, structure));
         secstruct_delete (structure);
         return EXIT_FAILURE;      
   }

   secstruct_delete (structure);
   rna_delete (rna);

   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}
