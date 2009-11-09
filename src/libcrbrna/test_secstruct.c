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
   unsigned long i1, j1, i2, j2, size1, size2;
   Rna* rna;
   SecStruct* structure;
   SecStructFtrs type;
   unsigned long idx;

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

secstruct_get_geometry_internal (&i1, &j1, &i2, &j2, &size1, &size2, 0,
                                 structure);

   if (   (i1 != 40)
       || (j1 != 71)
       || (size1 != 3))
   {
      THROW_ERROR_MSG ("Internal loop was stored with wrong geometry, "
                       "expected: start 40, end 71, size 3, is: %lu, %lu, %lu",
                       i1, j1, size1);
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

   if (secstruct_get_noof_stems_extloop (structure) != 1)
   {
      THROW_ERROR_MSG ("Failed to find right number of base pairs for the "
                       "external loop. Expected: 1, found: %lu",
                       secstruct_get_noof_stems_extloop (structure));
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

   mprintf ("Stacked base piars:\n");
   secstruct_fprintf_stacks (stdout, structure);
   mprintf ("\nHairpin loops:\n");
   secstruct_fprintf_hairpins (stdout, structure);
   mprintf ("\nBulge loops:\n");
   secstruct_fprintf_bulges (stdout, structure);
   mprintf ("\nInternal loops:\n");
   secstruct_fprintf_internals (stdout, structure);
   mprintf ("\nExternal loop:\n");
   secstruct_fprintf_external (stdout, structure);
   mprintf ("\nMulti loops:\n");
   secstruct_fprintf_multiloops (stdout, structure);
   mprintf ("\n\nSequence position to structural feature map:\n");
   /*secstruct_delete_element (13,
                             SCSTRCT_STACK,
                             structure);*/
   /*secstruct_delete_stack (2, structure);
     secstruct_delete_stack (4, structure);*/
   /*secstruct_delete_stack (12, structure);
     secstruct_delete_stack (13, structure);*/
   secstruct_fprintf_seqpos_map (stdout, structure);

   /* test seq pos -> feature technique */
   type = secstruct_get_feature_at_pos (0, &idx, structure);
   if (type != SCSTRCT_EXTERNAL)
   {
      THROW_ERROR_MSG("Sequence position 0 is not mapped to an external loop.");
      return EXIT_FAILURE;
   }
   type = secstruct_get_feature_at_pos (3, &idx, structure);
   if (type != SCSTRCT_STACK)
   {
      THROW_ERROR_MSG("Sequence position 3 is not mapped to a stacked base.");
      return EXIT_FAILURE;
   }
   type = secstruct_get_feature_at_pos (7, &idx, structure);
   if (type != SCSTRCT_MULTI)
   {
      THROW_ERROR_MSG("Sequence position 7 is not mapped to a multi loop.");
      return EXIT_FAILURE;
   }
   type = secstruct_get_feature_at_pos (70, &idx, structure);
   if (type != SCSTRCT_INTERNAL)
   {
      THROW_ERROR_MSG ("Sequence position 70 is not mapped to an internal "
                       "loop.");
      return EXIT_FAILURE;
   }
   type = secstruct_get_feature_at_pos (82, &idx, structure);
   if (type != SCSTRCT_HAIRPIN)
   {
      THROW_ERROR_MSG ("Sequence position 82 is not mapped to a hairpin loop.");
      return EXIT_FAILURE;
   }
   type = secstruct_get_feature_at_pos (62, &idx, structure);
   if (type != SCSTRCT_BULGE)
   {
      THROW_ERROR_MSG ("Sequence position 62 is not mapped to a bulge loop.");
      return EXIT_FAILURE;
   }
   type = secstruct_get_feature_at_pos (95, &idx, structure);
   if (type != SCSTRCT_MTO)
   {
      THROW_ERROR_MSG ("Sequence position 95 is not mapped as multiple site.");
      return EXIT_FAILURE;
   }
   else
   {
      type = secstruct_get_feature_multi_1st (95, &idx, structure);
      if (type != SCSTRCT_MULTI && type != SCSTRCT_STACK)
      {
         THROW_ERROR_MSG ("Sequence position 95 is not mapped as multi loop "
                          "or stacked base pair.");
         return EXIT_FAILURE;
      }
      type = secstruct_get_feature_multi_2nd (95, &idx, structure);
      if (type != SCSTRCT_MULTI && type != SCSTRCT_STACK)
      {
         THROW_ERROR_MSG ("Sequence position 95 is not mapped as multi loop "
                          "or stacked base pair.");
         return EXIT_FAILURE;
      }
   }

   /* now we try to delete some elements */
   /*error = secstruct_delete_element (0, 0, structure);*/

   secstruct_delete (structure);
   rna_delete (rna);

   FREE_MEMORY_MANAGER;

   return EXIT_SUCCESS;
}
/*
  0:   0   External loop
  1:   0   External loop
  2:   0   -> 0  External loop; 0  Stacked base pair;
  3:   1   Stacked base pair
  4:   0   Multi loop
  5:   0   Multi loop
  6:   0   Multi loop
  7:   0   Multi loop
  8:   2   -> 0  Multi loop; 2  Stacked base pair;
  9:   3   Stacked base pair
 10:   0   Hairpin loop
 11:   0   Hairpin loop
 12:   0   Hairpin loop
 13:   0   Hairpin loop
 14:   0   Hairpin loop
 15:   3   Stacked base pair
 16:   2   -> 0  Multi loop; 2  Stacked base pair;
 17:   0   Multi loop
 18:   0   Multi loop
 19:   0   Multi loop
 20:   4   -> 0  Multi loop; 4  Stacked base pair;
 21:   5   Stacked base pair
 22:   1   Multi loop
 23:   1   Multi loop
 24:   1   Multi loop
 25:   1   Multi loop
 26:   6   -> 1  Multi loop; 6  Stacked base pair;
 27:   7   Stacked base pair
 28:   1   Hairpin loop
 29:   1   Hairpin loop
 30:   1   Hairpin loop
 31:   1   Hairpin loop
 32:   1   Hairpin loop
 33:   7   Stacked base pair
 34:   6   -> 1  Multi loop; 6  Stacked base pair;
 35:   1   Multi loop
 36:   1   Multi loop
 37:   1   Multi loop
 38:   8   -> 1  Multi loop; 8  Stacked base pair;
 39:   9   Stacked base pair
 40:   0   Internal loop
 41:   0   Internal loop
 42:   0   Internal loop
 43:   0   Internal loop
 44:  10   -> 0  Internal loop; 10  Stacked base pair;
 45:  11   Stacked base pair
 46:   0   Bulge loop
 47:  12   -> 0  Bulge loop; 12  Stacked base pair;
 48:  13   Stacked base pair
 49:   2   Hairpin loop
 50:   2   Hairpin loop
 51:   2   Hairpin loop
 52:   2   Hairpin loop
 53:   2   Hairpin loop
 54:   2   Hairpin loop
 55:   2   Hairpin loop
 56:  13   Stacked base pair
 57:  12   -> 0  Bulge loop; 12  Stacked base pair;
 58:   0   Bulge loop
 59:   0   Bulge loop
 60:   0   Bulge loop
 61:   0   Bulge loop
 62:   0   Bulge loop
 63:   0   Bulge loop
 64:  11   Stacked base pair
 65:  10   -> 0  Internal loop; 10  Stacked base pair;
 66:   0   Internal loop
 67:   0   Internal loop
 68:   0   Internal loop
 69:   0   Internal loop
 70:   0   Internal loop
 71:   0   Internal loop
 72:   9   Stacked base pair
 73:   8   -> 1  Multi loop; 8  Stacked base pair;
 74:   1   Multi loop
 75:   1   Multi loop
 76:   1   Multi loop
 77:  14   -> 1  Multi loop; 14  Stacked base pair;
 78:  15   Stacked base pair
 79:   3   Hairpin loop
 80:   3   Hairpin loop
 81:   3   Hairpin loop
 82:   3   Hairpin loop
 83:   3   Hairpin loop
 84:  15   Stacked base pair
 85:  14   -> 1  Multi loop; 14  Stacked base pair;
 86:   1   Multi loop
 87:   0   Not assigned
 88:   1   Multi loop
 89:   1   Multi loop
 90:   5   Stacked base pair
 91:   4   -> 0  Multi loop; 4  Stacked base pair;
 92:   0   Multi loop
 93:   0   Multi loop
 94:   0   Multi loop
 95:  16   -> 0  Multi loop; 16  Stacked base pair;
 96:  17   Stacked base pair
 97:   4   Hairpin loop
 98:   4   Hairpin loop
 99:   4   Hairpin loop
100:   4   Hairpin loop
101:   4   Hairpin loop
102:  17   Stacked base pair
103:  16   -> 0  Multi loop; 16  Stacked base pair;
104:   0   Multi loop
105:   0   Not assigned
106:   0   Multi loop
107:   0   Multi loop
108:   1   Stacked base pair
109:   0   -> 0  External loop; 0  Stacked base pair;
 */
