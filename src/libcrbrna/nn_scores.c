/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbrna/nn_scores.c
 *
 *  @brief Neares neighbour Model for evaluating RNA secondary structures
 *
 *  Module: nn_scores
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
#include <math.h>
#include <libcrbbasic/crbbasic.h>
#include "alphabet.h"
#include "nn_scores.h"


struct NN_scores {
      long** G_stack;
      unsigned long G_stack_size;
};


/**********************   Constructors and destructors   **********************/

/** @brief Create a new Nearest Neighbour scoring scheme.
 *
 * The constructor for @c NN_scores objects. If compiled with enabled
 * memory checking, @c file and @c line should point to the position where the
 * function was called. Both parameters are automatically set by using the
 * macro @c NN_SCORES_NEW.\n
 * Returns @c NULL on error.
 *
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
NN_scores*
nn_scores_new (const char* file, const int line)
{
   /* allocate 1 object */
   NN_scores* this = XOBJ_MALLOC(sizeof (NN_scores), file, line);
   
   if (this != NULL)
   {
      this->G_stack      = NULL;
      this->G_stack_size = 0;
   }

   return this;
}

/** @brief Create a new Nearest Neighbour scoring scheme with standard values.
 *
 * The constructor for an initialised @c NN_scores objects. If compiled with
 * enabled memory checking, @c file and @c line should point to the position
 * where the function was called. Both parameters are automatically set by
 * using the macro @c NN_SCORES_NEW_INIT.\n
 * As parameters for the canonical Watson-Crick base pairs, plus the G-U wobble
 * base pair, the stacking energies (table "stack_energies") from the Vienna
 * RNA package are used. For stacks consisting of only one base pair and a
 * mismatch, we use ... mismatch_interior \n
 * Returns @c NULL on error.
 *
 * @param[in] sigma alphabet.
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
NN_scores*
nn_socres_new_init (Alphabet* sigma, const char* file, const int line)
{
   NN_scores* this;
   unsigned long i;
   char a, u, g, c;

   assert (sigma != NULL);

   this = nn_scores_new (file, line);

   if (this != NULL)
   {
      if (! alphabet_is_standard_rna (sigma))
      {
         nn_scores_delete (this);
         return NULL;
      }

      /* prepare table for stacking energies */
      this->G_stack_size = alphabet_size (sigma);
      for (i = 1; i < (alphabet_size (sigma) / 2); i++)
      {
         this->G_stack_size *= alphabet_size (sigma);
      }
      this->G_stack = (long**) XMALLOC_2D (this->G_stack_size,
                                           this->G_stack_size,
                                           sizeof (long));
      this->G_stack_size *= this->G_stack_size;

      for (i = 0; i < this->G_stack_size; i++)
      {
         this->G_stack[0][i] = 0;
      }

      /* set stacking energies (DG) */
      a = alphabet_base_2_no('A', sigma);
      u = alphabet_base_2_no('U', sigma);
      g = alphabet_base_2_no('G', sigma);
      c = alphabet_base_2_no('C', sigma);

      /* 5' - ip - 3' */
      /* 3' - jq - 5' --> i < p < q < j */
      /*    = ij qp    */

      /* regular pairs */
      /* AU AU */
      /* 5'- AU 
             UA -5' */
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (a, u, sigma)] = -110;

      /* AU UA */
      /* 5'- AA
             UU -5' */
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (u, a, sigma)] = -90;

      /* AU UG */
      /* 5'- AG
             UU -5' */
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (u, g, sigma)] = -60;
     

      /* AU GU */
      /* 5'- AU
             UG -5'*/
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (g, u, sigma)] = -140;      

      /* AU CG */
      /* 5'- AG
             UC -5'*/
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (c, g, sigma)] = -210;        

      /* AU GC */
      /* 5'- AC
             UG -5'*/
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (g, c, sigma)] = -220; 

      /* UA AU */
      /* 5'- UU 
             AA  -5'*/
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (a, u, sigma)] = -90;

      /* UA UA */
      /* 5'- UA
             AU -5'*/
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (u, a, sigma)] = -130;

      /* UA UG */
      /* 5'- UG
             AU -5'*/
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (u, g, sigma)] = -100;

      /* UA GU */
      /* 5'- UU
             AG -5'*/
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (g, u, sigma)] = -130;

      /* UA CG */
      /* 5'- UG
             AC -5'*/
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (c, g, sigma)] = -210;

      /* UA GC */
      /* 5'- UC
             AG -5'*/
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (g, c, sigma)] = -240;
      
      /* UG AU */
      /* 5'- UU
             GA -5'*/
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (a, u, sigma)] = -60;

      /* UG UA */
      /* 5'- UA 
             GU -5'*/
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (u, a, sigma)] = -100;

      /* UG UG */
      /* 5'- UG 
             GU -5'*/
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (u, g, sigma)] = 30;

      /* UG GU*/
      /* 5'- UU 
             GG -5'*/
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (g, u, sigma)] = -50;

      /* UG CG */
      /* 5'- UG 
             GC -5'*/
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (c, g, sigma)] = -140;

      /* UG GC */
      /* 5'- UC 
             GG -5'*/
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (g, c, sigma)] = -150;

      /* GU AU */
      /* 5'- GU 
             UA -5'*/
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (a, u, sigma)] = -140;

      /* GU UA */
      /* 5'- GA 
             UU -5'*/
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (u, a, sigma)] = -130;

      /* GU UG */
      /* 5'- GG 
             UU -5'*/
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (u, g, sigma)] = -50;

      /* GU GU */
      /* 5'- GU 
             UG -5'*/
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (g, u, sigma)] = 130;

      /* GU CG */
      /* 5'- GG 
             UC -5'*/
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (c, g, sigma)] = -210;

      /* GU GC */
      /* 5'- GC 
             UG -5'*/
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (g, c, sigma)] = -250;

      /* CG AU */
      /* 5'- C
             G -5'*/
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (a, u, sigma)] = -210;

      /* CG UA */
      /* 5'- CA
             GU -5'*/
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (u, a, sigma)] = -210;

      /* CG UG */
      /* 5'- CG
             GU -5'*/
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (u, g, sigma)] = -140;

      /* CG GU */
      /* 5'- CU
             GG -5'*/
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (g, u, sigma)] = -210;

      /* CG CG */
      /* 5'- CG
             GC -5'*/
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (c, g, sigma)] = -240;

      /* CG GC */
      /* 5'- CC
             GG -5'*/
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (g, c, sigma)] = -330;

      /* GC AU */
      /* 5'- GU
             CA -5'*/
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (a, u, sigma)] = -220;

      /* GC UA */
      /* 5'- GA
             CU -5'*/
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (u, a, sigma)] = -240;

      /* GC UG */
      /* 5'- GG
             CU -5'*/
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (u, g, sigma)] = -150;

      /* GC GU */
      /* 5'- GU
             CG -5'*/
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (g, u, sigma)] = -250;

      /* GC CG */
      /* 5'- GG
             CC -5'*/
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (c, g, sigma)] = -330;

      /* GC GC */
      /* 5'- GC
             CG -5'*/
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (g, c, sigma)] = -340;

      /* stacks containing a mismatch */
      /* mi: param from mismatch_interior table */
      /* mh: param from mismatch_hairpin */
      /* for testing at start: each pair to be set assert if not 0!!! */

      /* CG AA */
      /* mi: 0 mh: -150 */
      assert (this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                           [nn_scores_bp_2_idx (a, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (a, a, sigma)] = -75;

      /* CG AC */
      /* mi: 0 mh: -150 */
      assert (this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                           [nn_scores_bp_2_idx (a, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (a, c, sigma)] = -75;

      /* CG AG */
      /* mi: -110 mh: -140 */
      assert (this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                           [nn_scores_bp_2_idx (a, g, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (a, g, sigma)] = -125;

      /* CG CA */
      /* mi: 0 mh: -100 */
      assert (this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                           [nn_scores_bp_2_idx (c, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (c, a, sigma)] = -50;

      /* CG CC */
      /* mi: 0 mh: -90 */
      assert (this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                           [nn_scores_bp_2_idx (c, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (c, c, sigma)] = -45;

      /* CG CU */
      /* mi: 0 mh: -80 */
      assert (this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                           [nn_scores_bp_2_idx (c, u, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (c, u, sigma)] = -40;

      /* CG GA */
      /* mi: -110 mh: -220 */
      assert (this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                           [nn_scores_bp_2_idx (g, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (g, a, sigma)] = -165;

      /* CG GG */
      /* mi: 0 mh: -160 */
      assert (this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                           [nn_scores_bp_2_idx (g, g, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (g, g, sigma)] = -80;

      /* CG UC */
      /* mi: 0 mh: -140 */
      assert (this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                           [nn_scores_bp_2_idx (u, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (u, c, sigma)] = -70;

      /* CG UU */
      /* mi: -70 mh: -200 */
      assert (this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                           [nn_scores_bp_2_idx (u, u, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (c, g, sigma)]
                   [nn_scores_bp_2_idx (u, u, sigma)] = -135;

      /* GC AA */
      /* mi: 0 mh: -110 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                           [nn_scores_bp_2_idx (a, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (a, a, sigma)] = -55;

      /* GC AC */
      /* mi: 0 mh: -150 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                           [nn_scores_bp_2_idx (a, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (a, c, sigma)] = -75;

      /* GC AG */
      /* mi: -110 mh: -130 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                           [nn_scores_bp_2_idx (a, g, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (a, g, sigma)] = -120;

      /* GC CA */
      /* mi: 0 mh: -110 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                           [nn_scores_bp_2_idx (c, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (c, a, sigma)] = -55;
     
      /* GC CC */
      /* mi: 0 mh: -70 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                           [nn_scores_bp_2_idx (c, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (c, c, sigma)] = -35;

      /* GC CU */
      /* mi: 0 mh: -50 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                           [nn_scores_bp_2_idx (c, u, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (c, u, sigma)] = -25;

      /* GC GA */
      /* mi: -110 mh: -240 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                           [nn_scores_bp_2_idx (g, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (g, a, sigma)] = -175;

      /* GC GG */
      /* mi: 0 mh: -140 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                           [nn_scores_bp_2_idx (g, g, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (g, g, sigma)] = -70;

      /* GC UC */
      /* mi: 0 mh: -100 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                           [nn_scores_bp_2_idx (u, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (u, c, sigma)] = -50;

      /* GC UU */
      /* mi: -70 mh: -150 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                           [nn_scores_bp_2_idx (u, u, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, c, sigma)]
                   [nn_scores_bp_2_idx (u, u, sigma)] = -110;

      /* GU AA */
      /* mi: 70 mh: 20 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                           [nn_scores_bp_2_idx (a, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (a, a, sigma)] = 45;

      /* GU AC */
      /* mi: 70 mh: -50 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                           [nn_scores_bp_2_idx (a, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (a, c, sigma)] = 10;

      /* GU AG */
      /* mi: -40 mh: -30 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                           [nn_scores_bp_2_idx (a, g, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (a, g, sigma)] = -35;

      /* GU CA */
      /* mi: 70 mh: -10 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                           [nn_scores_bp_2_idx (c, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (c, a, sigma)] = 30;

      /* GU CC */
      /* mi: 70 mh: -20 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                           [nn_scores_bp_2_idx (c, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (c, c, sigma)] = 25;

      /* GU CU */
      /* mi: 70 mh: -20 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                           [nn_scores_bp_2_idx (c, u, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (c, u, sigma)] = 25;

      /* GU GA */
      /* mi: -40 mh: -90 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                           [nn_scores_bp_2_idx (g, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (g, a, sigma)] = -65;

      /* GU GG */
      /* mi: 70 mh: -30 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                           [nn_scores_bp_2_idx (g, g, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (g, g, sigma)] = 20;

      /* GU UC */
      /* mi: 70 mh: -30 */     
      assert (this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                           [nn_scores_bp_2_idx (u, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (u, c, sigma)] = 20;

      /* GU UU */
      /* mi: 0 mh: -110 */
      assert (this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                           [nn_scores_bp_2_idx (u, u, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (g, u, sigma)]
                   [nn_scores_bp_2_idx (u, u, sigma)] = -55;

      /* UG AA */
      /* mi: 70 mh: -50 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                           [nn_scores_bp_2_idx (a, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (a, a, sigma)] = 10;

      /* UG AC */
      /* mi: 70 mh: -30 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                           [nn_scores_bp_2_idx (a, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (a, c, sigma)] = 20;

      /* UG AG */
      /* mi: -40 mh: -60 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                           [nn_scores_bp_2_idx (a, g, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (a, g, sigma)] = -50;

      /* UG CA */
      /* mi: 70 mh: -20 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                           [nn_scores_bp_2_idx (c, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (c, a, sigma)] = 25;      

      /* UG CC */
      /* mi: 70 mh: -10 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                           [nn_scores_bp_2_idx (c, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (c, c, sigma)] = 30;   

      /* UG CU */
      /* mi: 70 mh: 0 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                           [nn_scores_bp_2_idx (c, u, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (c, u, sigma)] = 35;   

      /* UG GA */
      /* mi: -40 mh: -80 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                           [nn_scores_bp_2_idx (g, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (g, a, sigma)] = -60; 

      /* UG GG */
      /* mi: 70 mh: -30 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                           [nn_scores_bp_2_idx (g, g, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (g, g, sigma)] = 20; 

      /* UG UC */
      /* mi: 70 mh: -10 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                           [nn_scores_bp_2_idx (u, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (u, c, sigma)] = 30;

      /* UG UU */
      /* mi: 0 mh: -80 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                           [nn_scores_bp_2_idx (u, u, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, g, sigma)]
                   [nn_scores_bp_2_idx (u, u, sigma)] = -40;

      /* AU AA */
      /* mi: 70 mh: -30 */
      assert (this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                           [nn_scores_bp_2_idx (a, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (a, a, sigma)] = 20;

      /* AU AC */
      /* mi: 70 mh: -50 */
      assert (this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                           [nn_scores_bp_2_idx (a, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (a, c, sigma)] = 10;

      /* AU AG */
      /* mi: -40 mh: -30 */
      assert (this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                           [nn_scores_bp_2_idx (a, g, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (a, g, sigma)] = -35;

      /* AU CA */
      /* mi: 70 mh: -10 */
      assert (this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                           [nn_scores_bp_2_idx (c, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (c, a, sigma)] = 30;

      /* AU CC */
      /* mi: 70 mh: -20 */
      assert (this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                           [nn_scores_bp_2_idx (c, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (c, c, sigma)] = 25;

      /* AU CU */
      /* mi: 70 mh: -20 */
      assert (this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                           [nn_scores_bp_2_idx (c, u, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (c, u, sigma)] = 25;

      /* AU GA */
      /* mi: -40 mh: -110 */
      assert (this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                           [nn_scores_bp_2_idx (g, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (g, a, sigma)] = -75;

      /* AU GG */
      /* mi: 70 mh: -20 */
      assert (this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                           [nn_scores_bp_2_idx (g, g, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (g, g, sigma)] = 25;

      /* AU UC */
      /* mi: 70 mh: -30 */
      assert (this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                           [nn_scores_bp_2_idx (u, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (u, c, sigma)] = 20;

      /* AU UU */
      /* mi: 0 mh: -110 */
      assert (this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                           [nn_scores_bp_2_idx (u, u, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (a, u, sigma)]
                   [nn_scores_bp_2_idx (u, u, sigma)] = -55;

      /* UA AA */
      /* mi: 70 mh: -50 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                           [nn_scores_bp_2_idx (a, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (a, a, sigma)] = 10;

      /* UA AC */
      /* mi: 70 mh: -30 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                           [nn_scores_bp_2_idx (a, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (a, c, sigma)] = 20;

      /* UA AG */
      /* mi: -40 mh: -60 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                           [nn_scores_bp_2_idx (a, g, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (a, g, sigma)] = -50;

      /* UA CA */
      /* mi: 70 mh: -20 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                           [nn_scores_bp_2_idx (c, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (c, a, sigma)] = 25;

      /* UA CC */
      /* mi: 70 mh: -10 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                           [nn_scores_bp_2_idx (c, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (c, c, sigma)] = 30;

      /* UA CU */
      /* mi: 70 mh: 0 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                           [nn_scores_bp_2_idx (c, u, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (c, u, sigma)] = 35;

      /* UA GA */
      /* mi: -40 mh: -140 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                           [nn_scores_bp_2_idx (g, a, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (g, a, sigma)] = -90;

      /* UA GG */
      /* mi: 70 mh: -70 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                           [nn_scores_bp_2_idx (g, g, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (g, g, sigma)] = 0;

      /* UA UC */
      /* mi: 70 mh: -10 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                           [nn_scores_bp_2_idx (u, c, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (u, c, sigma)] = 30;

      /* UA UU */
      /* mi: 0 mh: -80 */
      assert (this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                           [nn_scores_bp_2_idx (u, u, sigma)] == 0);
      this->G_stack[nn_scores_bp_2_idx (u, a, sigma)]
                   [nn_scores_bp_2_idx (u, u, sigma)] = -40;
   }

   return this;
}

/** @brief Delete a Nearest Neighbour scoring scheme.
 *
 * The destructor for @c NN_scores objects.
 *
 * @param[in] this object to be freed.
 */
void
nn_scores_delete (NN_scores* this)
{
   if (this != NULL)
   {
     XFREE_2D ((void**)this->G_stack);
     XFREE (this);
   }
}


/*********************************   Output   *********************************/

/** @brief Print the stacking energies of a scoring scheme to a stream.
 *
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] scheme The scoring scheme.
 * @params[in] sigma Alphabet of the scheme.
 */
void
nn_scores_fprintf_G_stack (FILE* stream,
                           const NN_scores* scheme,
                           const Alphabet* sigma)
{
   unsigned long i, j, k, l;
   long tmp;
   int rprec;
   unsigned long pline_width = 0;
   unsigned long matrix_width;
   unsigned long alpha_size = alphabet_size (sigma);
   char* string;
   char* string_start;

   assert (scheme != NULL);
   assert (scheme->G_stack != NULL);
   assert (sigma != NULL);

   matrix_width = alpha_size * alpha_size;

   /* dermine widest cell */
   for (i = 0; i < matrix_width; i++)
   {
      for (j = 0; j < matrix_width; j++)
      {
         rprec = 0;
         tmp = scheme->G_stack[i][j];
         if (tmp < 0)
         {
            tmp *= (-1);
            rprec++;
         }
         
         if (tmp > 0)
         {
            rprec += floor (log10 (tmp) + 1);
         }
         else
         {
            rprec += 1;
         }

         if ((unsigned) rprec > pline_width) 
         {
            pline_width = rprec;
         } 
      } 
   }
   rprec = pline_width;

   /* add up components of a line */
   pline_width += 3;            /*\s|\s*/
   pline_width *= matrix_width;
   pline_width += 2;            /* + NN */
   pline_width += 1;            /* + \n */

   /* alloc memory for the string */
   string = XMALLOC (sizeof (char) * ((pline_width * (matrix_width + 1)) + 1));
   if (string == NULL)
   {
      return;
   }

   string_start = string;
   
   /* print base pairs horizontally */
   msprintf (string, "  "); /*skip first "NN"*/
   string += 2;

   if (rprec > 0)
   {
      rprec--; /* lower cell width because we print 2 nucleotides */
   }

   for (i = 0; i < alpha_size; i++)
   {
      for (j = 0; j < alpha_size; j++)
      {
         msprintf (string, " | %*c%c", rprec, alphabet_no_2_base (i, sigma),
                                              alphabet_no_2_base (j, sigma));
         string += (rprec + 4);
      }
   }
   string[0] = '\n';
   string++;
   rprec++; /* restore cell width */

   /* print matrix */
   for (i = 0; i < alpha_size; i++)
   {   
      for (j = 0; j < alpha_size; j++)
      {
         /* column id */
         msprintf (string, "%c%c", alphabet_no_2_base (i, sigma),
                                   alphabet_no_2_base (j, sigma));
         string += 2;

         /* values for (i,j) (k, l) */
         for (k = 0; k < alpha_size; k++)
         {   
            for (l = 0; l < alpha_size; l++)
            {
               msprintf (string, " | ");
               string += 3;

               msprintf (string, "%*ld", rprec,
                         scheme->G_stack[nn_scores_bp_2_idx (i, j, sigma)]
                                        [nn_scores_bp_2_idx (k, l, sigma)]);
               string += rprec;
            }
         }         
         string[0] = '\n';
         string++;
      }
   }
  
   /* print matrix */
   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   XFREE (string_start);
}


/******************************   Miscellaneous   *****************************/

/** @brief Create a new Nearest Neighbour scoring scheme with standard values.
 *
 * The constructor for an initialised @c NN_scores objects. If compiled with
 * enabled memory checking, @c file and @c line should point to the position
 * where the function was called. Both parameters are automatically set by
 * using the macro @c NN_SCORES_NEW_INIT.\n
 * Returns @c NULL on error.
 *
 * @param[in] sigma alphabet.
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
unsigned long
nn_scores_bp_2_idx (const char base1, const char base2, const Alphabet* sigma)
{
   assert (sigma != NULL);

   /*mprintf ("(%c,%c) = %lu\n", alphabet_no_2_base(base1, sigma),
                               alphabet_no_2_base(base2, sigma),
                               base1 * alphabet_size (sigma) + base2);*/

   return base1 * alphabet_size (sigma) + base2;
}
