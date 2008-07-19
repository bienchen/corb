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

/* no. of canonical base pairs + whobble GU */
#define NO_ALLOWED_BP 6

struct NN_scores {
      long** G_stack;                /* stacking energies */
      long** G_mm_stack;             /* stacks with one mismatch */
      unsigned long G_stack_size;
      unsigned long G_mm_stack_size;
      char** bp_allowed;             /* canonical base pairs + whobble GU */
      unsigned long bp_allowed_size;
      char** bp_idx;                 /* indeces for base pairs */
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
      this->G_stack         = NULL;
      this->G_stack_size    = 0;
      this->G_mm_stack      = NULL;
      this->G_mm_stack_size = 0;
      this->bp_idx          = NULL;
      this->bp_allowed      = NULL;
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

      /* create list of allowed base pairs */
      this->bp_allowed_size = NO_ALLOWED_BP;
      this->bp_allowed = (char**) XMALLOC_2D (this->bp_allowed_size, 2,
                                              sizeof (char));
      if (this->bp_allowed == NULL)
      {
         nn_scores_delete (this);
         return NULL;        
      }

      a = alphabet_base_2_no('A', sigma);
      u = alphabet_base_2_no('U', sigma);
      g = alphabet_base_2_no('G', sigma);
      c = alphabet_base_2_no('C', sigma);

      this->bp_allowed[0][0] = c; this->bp_allowed[0][1] = g; /* CG */
      this->bp_allowed[1][0] = g; this->bp_allowed[1][1] = c; /* GC */
      this->bp_allowed[2][0] = g; this->bp_allowed[2][1] = u; /* GU */
      this->bp_allowed[3][0] = u; this->bp_allowed[3][1] = g; /* UG */
      this->bp_allowed[4][0] = a; this->bp_allowed[4][1] = u; /* AU */
      this->bp_allowed[5][0] = u; this->bp_allowed[5][1] = a; /* UA */

      /* create base pair indeces */
      this->bp_idx = (char**) XMALLOC_2D (alphabet_size (sigma),
                                              alphabet_size (sigma),
                                              sizeof (char));
      if (this->bp_idx == NULL)
      {
         nn_scores_delete (this);
         return NULL;        
      }

      for (i = 0; i < this->bp_allowed_size; i++)
      {
         this->bp_idx[(int) this->bp_allowed[i][0]]
                     [(int) this->bp_allowed[i][1]] = i;
      }

      this->bp_idx[(int) a][(int) a] = i; /* AA */
      i++;
      this->bp_idx[(int) a][(int) g] = i; /* AG */
      i++;
      this->bp_idx[(int) a][(int) c] = i; /* AC */
      i++;
      this->bp_idx[(int) u][(int) u] = i; /* UU */
      i++;
      this->bp_idx[(int) u][(int) c] = i; /* UC */
      i++;
      this->bp_idx[(int) g][(int) a] = i; /* GA */
      i++;
      this->bp_idx[(int) g][(int) g] = i; /* GG */
      i++;
      this->bp_idx[(int) c][(int) a] = i; /* CA */
      i++;
      this->bp_idx[(int) c][(int) u] = i; /* CU */
      i++;
      this->bp_idx[(int) c][(int) c] = i; /* CC */
      i++;

      /* prepare table for stacking energies */
      this->G_stack_size = this->bp_allowed_size;

      this->G_stack = (long**) XMALLOC_2D (this->G_stack_size,
                                           this->G_stack_size,
                                           sizeof (long));
      this->G_stack_size *= this->G_stack_size;

      if (this->G_stack == NULL)
      {
         nn_scores_delete (this);
         return NULL;        
      }

      for (i = 0; i < this->G_stack_size; i++)
      {
         this->G_stack[0][i] = 0;
      }

      /* set stacking energies (DG) */
      /* 5' - ip - 3' */
      /* 3' - jq - 5' --> i < p < q < j */
      /*    = ij qp    */

      /* regular pairs */
      /* AU AU */
      /* 5'- AU 
             UA -5' */
      this->G_stack[(int) this->bp_idx[(int)a][(int)u]]
                   [(int) this->bp_idx[(int)a][(int)u]] = -110;

      /* AU UA */
      /* 5'- AA
             UU -5' */
      this->G_stack[(int) this->bp_idx[(int)a][(int)u]]
                   [(int) this->bp_idx[(int)u][(int)a]] = -90;

      /* AU UG */
      /* 5'- AG
             UU -5' */
      this->G_stack[(int) this->bp_idx[(int)a][(int)u]]
                   [(int) this->bp_idx[(int)u][(int)g]] = -60;
     

      /* AU GU */
      /* 5'- AU
             UG -5'*/
      this->G_stack[(int) this->bp_idx[(int)a][(int)u]]
                   [(int) this->bp_idx[(int)g][(int)u]] = -140;      

      /* AU CG */
      /* 5'- AG
             UC -5'*/
      this->G_stack[(int) this->bp_idx[(int)a][(int)u]]
                   [(int) this->bp_idx[(int)c][(int)g]] = -210;        

      /* AU GC */
      /* 5'- AC
             UG -5'*/
      this->G_stack[(int) this->bp_idx[(int)a][(int)u]]
                   [(int) this->bp_idx[(int)g][(int)c]] = -220; 

      /* UA AU */
      /* 5'- UU 
             AA  -5'*/
      this->G_stack[(int) this->bp_idx[(int)u][(int)a]]
                   [(int) this->bp_idx[(int)a][(int)u]] = -90;

      /* UA UA */
      /* 5'- UA
             AU -5'*/
      this->G_stack[(int) this->bp_idx[(int)u][(int)a]]
                   [(int) this->bp_idx[(int)u][(int)a]] = -130;

      /* UA UG */
      /* 5'- UG
             AU -5'*/
      this->G_stack[(int) this->bp_idx[(int)u][(int)a]]
                   [(int) this->bp_idx[(int)u][(int)g]] = -100;

      /* UA GU */
      /* 5'- UU
             AG -5'*/
      this->G_stack[(int) this->bp_idx[(int)u][(int)a]]
                   [(int) this->bp_idx[(int)g][(int)u]] = -130;

      /* UA CG */
      /* 5'- UG
             AC -5'*/
      this->G_stack[(int) this->bp_idx[(int)u][(int)a]]
                   [(int) this->bp_idx[(int)c][(int)g]] = -210;

      /* UA GC */
      /* 5'- UC
             AG -5'*/
      this->G_stack[(int) this->bp_idx[(int)u][(int)a]]
                   [(int) this->bp_idx[(int)g][(int)c]] = -240;
    
      /* UG AU */
      /* 5'- UU
             GA -5'*/
      this->G_stack[(int) this->bp_idx[(int)u][(int)g]]
                   [(int) this->bp_idx[(int)a][(int)u]] = -60;

      /* UG UA */
      /* 5'- UA 
             GU -5'*/
      this->G_stack[(int) this->bp_idx[(int)u][(int)g]]
                   [(int) this->bp_idx[(int)u][(int)a]] = -100;

      /* UG UG */
      /* 5'- UG 
             GU -5'*/
      this->G_stack[(int) this->bp_idx[(int)u][(int)g]]
                   [(int) this->bp_idx[(int)u][(int)g]] = 30;

      /* UG GU*/
      /* 5'- UU 
             GG -5'*/
      this->G_stack[(int) this->bp_idx[(int)u][(int)g]]
                   [(int) this->bp_idx[(int)g][(int)u]] = -50;

      /* UG CG */
      /* 5'- UG 
             GC -5'*/
      this->G_stack[(int) this->bp_idx[(int)u][(int)g]]
                   [(int) this->bp_idx[(int)c][(int)g]] = -140;

      /* UG GC */
      /* 5'- UC 
             GG -5'*/
      this->G_stack[(int) this->bp_idx[(int)u][(int)g]]
                   [(int) this->bp_idx[(int)g][(int)c]] = -150;

      /* GU AU */
      /* 5'- GU 
             UA -5'*/
      this->G_stack[(int) this->bp_idx[(int)g][(int)u]]
                   [(int) this->bp_idx[(int)a][(int)u]] = -140;

      /* GU UA */
      /* 5'- GA 
             UU -5'*/
      this->G_stack[(int) this->bp_idx[(int)g][(int)u]]
                   [(int) this->bp_idx[(int)u][(int)a]] = -130;

      /* GU UG */
      /* 5'- GG 
             UU -5'*/
      this->G_stack[(int) this->bp_idx[(int)g][(int)u]]
                   [(int) this->bp_idx[(int)u][(int)g]] = -50;

      /* GU GU */
      /* 5'- GU 
             UG -5'*/
      this->G_stack[(int) this->bp_idx[(int)g][(int)u]]
                   [(int) this->bp_idx[(int)g][(int)u]] = 130;

      /* GU CG */
      /* 5'- GG 
             UC -5'*/
      this->G_stack[(int) this->bp_idx[(int)g][(int)u]]
                   [(int) this->bp_idx[(int)c][(int)g]] = -210;

      /* GU GC */
      /* 5'- GC 
             UG -5'*/
      this->G_stack[(int) this->bp_idx[(int)g][(int)u]]
                   [(int) this->bp_idx[(int)g][(int)c]] = -250;

      /* CG AU */
      /* 5'- C
             G -5'*/
      this->G_stack[(int) this->bp_idx[(int)c][(int)g]]
                   [(int) this->bp_idx[(int)a][(int)u]] = -210;

      /* CG UA */
      /* 5'- CA
             GU -5'*/
      this->G_stack[(int) this->bp_idx[(int)c][(int)g]]
                   [(int) this->bp_idx[(int)u][(int)a]] = -210;

      /* CG UG */
      /* 5'- CG
             GU -5'*/
      this->G_stack[(int) this->bp_idx[(int)c][(int)g]]
                   [(int) this->bp_idx[(int)u][(int)g]] = -140;

      /* CG GU */
      /* 5'- CU
             GG -5'*/
      this->G_stack[(int) this->bp_idx[(int)c][(int)g]]
                   [(int) this->bp_idx[(int)g][(int)u]] = -210;

      /* CG CG */
      /* 5'- CG
             GC -5'*/
      this->G_stack[(int) this->bp_idx[(int)c][(int)g]]
                   [(int) this->bp_idx[(int)c][(int)g]] = -240;

      /* CG GC */
      /* 5'- CC
             GG -5'*/
      this->G_stack[(int) this->bp_idx[(int)c][(int)g]]
                   [(int) this->bp_idx[(int)g][(int)c]] = -330;

      /* GC AU */
      /* 5'- GU
             CA -5'*/
      this->G_stack[(int) this->bp_idx[(int)g][(int)c]]
                   [(int) this->bp_idx[(int)a][(int)u]] = -220;

      /* GC UA */
      /* 5'- GA
             CU -5'*/
      this->G_stack[(int) this->bp_idx[(int)g][(int)c]]
                   [(int) this->bp_idx[(int)u][(int)a]] = -240;

      /* GC UG */
      /* 5'- GG
             CU -5'*/
      this->G_stack[(int) this->bp_idx[(int)g][(int)c]]
                   [(int) this->bp_idx[(int)u][(int)g]] = -150;

      /* GC GU */
      /* 5'- GU
             CG -5'*/
      this->G_stack[(int) this->bp_idx[(int)g][(int)c]]
                   [(int) this->bp_idx[(int)g][(int)u]] = -250;

      /* GC CG */
      /* 5'- GG
             CC -5'*/
      this->G_stack[(int) this->bp_idx[(int)g][(int)c]]
                   [(int) this->bp_idx[(int)c][(int)g]] = -330;

      /* GC GC */
      /* 5'- GC
             CG -5'*/
      this->G_stack[(int) this->bp_idx[(int)g][(int)c]]
                   [(int) this->bp_idx[(int)g][(int)c]] = -340;

      /* prepare table for mismatch stacking energies */
      this->G_mm_stack_size = alphabet_size (sigma) * alphabet_size (sigma);

      this->G_mm_stack = (long**) XMALLOC_2D (this->bp_allowed_size,
                                              this->G_mm_stack_size,
                                              sizeof (long));
      this->G_mm_stack_size *= this->bp_allowed_size;

      if (this->G_mm_stack == NULL)
      {
         nn_scores_delete (this);
         return NULL;        
      }

      for (i = 0; i < this->G_mm_stack_size; i++)
      {
         this->G_mm_stack[0][i] = 0;
      }

      /* stacks containing a mismatch */
      /* mi: param from mismatch_interior table */
      /* mh: param from mismatch_hairpin */

      /* AU AA */
      /* mi: 70 mh: -30 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)a][(int)a]] = 20;
      /* AU AU */
      /* mi: 70 mh: -30 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)a][(int)u]] = 20;
      /* AU AG */
      /* mi: -40 mh: -30 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)a][(int)g]] = -35;
      /* AU AC */
      /* mi: 70 mh: -50 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)a][(int)c]] = 10;

      /* AU UA */
      /* mi: 70 mh: -30 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)u][(int)a]] = 20;
      /* AU UU */
      /* mi: 0 mh: -110 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)u][(int)u]] = -55;
      /* AU UG */
      /* mi: 70 mh: -60 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)u][(int)g]] = 5;
      /* AU UC */
      /* mi: 70 mh: -30 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)u][(int)c]] = 20;

      /* AU GA */
      /* mi: -40 mh: -110 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)g][(int)a]] = -75;
      /* AU GU */
      /* mi: 70 mh: 20 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)g][(int)u]] = 45;
      /* AU GG */
      /* mi: 70 mh: -20 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)g][(int)g]] = 25;
      /* AU GC */
      /* mi: 70 mh: -120 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)g][(int)c]] = -25;

      /* AU CA */
      /* mi: 70 mh: -10 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)c][(int)a]] = 30;
      /* AU CU */
      /* mi: 70 mh: -20 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)c][(int)u]] = 25;
      /* AU CG */
      /* mi: 70 mh: -150 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)c][(int)g]] = -40;
      /* AU CC */
      /* mi: 70 mh: -20 */
      this->G_mm_stack[(int) this->bp_idx[(int)a][(int)u]]
                      [(int) this->bp_idx[(int)c][(int)c]] = 25;


      /* UA AA */
      /* mi: 70 mh: -50 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)a][(int)a]] = 10;
      /* UA AU */
      /* mi: 70 mh: -50 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)a][(int)u]] = 10;
      /* UA AG */
      /* mi: -40 mh: -60 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)a][(int)g]] = -50;
      /* UA AC */
      /* mi: 70 mh: -30 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)a][(int)c]] = 20;

      /* UA UA */
      /* mi: 70 mh: -30 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)u][(int)a]] = 20;
      /* UA UU */
      /* mi: 0 mh: -80 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)u][(int)u]] = -40;
      /* UA UG */
      /* mi: 70 mh: -50 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)u][(int)g]] = 10;
      /* UA UC */
      /* mi: 70 mh: -10 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)u][(int)c]] = 30;

      /* UA GA */
      /* mi: -40 mh: -140 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)g][(int)a]] = -90;
      /* UA GU */
      /* mi: 70 mh: -20 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)g][(int)u]] = 25;
      /* UA GG */
      /* mi: 70 mh: -70 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)g][(int)g]] = 0;
      /* UA GC */
      /*mi: 70 mh: -120 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)g][(int)c]] = -25;

      /* UA CA */
      /* mi: 70 mh: -20 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)c][(int)a]] = 25;
      /* UA CU */
      /* mi: 70 mh: 0 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)c][(int)u]] = 35;
      /* UA CG */
      /* mi: 70 mh: -120 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)c][(int)g]] = -25;
      /* UA CC */
      /* mi: 70 mh: -10 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)a]]
                      [(int) this->bp_idx[(int)c][(int)c]] = 30;

      /* UG AA */
      /* mi: 70 mh: -50 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)a][(int)a]] = 10;
      /* UG AU */
      /* mi: 70 mh: -50 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)a][(int)u]] = 10;
      /* UG AG */
      /* mi: -40 mh: -60 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)a][(int)g]] = -50;
      /* UG AC */
      /* mi: 70 mh: -30 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)a][(int)c]] = 20;

      /* UG UA */
      /* mi: 70 mh: -60 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)u][(int)a]] = 5;
      /* UG UU */
      /* mi: 0 mh: -80 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)u][(int)u]] = -40;
      /* UG UG */
      /* mi: 70 mh: -60 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)u][(int)g]] = 5;
      /* UG UC */
      /* mi: 70 mh: -10 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)u][(int)c]] = 30;

      /* UG GA */
      /* mi: -40 mh: -80 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)g][(int)a]] = -60;
      /* UG GU */
      /* mi: 70 mh: -70 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)g][(int)u]] = 0;
      /* UG GG */
      /* mi: 70 mh: -30 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)g][(int)g]] = 20;
      /* UG GC */
      /* mi: 70 mh: -120 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)g][(int)c]] = -25;

      /* UG CA */
      /* mi: 70 mh: -20 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)c][(int)a]] = 25;
      /* UG CU */
      /* mi: 70 mh: 0 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)c][(int)u]] = 35;
      /* UG CG */
      /* mi: 70 mh: -170 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)c][(int)g]] = -50;
      /* UG CC */
      /* mi: 70 mh: -10 */
      this->G_mm_stack[(int) this->bp_idx[(int)u][(int)g]]
                      [(int) this->bp_idx[(int)c][(int)c]] = 30;


      /* GU AA */
      /* mi: 70 mh: 20 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)a][(int)a]] = 45;
      /* GU AU */
      /* mi: 70 mh: -30 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)a][(int)u]] = 20;
      /* GU AG */
      /* mi: -40 mh: -30 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)a][(int)g]] = -35;
      /* GU AC */
      /* mi: 70 mh: -50 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)a][(int)c]] = 10;

      /* GU UA */
      /* mi: 70 mh: -30 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)u][(int)a]] = 20;
      /* GU UU */
      /* mi: 0 mh: -110 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)u][(int)u]] = -55;
      /* GU UG */
      /* mi: 70 mh: -30 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)u][(int)g]] = 20;
      /* GU UC */
      /* mi: 70 mh: -30 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)u][(int)c]] = 20;

      /* GU GA */
      /* mi: -40 mh: -90 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)g][(int)a]] = -65;
      /* GU GU */
      /* mi: 70 mh: 0 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)g][(int)u]] = 35;
      /* GU GG */
      /* mi: 70 mh: -30 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)g][(int)g]] = 20;
      /* GU GC */
      /* mi: 70 mh: -110*/
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)g][(int)c]] = -20;

      /* GU CA */
      /* mi: 70 mh: -10 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)c][(int)a]] = 30;
      /* GU CU */
      /* mi: 70 mh: -20 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)c][(int)u]] = 25;
      /* GU CG */
      /* mi: 70 mh: -150 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)c][(int)g]] = -40;
      /* GU CC */
      /* mi: 70 mh: -20 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)u]]
                      [(int) this->bp_idx[(int)c][(int)c]] = 25;


      /* GC AA */
      /* mi: 0 mh: -110 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)a][(int)a]] = -55;
      /* GC AU */
      /* mi: 0 mh: -210 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)a][(int)u]] = -105;
      /* GC AG */
      /* mi: -110 mh: -130 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)a][(int)g]] = -120;
      /* GC AC */
      /* mi: 0 mh: -150 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)a][(int)c]] = -75;

      /* GC UA */
      /* mi: 0 mh: -190 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)u][(int)a]] = -145;
      /* GC UU */
      /* mi: -70 mh: -150 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)u][(int)u]] = -110;
      /* GC UG */
      /* mi: 0 mh: -220 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)u][(int)g]] = -110;
      /* GC UC */
      /* mi: 0 mh: -100 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)u][(int)c]] = -50;

      /* GC GA */
      /* mi: -110 mh: -240 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)g][(int)a]] = -175;
      /* GC GU */
      /* mi: 0 mh: -120 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)g][(int)u]] = -60;
      /* GC GG */
      /* mi: 0 mh: -140 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)g][(int)g]] = -70;
      /* GC GC */
      /* mi: 0 mh: -290 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)g][(int)c]] = -145;

      /* GC CA */
      /* mi: 0 mh: -110 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)c][(int)a]] = -55;
      /* GC CU */
      /* mi: 0 mh: -50 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)c][(int)u]] = -25;
      /* GC CG */
      /* mi: 0 mh: -240 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)c][(int)g]] = -120;
      /* GC CC */
      /* mi: 0 mh: -70 */
      this->G_mm_stack[(int) this->bp_idx[(int)g][(int)c]]
                      [(int) this->bp_idx[(int)c][(int)c]] = -35;


      /* CG AA */
      /* mi: 0 mh: -150 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                      [(int) this->bp_idx[(int)a][(int)a]] = -75;
      /* CG AU */
      /* mi: 0 mh: -180 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                      [(int) this->bp_idx[(int)a][(int)u]] = -90;
      /* CG AG */
      /* mi: -110 mh: -140 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                      [(int) this->bp_idx[(int)a][(int)g]] = -125;
      /* CG AC */
      /* mi: 0 mh: -150 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                      [(int) this->bp_idx[(int)a][(int)c]] = -75;

      /* CG UA */
      /* mi: 0 mh: -170 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                      [(int) this->bp_idx[(int)u][(int)a]] = -85;
      /* CG UU */
      /* mi: -70 mh: -200 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                      [(int) this->bp_idx[(int)u][(int)u]] = -135;
      /* CG UG */
      /* mi: 0 mh: -180 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                      [(int) this->bp_idx[(int)u][(int)g]] = -90;
      /* CG UC */
      /* mi: 0 mh: -140 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                      [(int) this->bp_idx[(int)u][(int)c]] = -70;

      /* CG GA */
      /* mi: -110 mh: -220 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                      [(int) this->bp_idx[(int)g][(int)a]] = -165;
      /* CG GU */
      /* mi: 0 mh: -110 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                      [(int) this->bp_idx[(int)g][(int)u]] = -55;
      /* CG GG */
      /* mi: 0 mh: -160 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                      [(int) this->bp_idx[(int)g][(int)g]] = -80;
      /* CG GC */
      /* mi: 0 mh: -200 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                      [(int) this->bp_idx[(int)g][(int)c]] = -100;

      /* CG CA */
      /* mi: 0 mh: -100 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                   [(int) this->bp_idx[(int)c][(int)a]] = -50;
      /* CG CU */
      /* mi: 0 mh: -80 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                      [(int) this->bp_idx[(int)c][(int)u]] = -40;
      /* CG CG */
      /* mi: 0 mh: -290 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                      [(int) this->bp_idx[(int)c][(int)g]] = -145;

      /* CG CC */
      /* mi: 0 mh: -90 */
      this->G_mm_stack[(int) this->bp_idx[(int)c][(int)g]]
                      [(int) this->bp_idx[(int)c][(int)c]] = -45;
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
     XFREE_2D ((void**)this->G_mm_stack);
     XFREE_2D ((void**)this->bp_idx);
     XFREE_2D ((void**)this->bp_allowed);
     XFREE (this);
   }
}


/*********************************   Access   *********************************/

/** @brief Return number of allowed base pairs in a schoring scheme.
 *
 * @params[in]  i position of base pair.
 * @params[out] b5 5' base of the pair.
 * @params[out] b3 3' base of the pair.
 * @params[in] scheme The scoring scheme.
 */
void
nn_scores_get_allowed_basepair (const unsigned i,
                                char* b5,
                                char* b3,
                                const NN_scores* scheme)
{
   assert (scheme);
   assert (i < scheme->bp_allowed_size);

   b5[0] = scheme->bp_allowed[i][0];
   b3[0] = scheme->bp_allowed[i][1];
}


/*********************************    Size    *********************************/

/** @brief Return number of allowed base pairs in a schoring scheme.
 *
 * @params[in] scheme The scoring scheme.
 */
unsigned long
nn_scores_no_allowed_basepairs (const NN_scores* scheme)
{
   assert (scheme);

   return scheme->bp_allowed_size;
}

/*********************************   Output   *********************************/

/** @brief Print the allowed base pairs of a scoring scheme to a stream.
 *
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] scheme The scoring scheme.
 * @params[in] sigma Alphabet of the scheme.
 */
void
nn_scores_fprintf_bp_allowed (FILE* stream,
                              const NN_scores* scheme,
                              const Alphabet* sigma)
{
   unsigned long i;
   char* string;
   char* string_start;

   assert (scheme != NULL);
   assert (scheme->bp_allowed != NULL);
   assert (sigma != NULL);

   /* alloc memory for the string */
   string = XMALLOC (sizeof (char) * ((scheme->bp_allowed_size * 3) + 1));
   if (string == NULL)
   {
      return;
   }   
   string_start = string;

   for (i = 0; i < scheme->bp_allowed_size; i++)
   {
      msprintf (string, "%c%c\n",
                alphabet_no_2_base (scheme->bp_allowed[i][0], sigma),
                alphabet_no_2_base (scheme->bp_allowed[i][1], sigma));
      string += 3;
   }

   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   XFREE (string_start);
}

/** @brief Print the indeces of base pairs of a scoring scheme to a stream.
 *
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] scheme The scoring scheme.
 * @params[in] sigma Alphabet of the scheme.
 */
void
nn_scores_fprintf_bp_idx (FILE* stream,
                          const NN_scores* scheme,
                          const Alphabet* sigma)
{
   unsigned long pline_width = 0;
   char* string;
   char* string_start;
   unsigned long alpha_size = alphabet_size (sigma);
   unsigned long i, j;
   int rprec;

   assert (scheme != NULL);
   assert (scheme->bp_idx != NULL);
   assert (sigma != NULL);

   /* widest cell can be determined by the size of the alphabet */
   rprec = alpha_size * alpha_size;
   rprec = floor (log10 (rprec) + 1);
   
   /* calculate linewidth */
   pline_width = rprec;
   pline_width += 3;            /* \s|\s */
   pline_width *= alpha_size;
   pline_width += 2;            /* N\n */

   string = XMALLOC (sizeof (char) *
                     ((pline_width * (alpha_size + 1)) + 1));
   if (string == NULL)
   {
      return;
   }

   string_start = string;

   /* print base pairs horizontally */
   msprintf (string, " "); /*skip first "N"*/
   string += 1;

   for (i = 0; i < alpha_size; i++)
   {
      msprintf (string, " | %*c", rprec, alphabet_no_2_base (i, sigma));
      string += (rprec + 3);      
   }
   string[0] = '\n';
   string++;

   /* print index lines */
   for (i = 0; i < alpha_size; i++)
   {
      msprintf (string, "%c", alphabet_no_2_base (i, sigma));
      string++;
      
      for (j = 0; j < alpha_size; j++)
      {
         msprintf (string, " | %*d", rprec, scheme->bp_idx[i][j]);
         string += (rprec + 3);         
      }

      string[0] = '\n';
      string++;
   }

   /* print matrix */
   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   XFREE (string_start);
}

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
   unsigned long i, j;
   int d5, u5, d3, u3;
   long tmp;
   int rprec;
   unsigned long pline_width = 0;
   unsigned long matrix_edge;
   char* string;
   char* string_start;

   assert (scheme != NULL);
   assert (scheme->G_stack != NULL);
   assert (scheme->bp_allowed != NULL);
   assert (sigma != NULL);

   matrix_edge = scheme->bp_allowed_size;

   /* dermine widest cell */
   for (i = 0; i < matrix_edge; i++)
   {
      for (j = 0; j < matrix_edge; j++)
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
   pline_width *= matrix_edge;
   pline_width += 2;            /* + NN */
   pline_width += 1;            /* + \n */

   /* alloc memory for the string */
   string = XMALLOC (sizeof (char) * ((pline_width * (matrix_edge + 1)) + 1));
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

   for (i = 0; i < scheme->bp_allowed_size; i++)
   {
      msprintf (string, " | %*c%c", rprec,
                alphabet_no_2_base (scheme->bp_allowed[i][0], sigma),
                alphabet_no_2_base (scheme->bp_allowed[i][1], sigma));
      string += (rprec + 4);
   }

   string[0] = '\n';
   string++;
   rprec++; /* restore cell width */

   /* print matrix */
   for (i = 0; i < scheme->bp_allowed_size; i++)
   {   
      /* column id */
      msprintf (string, "%c%c", 
                alphabet_no_2_base (scheme->bp_allowed[i][0], sigma),
                alphabet_no_2_base (scheme->bp_allowed[i][1], sigma));
      string += 2;

      u5 = scheme->bp_allowed[i][0];
      d5 = scheme->bp_allowed[i][1];
      
      /* values for (i,j) */
      for (j = 0; j < scheme->bp_allowed_size; j++)
      {   
            msprintf (string, " | ");
            string += 3;

            u3 = scheme->bp_allowed[j][0];
            d3 = scheme->bp_allowed[j][1];

            msprintf (string, "%*ld", rprec,
                      scheme->G_stack[(int) scheme->bp_idx[u5][d5]]
                                     [(int) scheme->bp_idx[u3][d3]]);
            string += rprec;
      }       
      string[0] = '\n';
      string++;
   }
  
   /* print matrix */
   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   XFREE (string_start);
}

/** @brief Print the mismatch stacking energies of a scoring scheme to a stream.
 *
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] scheme The scoring scheme.
 * @params[in] sigma Alphabet of the scheme.
 */
void
nn_scores_fprintf_mm_G_stack (FILE* stream,
                              const NN_scores* scheme,
                              const Alphabet* sigma)
{
   unsigned long i, j, k;
   int d, u;
   long tmp;
   int rprec;
   unsigned long pline_width = 2;
   unsigned long matrix_rows;
   unsigned long matrix_cols;
   char* string;
   char* string_start;

   assert (scheme != NULL);
   assert (scheme->G_mm_stack != NULL);
   assert (scheme->bp_allowed != NULL);
   assert (sigma != NULL);

   matrix_rows = scheme->bp_allowed_size;
   matrix_cols = scheme->G_mm_stack_size / scheme->bp_allowed_size;

   /* dermine widest cell */
   for (i = 0; i < matrix_rows; i++)
   {
      for (j = 0; j < matrix_cols; j++)
      {
         rprec = 0;
         tmp = scheme->G_mm_stack[i][j];
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
   pline_width *= matrix_cols;
   pline_width += 2;            /* + NN */
   pline_width += 1;            /* + \n */

   /* alloc memory for the string */
   string = XMALLOC (sizeof (char) * ((pline_width * (matrix_rows + 1)) + 1));
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

   for (i = 0; i < alphabet_size (sigma); i++)
   {
      for (j = 0; j < alphabet_size (sigma); j++)
      {
         msprintf (string, " | %*c%c", rprec,
                   alphabet_no_2_base (i, sigma),
                   alphabet_no_2_base (j, sigma));
         string += (rprec + 4);
      }
   }
   string[0] = '\n';
   string++;
   rprec++; /* restore cell width */

   /* print matrix */
   for (i = 0; i < scheme->bp_allowed_size; i++)
   {   
      /* column id */
      msprintf (string, "%c%c", 
                alphabet_no_2_base (scheme->bp_allowed[i][0], sigma),
                alphabet_no_2_base (scheme->bp_allowed[i][1], sigma));
      string += 2;

      u = scheme->bp_allowed[i][0];
      d = scheme->bp_allowed[i][1];
      
      /* values for (i,j) */
      for (j = 0; j < alphabet_size (sigma); j++)
      {   
         for (k = 0; k < alphabet_size (sigma); k++)
         {
            msprintf (string, " | ");
            string += 3;
            
            msprintf (string, "%*ld", rprec,
                      scheme->G_mm_stack[(int) scheme->bp_idx[u][d]]
                                     [(int) scheme->bp_idx[j][k]]);
            string += rprec;
         }
      }    
      string[0] = '\n';
      string++;
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
nn_scores_bp_2_idx (const char base1, const char base2, const NN_scores* scheme)
{
   assert (scheme != NULL);

   /*mprintf ("(%c,%c) = %lu\n", alphabet_no_2_base(base1, sigma),
                               alphabet_no_2_base(base2, sigma),
                               base1 * alphabet_size (sigma) + base2);*/

   return scheme->bp_idx[(int)base1][(int)base2];
}
