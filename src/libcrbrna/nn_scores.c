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
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <libcrbbasic/crbbasic.h>
/*#include "alphabet.h"*/
#include "nn_scores.h"

#define NO_ALLOWED_BP 6         /* no. of canonical base pairs + whobble GU */
#define D_MM_H 3                /* dimensions of the G_mismatch_hairpin table */
#define D_MM_I 3                /* dim. of the G_mismatch_interior table */
#define D_INT11 4               /* dimensionality of G_int11 table */
#define D_INT21 5               /* dimensionality of G_int21 table */
#define D_INT22 6               /* dimensionality of G_int22 table */
#define D_TL   6                /* size of a tetraloop + closing bp */
#define NN_LXC37 107.856        /* ask marco about this */
#define NN_ML_OFFSET 340
#define NN_ML_UNPAIRED 0
#define NN_ML_STEMS 40
#define NN_NINIO_M 50
#define NN_NINIO_MAX 300

struct NN_scores {
      long** G_stack;                        /* stacking energies */
      unsigned long G_stack_size;
      long** G_mm_stack;                     /* stacks with one mismatch */
      unsigned long G_mm_stack_size;
      int* G_hairpin_loop;                   /* hairpin loops */
      unsigned long G_hairpin_loop_size;
      int*** G_mismatch_hairpin;             /* hairpin loop closing bp */
      unsigned long G_mismatch_hairpin_size;
      int* non_gc_penalty_for_bp;            /* penalty for closing non-GC */
      char** tetra_loop;                     /* sorted list of possible loops */
      int* G_tetra_loop;                     /* scores */
      unsigned long tetra_loop_size;
      int* G_bulge_loop;                     /* bulge loops */
      unsigned long G_bulge_loop_size;
      /* internal loops */
      int* G_internal_loop;                  /* generic loops */
      unsigned long G_internal_loop_size;
      int**** G_int11;                       /* 1x1 loops */
      unsigned long G_int11_size;
      int***** G_int21;                      /* 2x1 loops */
      unsigned long G_int21_size;
      int****** G_int22;                      /* 2x1 loops */
      unsigned long G_int22_size;
      int*** G_mismatch_interior;             /* interior loop closing bp */
      unsigned long G_mismatch_interior_size;
      /* internal loops */
      int** G_dangle5;                       /* 5' dangling end, bp + base */
      unsigned long G_dangle5_size;
      int** G_dangle3;                       /* 3' dangling end, bp + base */
      unsigned long G_dangle3_size;
      char** bp_allowed;                     /* WC base pairs + whobble GU */
      unsigned long bp_allowed_size;
      char** bp_idx;                     /* indeces for base pairs */
      unsigned long bp_idx_size;
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
      this->G_stack                  = NULL;
      this->G_stack_size             = 0;
      this->G_mm_stack               = NULL;
      this->G_mm_stack_size          = 0;
      this->G_hairpin_loop           = NULL;
      this->G_hairpin_loop_size      = 0;
      this->G_mismatch_hairpin       = NULL;
      this->G_mismatch_hairpin_size  = 0;
      this->G_bulge_loop             = NULL;
      this->G_bulge_loop_size        = 0;
      this->non_gc_penalty_for_bp    = NULL;
      this->tetra_loop               = NULL;
      this->G_tetra_loop             = NULL;
      this->tetra_loop_size          = 0;
      this->G_internal_loop          = NULL;
      this->G_internal_loop_size     = 0;
      this->G_int11                  = NULL;
      this->G_int11_size             = 0;
      this->G_int21                  = NULL;
      this->G_int21_size             = 0;
      this->G_int22                  = NULL;
      this->G_int22_size             = 0;
      this->G_mismatch_interior      = NULL;
      this->G_mismatch_interior_size = 0;
      this->G_dangle5                = NULL;
      this->G_dangle5_size           = 0;
      this->G_dangle3                = NULL;
      this->G_dangle3_size           = 0;
      this->bp_idx                   = NULL;
      this->bp_allowed               = NULL;
      this->bp_allowed_size          = 0;
   }

   return this;
}

static int
allocate_init_bp_allowed (char a, char u, char g, char c,
                          NN_scores* this,
                          const char* file, const int line)
{
   this->bp_allowed_size = NO_ALLOWED_BP;
   this->bp_allowed = (char**) XOBJ_MALLOC_2D (this->bp_allowed_size, 2,
                                               sizeof (char),
                                               file, line);
   if (this->bp_allowed == NULL)
   {
      return 1;        
   }
   
   this->bp_allowed[0][0] = c; this->bp_allowed[0][1] = g; /* CG */
   this->bp_allowed[1][0] = g; this->bp_allowed[1][1] = c; /* GC */
   this->bp_allowed[2][0] = g; this->bp_allowed[2][1] = u; /* GU */
   this->bp_allowed[3][0] = u; this->bp_allowed[3][1] = g; /* UG */
   this->bp_allowed[4][0] = a; this->bp_allowed[4][1] = u; /* AU */
   this->bp_allowed[5][0] = u; this->bp_allowed[5][1] = a; /* UA */

   return 0;
}

static int
allocate_init_bp_idx (unsigned long size,
                      char a, char u, char g, char c, NN_scores* this,
                      const char* file, const int line)
{
   unsigned long i;

   this->bp_idx = (char**) XOBJ_MALLOC_2D (size, size, sizeof (char),
                                           file, line);
   if (this->bp_idx == NULL)
   {
      return 1;        
   }
   this->bp_idx_size = size * size;

   /* place indeces of allowed base pairs at the begining og the table */
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

   return 0;
}

static int
allocate_init_G_stack (char a, char u, char g, char c, NN_scores* this,
                       const char* file, const int line)
{
   /*unsigned long i;*/

   this->G_stack_size = this->bp_allowed_size;
   
   /* allocate matrix */
   this->G_stack = (long**) XOBJ_MALLOC_2D (this->G_stack_size,
                                            this->G_stack_size,
                                            sizeof (long),
                                            file, line);
   this->G_stack_size *= this->G_stack_size;
   
   if (this->G_stack == NULL)
   {
      return 1;
   }

   /*for (i = 0; i < this->G_stack_size; i++)
   {
      this->G_stack[0][i] = 0;
      }*/

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
   

   return 0;
}

static int
allocate_init_G_mm_stack_size (char a, char u, char g, char c,
                               unsigned long size, NN_scores* this,
                               const char* file, const int line)
{
   this->G_mm_stack_size = size * size;
   
   this->G_mm_stack = (long**) XOBJ_MALLOC_2D (this->bp_allowed_size,
                                               this->G_mm_stack_size,
                                               sizeof (long),
                                               file, line);
   this->G_mm_stack_size *= this->bp_allowed_size;

   if (this->G_mm_stack == NULL)
   {
      return 1;
   }
   
   /*for (i = 0; i < this->G_mm_stack_size; i++)
   {
      this->G_mm_stack[0][i] = 0;
      }*/

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

   return 0;
}

static int
allocate_init_hairpin_loop (NN_scores* this, const char* file, const int line)
{

   this->G_hairpin_loop_size = 31;
   this->G_hairpin_loop = XOBJ_MALLOC (  this->G_hairpin_loop_size
                                       * sizeof (this->G_hairpin_loop[0]),
                                         file, line);
   if (this->G_hairpin_loop == NULL)
   {
      return 1;        
   }

   this->G_hairpin_loop[ 0] = INT_UNDEF; /* min. loop length is 3 */
   this->G_hairpin_loop[ 1] = INT_UNDEF;
   this->G_hairpin_loop[ 2] = INT_UNDEF;
   this->G_hairpin_loop[ 3] = 570;
   this->G_hairpin_loop[ 4] = 560;
   this->G_hairpin_loop[ 5] = 560;
   this->G_hairpin_loop[ 6] = 540;
   this->G_hairpin_loop[ 7] = 590;
   this->G_hairpin_loop[ 8] = 560;
   this->G_hairpin_loop[ 9] = 640;
   this->G_hairpin_loop[10] = 650;
   this->G_hairpin_loop[11] = 660;
   this->G_hairpin_loop[12] = 670;
   this->G_hairpin_loop[13] = 678;
   this->G_hairpin_loop[14] = 686;
   this->G_hairpin_loop[15] = 694;
   this->G_hairpin_loop[16] = 701;
   this->G_hairpin_loop[17] = 707;
   this->G_hairpin_loop[18] = 713;
   this->G_hairpin_loop[19] = 719;
   this->G_hairpin_loop[20] = 725;
   this->G_hairpin_loop[21] = 730;
   this->G_hairpin_loop[22] = 735;
   this->G_hairpin_loop[23] = 740;
   this->G_hairpin_loop[24] = 744;
   this->G_hairpin_loop[25] = 749;
   this->G_hairpin_loop[26] = 753;
   this->G_hairpin_loop[27] = 757;
   this->G_hairpin_loop[28] = 761;
   this->G_hairpin_loop[29] = 765;
   this->G_hairpin_loop[30] = 769;

   return 0;
}

static int
allocate_init_mismatch_hairpin (int a, int u, int g, int c,
                                const unsigned long no_of_b,
                                NN_scores* this,
                                const char* file, const int line)
{
   /* allocate memory */
   this->G_mismatch_hairpin
      = (int***) XOBJ_MALLOC_ND(sizeof (***this->G_mismatch_hairpin),
                                D_MM_H,
                                file, line,
                                this->bp_allowed_size, no_of_b, no_of_b);
   if (this->G_mismatch_hairpin == NULL)
   {
      return 1;        
   }
   this->G_mismatch_hairpin_size = this->bp_allowed_size * no_of_b * no_of_b;

   /* store values this->bp_idx[][] */
   /* CG */
   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][a][a] = -150; /* AA */
   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][a][c] = -150; /* AC */
   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][a][g] = -140; /* AG */
   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][a][u] = -180; /* AU */

   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][c][a] = -100; /* CA */
   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][c][c] =  -90; /* CC */
   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][c][g] = -290; /* CG */
   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][c][u] =  -80; /* CU */

   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][g][a] = -220; /* GA */
   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][g][c] = -200; /* GC */
   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][g][g] = -160; /* GG */
   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][g][u] = -110; /* GU */

   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][u][a] = -170; /* UA */
   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][u][c] = -140; /* UC */
   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][u][g] = -180; /* UG */
   this->G_mismatch_hairpin[(int)this->bp_idx[c][g]][u][u] = -200; /* UU */

   /* GC */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][a][a] = -110; /* AA */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][a][c] = -150; /* AC */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][a][g] = -130; /* AG */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][a][u] = -210; /* AU */

   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][c][a] = -110; /* CA */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][c][c] =  -70; /* CC */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][c][g] = -240; /* CG */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][c][u] =  -50; /* CU */

   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][g][a] = -240; /* GA */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][g][c] = -290; /* GC */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][g][g] = -140; /* GG */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][g][u] = -120; /* GU */

   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][u][a] = -190; /* UA */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][u][c] = -100; /* UC */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][u][g] = -220; /* UG */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][c]][u][u] = -150; /* UU */

   /* GU */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][a][a] =  20; /* AA */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][a][c] = -50; /* AC */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][a][g] = -30; /* AG */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][a][u] = -30; /* AU */

   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][c][a] =  -10; /* CA */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][c][c] =  -20; /* CC */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][c][g] = -150; /* CG */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][c][u] =  -20; /* CU */

   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][g][a] =  -90; /* GA */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][g][c] = -110; /* GC */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][g][g] =  -30; /* GG */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][g][u] =    0; /* GU */

   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][u][a] =  -30; /* UA */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][u][c] =  -30; /* UC */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][u][g] =  -40; /* UG */
   this->G_mismatch_hairpin[(int)this->bp_idx[g][u]][u][u] = -110; /* UU */

   /* UG */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][a][a] = -50; /* AA */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][a][c] = -30; /* AC */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][a][g] = -60; /* AG */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][a][u] = -50; /* AU */

   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][c][a] =  -20; /* CA */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][c][c] =  -10; /* CC */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][c][g] = -170; /* CG */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][c][u] =    0; /* CU */

   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][g][a] =  -80; /* GA */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][g][c] = -120; /* GC */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][g][g] =  -30; /* GG */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][g][u] =  -70; /* GU */

   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][u][a] = -60; /* UA */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][u][c] = -10; /* UC */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][u][g] = -60; /* UG */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][g]][u][u] = -80; /* UU */

   /* AU */
   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][a][a] = -30; /* AA */
   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][a][c] = -50; /* AC */
   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][a][g] = -30; /* AG */
   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][a][u] = -30; /* AU */

   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][c][a] =  -10; /* CA */
   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][c][c] =  -20; /* CC */
   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][c][g] = -150; /* CG */
   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][c][u] =  -20; /* CU */

   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][g][a] = -110; /* GA */
   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][g][c] = -120; /* GC */
   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][g][g] =  -20; /* GG */
   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][g][u] =   20; /* GU */

   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][u][a] =  -30; /* UA */
   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][u][c] =  -30; /* UC */
   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][u][g] =  -60; /* UG */
   this->G_mismatch_hairpin[(int)this->bp_idx[a][u]][u][u] = -110; /* UU */

   /* UA */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][a][a] = -50; /* AA */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][a][c] = -30; /* AC */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][a][g] = -60; /* AG */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][a][u] = -50; /* AU */

   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][c][a] =  -20; /* CA */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][c][c] =  -10; /* CC */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][c][g] = -120; /* CG */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][c][u] =    0; /* CU */

   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][g][a] = -140; /* GA */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][g][c] = -120; /* GC */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][g][g] =  -70; /* GG */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][g][u] =  -20; /* GU */

   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][u][a] = -30; /* UA */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][u][c] = -10; /* UC */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][u][g] = -50; /* UG */
   this->G_mismatch_hairpin[(int)this->bp_idx[u][a]][u][u] = -80; /* UU */

   return 0;
}

static int
allocate_init_mismatch_interior (int a, int u, int g, int c,
                                 const unsigned long no_of_b,
                                 NN_scores* this,
                                 const char* file, const int line)
{
   unsigned long bp1;

   /* allocate memory */
   this->G_mismatch_interior
      = (int***) XOBJ_MALLOC_ND(sizeof (***this->G_mismatch_interior),
                                D_MM_I,
                                file, line,
                                this->bp_allowed_size, no_of_b, no_of_b);
   if (this->G_mismatch_interior == NULL)
   {
      return 1;        
   }
   this->G_mismatch_interior_size = this->bp_allowed_size * no_of_b * no_of_b;

   /* CG */
   bp1 = this->bp_idx[c][g];
   this->G_mismatch_interior[bp1][a][a] =    0; /* AA */
   this->G_mismatch_interior[bp1][a][c] =    0; /* AC */
   this->G_mismatch_interior[bp1][a][g] = -110; /* AG */
   this->G_mismatch_interior[bp1][a][u] =    0; /* AU */

   this->G_mismatch_interior[bp1][c][a] =    0; /* CA */
   this->G_mismatch_interior[bp1][c][c] =    0; /* CC */
   this->G_mismatch_interior[bp1][c][g] =    0; /* CG */
   this->G_mismatch_interior[bp1][c][u] =    0; /* CU */

   this->G_mismatch_interior[bp1][g][a] = -110; /* GA */
   this->G_mismatch_interior[bp1][g][c] =    0; /* GC */
   this->G_mismatch_interior[bp1][g][g] =    0; /* GG */
   this->G_mismatch_interior[bp1][g][u] =    0; /* GU */

   this->G_mismatch_interior[bp1][u][a] =    0; /* UA */
   this->G_mismatch_interior[bp1][u][c] =    0; /* UC */
   this->G_mismatch_interior[bp1][u][g] =    0; /* UG */
   this->G_mismatch_interior[bp1][u][u] =  -70; /* UU */

   /* GC */
   bp1 = this->bp_idx[g][c];
   this->G_mismatch_interior[bp1][a][a] =    0; /* AA */
   this->G_mismatch_interior[bp1][a][c] =    0; /* AC */
   this->G_mismatch_interior[bp1][a][g] = -110; /* AG */
   this->G_mismatch_interior[bp1][a][u] =    0; /* AU */

   this->G_mismatch_interior[bp1][c][a] =    0; /* CA */
   this->G_mismatch_interior[bp1][c][c] =    0; /* CC */
   this->G_mismatch_interior[bp1][c][g] =    0; /* CG */
   this->G_mismatch_interior[bp1][c][u] =    0; /* CU */

   this->G_mismatch_interior[bp1][g][a] = -110; /* GA */
   this->G_mismatch_interior[bp1][g][c] =    0; /* GC */
   this->G_mismatch_interior[bp1][g][g] =    0; /* GG */
   this->G_mismatch_interior[bp1][g][u] =    0; /* GU */

   this->G_mismatch_interior[bp1][u][a] =    0; /* UA */
   this->G_mismatch_interior[bp1][u][c] =    0; /* UC */
   this->G_mismatch_interior[bp1][u][g] =    0; /* UG */
   this->G_mismatch_interior[bp1][u][u] =  -70; /* UU */

   /* GU */
   bp1 = this->bp_idx[g][u];
   this->G_mismatch_interior[bp1][a][a] =   70; /* AA */
   this->G_mismatch_interior[bp1][a][c] =   70; /* AC */
   this->G_mismatch_interior[bp1][a][g] =  -40; /* AG */
   this->G_mismatch_interior[bp1][a][u] =   70; /* AU */

   this->G_mismatch_interior[bp1][c][a] =   70; /* CA */
   this->G_mismatch_interior[bp1][c][c] =   70; /* CC */
   this->G_mismatch_interior[bp1][c][g] =   70; /* CG */
   this->G_mismatch_interior[bp1][c][u] =   70; /* CU */

   this->G_mismatch_interior[bp1][g][a] =  -40; /* GA */
   this->G_mismatch_interior[bp1][g][c] =   70; /* GC */
   this->G_mismatch_interior[bp1][g][g] =   70; /* GG */
   this->G_mismatch_interior[bp1][g][u] =   70; /* GU */

   this->G_mismatch_interior[bp1][u][a] =   70; /* UA */
   this->G_mismatch_interior[bp1][u][c] =   70; /* UC */
   this->G_mismatch_interior[bp1][u][g] =   70; /* UG */
   this->G_mismatch_interior[bp1][u][u] =    0; /* UU */

   /* UG */
   bp1 = this->bp_idx[u][g];
   this->G_mismatch_interior[bp1][a][a] =   70; /* AA */
   this->G_mismatch_interior[bp1][a][c] =   70; /* AC */
   this->G_mismatch_interior[bp1][a][g] =  -40; /* AG */
   this->G_mismatch_interior[bp1][a][u] =   70; /* AU */

   this->G_mismatch_interior[bp1][c][a] =   70; /* CA */
   this->G_mismatch_interior[bp1][c][c] =   70; /* CC */
   this->G_mismatch_interior[bp1][c][g] =   70; /* CG */
   this->G_mismatch_interior[bp1][c][u] =   70; /* CU */

   this->G_mismatch_interior[bp1][g][a] =  -40; /* GA */
   this->G_mismatch_interior[bp1][g][c] =   70; /* GC */
   this->G_mismatch_interior[bp1][g][g] =   70; /* GG */
   this->G_mismatch_interior[bp1][g][u] =   70; /* GU */

   this->G_mismatch_interior[bp1][u][a] =   70; /* UA */
   this->G_mismatch_interior[bp1][u][c] =   70; /* UC */
   this->G_mismatch_interior[bp1][u][g] =   70; /* UG */
   this->G_mismatch_interior[bp1][u][u] =    0; /* UU */

   /* AU */
   bp1 = this->bp_idx[a][u];
   this->G_mismatch_interior[bp1][a][a] =   70; /* AA */
   this->G_mismatch_interior[bp1][a][c] =   70; /* AC */
   this->G_mismatch_interior[bp1][a][g] =  -40; /* AG */
   this->G_mismatch_interior[bp1][a][u] =   70; /* AU */

   this->G_mismatch_interior[bp1][c][a] =   70; /* CA */
   this->G_mismatch_interior[bp1][c][c] =   70; /* CC */
   this->G_mismatch_interior[bp1][c][g] =   70; /* CG */
   this->G_mismatch_interior[bp1][c][u] =   70; /* CU */

   this->G_mismatch_interior[bp1][g][a] =  -40; /* GA */
   this->G_mismatch_interior[bp1][g][c] =   70; /* GC */
   this->G_mismatch_interior[bp1][g][g] =   70; /* GG */
   this->G_mismatch_interior[bp1][g][u] =   70; /* GU */

   this->G_mismatch_interior[bp1][u][a] =   70; /* UA */
   this->G_mismatch_interior[bp1][u][c] =   70; /* UC */
   this->G_mismatch_interior[bp1][u][g] =   70; /* UG */
   this->G_mismatch_interior[bp1][u][u] =    0; /* UU */

   /* UA */
   bp1 = this->bp_idx[u][a];
   this->G_mismatch_interior[bp1][a][a] =   70; /* AA */
   this->G_mismatch_interior[bp1][a][c] =   70; /* AC */
   this->G_mismatch_interior[bp1][a][g] =  -40; /* AG */
   this->G_mismatch_interior[bp1][a][u] =   70; /* AU */

   this->G_mismatch_interior[bp1][c][a] =   70; /* CA */
   this->G_mismatch_interior[bp1][c][c] =   70; /* CC */
   this->G_mismatch_interior[bp1][c][g] =   70; /* CG */
   this->G_mismatch_interior[bp1][c][u] =   70; /* CU */

   this->G_mismatch_interior[bp1][g][a] =  -40; /* GA */
   this->G_mismatch_interior[bp1][g][c] =   70; /* GC */
   this->G_mismatch_interior[bp1][g][g] =   70; /* GG */
   this->G_mismatch_interior[bp1][g][u] =   70; /* GU */

   this->G_mismatch_interior[bp1][u][a] =   70; /* UA */
   this->G_mismatch_interior[bp1][u][c] =   70; /* UC */
   this->G_mismatch_interior[bp1][u][g] =   70; /* UG */
   this->G_mismatch_interior[bp1][u][u] =    0; /* UU */

   return 0;
}

static int
allocate_init_internal_loop (NN_scores* this, const char* file, const int line)
{

   this->G_internal_loop_size = 31;
   this->G_internal_loop = XOBJ_MALLOC (  this->G_internal_loop_size
                                          * sizeof (this->G_internal_loop[0]),
                                          file, line);
   if (this->G_internal_loop == NULL)
   {
      return 1;
   }

   this->G_internal_loop[ 0] = INT_UNDEF; /* min. loop length is 2 */
   this->G_internal_loop[ 1] = INT_UNDEF;
   this->G_internal_loop[ 2] = 410;
   this->G_internal_loop[ 3] = 510;
   this->G_internal_loop[ 4] = 170;
   this->G_internal_loop[ 5] = 180;
   this->G_internal_loop[ 6] = 200;
   this->G_internal_loop[ 7] = 220;
   this->G_internal_loop[ 8] = 230;
   this->G_internal_loop[ 9] = 240;
   this->G_internal_loop[10] = 250;
   this->G_internal_loop[11] = 260;
   this->G_internal_loop[12] = 270;
   this->G_internal_loop[13] = 278;
   this->G_internal_loop[14] = 286;
   this->G_internal_loop[15] = 294;
   this->G_internal_loop[16] = 301;
   this->G_internal_loop[17] = 307;
   this->G_internal_loop[18] = 313;
   this->G_internal_loop[19] = 319;
   this->G_internal_loop[20] = 325;
   this->G_internal_loop[21] = 330;
   this->G_internal_loop[22] = 335;
   this->G_internal_loop[23] = 340;
   this->G_internal_loop[24] = 345;
   this->G_internal_loop[25] = 349;
   this->G_internal_loop[26] = 353;
   this->G_internal_loop[27] = 357;
   this->G_internal_loop[28] = 361;
   this->G_internal_loop[29] = 365;
   this->G_internal_loop[30] = 369;

   return 0;
}

static int
allocate_init_int11 (const int a, const int u, const int g, const int c,
                     const unsigned long no_of_b,
                     NN_scores* this,
                     const char* file, const int line)
{
   unsigned long bp1, bp2;

   /* allocate memory */
   this->G_int11 = (int****) XOBJ_MALLOC_ND(sizeof (****this->G_int11),
                                            D_INT11,
                                            file, line,
                                            this->bp_allowed_size,
                                            this->bp_allowed_size,
                                            no_of_b,
                                            no_of_b);
   if (this->G_int11 == NULL)
   {
      return 1;        
   }
   this->G_int11_size = this->bp_allowed_size
                      * this->bp_allowed_size
                      * no_of_b
                      * no_of_b;

   /* CG */
   bp1 = this->bp_idx[c][g];
   /*    CG */
   bp2 = this->bp_idx[c][g];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =   40; /*   C */
   this->G_int11[bp1][bp2][a][g] =   40; /*   G */
   this->G_int11[bp1][bp2][a][u] =   40; /*   U */
   this->G_int11[bp1][bp2][c][a] =   40; /* C A */
   this->G_int11[bp1][bp2][c][c] =   40; /*   C */
   this->G_int11[bp1][bp2][c][g] =   40; /*   G */
   this->G_int11[bp1][bp2][c][u] =   40; /*   U */
   this->G_int11[bp1][bp2][g][a] =   40; /* G A */
   this->G_int11[bp1][bp2][g][c] =   40; /*   C */
   this->G_int11[bp1][bp2][g][g] = -140; /*   G */
   this->G_int11[bp1][bp2][g][u] =   40; /*   U */
   this->G_int11[bp1][bp2][u][a] =   40; /* U A */
   this->G_int11[bp1][bp2][u][c] =   40; /*   C */
   this->G_int11[bp1][bp2][u][g] =   40; /*   G */
   this->G_int11[bp1][bp2][u][u] =   40; /*   U */
   /*    GC */
   bp2 = this->bp_idx[g][c];
   this->G_int11[bp1][bp2][a][a] =   40; /* A A */
   this->G_int11[bp1][bp2][a][c] =  -40; /*   C */
   this->G_int11[bp1][bp2][a][g] =   40; /*   G */
   this->G_int11[bp1][bp2][a][u] =   40; /*   U */
   this->G_int11[bp1][bp2][c][a] =   30; /* C A */
   this->G_int11[bp1][bp2][c][c] =   50; /*   C */
   this->G_int11[bp1][bp2][c][g] =   40; /*   G */
   this->G_int11[bp1][bp2][c][u] =   50; /*   U */
   this->G_int11[bp1][bp2][g][a] =  -10; /* G A */
   this->G_int11[bp1][bp2][g][c] =   40; /*   C */
   this->G_int11[bp1][bp2][g][g] = -170; /*   G */
   this->G_int11[bp1][bp2][g][u] =   40; /*   U */
   this->G_int11[bp1][bp2][u][a] =   40; /* U A */
   this->G_int11[bp1][bp2][u][c] =    0; /*   C */
   this->G_int11[bp1][bp2][u][g] =   40; /*   G */
   this->G_int11[bp1][bp2][u][u] =  -30; /*   U */
   /*    GU */
   bp2 = this->bp_idx[g][u];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110; /*   U */
   /*    UG */
   bp2 = this->bp_idx[u][g];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 ; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110; /*   U */
   /*    AU */
   bp2 = this->bp_idx[a][u];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110; /*   U */
   /*    UA */
   bp2 = this->bp_idx[u][a];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110; /*   U */

   /* GC */
   bp1 = this->bp_idx[g][c];
   /*    CG */
   bp2 = this->bp_idx[c][g];
   this->G_int11[bp1][bp2][a][a] =   40; /* A A */
   this->G_int11[bp1][bp2][a][c] =   30; /*   C */
   this->G_int11[bp1][bp2][a][g] =  -10; /*   G */
   this->G_int11[bp1][bp2][a][u] =   40; /*   U */
   this->G_int11[bp1][bp2][c][a] =  -40; /* C A */
   this->G_int11[bp1][bp2][c][c] =   50; /*   C */
   this->G_int11[bp1][bp2][c][g] =   40; /*   G */
   this->G_int11[bp1][bp2][c][u] =    0; /*   U */
   this->G_int11[bp1][bp2][g][a] =   40; /* G A */
   this->G_int11[bp1][bp2][g][c] =   40; /*   C */
   this->G_int11[bp1][bp2][g][g] = -170; /*   G */
   this->G_int11[bp1][bp2][g][u] =   40; /*   U */
   this->G_int11[bp1][bp2][u][a] =   40; /* U A */
   this->G_int11[bp1][bp2][u][c] =   50; /*   C */
   this->G_int11[bp1][bp2][u][g] =   40; /*   G */
   this->G_int11[bp1][bp2][u][u] =  -30; /*   U */
   /*    GC */
   bp2 = this->bp_idx[g][c];
   this->G_int11[bp1][bp2][a][a] =   80; /* A A */
   this->G_int11[bp1][bp2][a][c] =   40; /*   C */
   this->G_int11[bp1][bp2][a][g] =   40; /*   G */
   this->G_int11[bp1][bp2][a][u] =   40; /*   U */
   this->G_int11[bp1][bp2][c][a] =   40; /* C A */
   this->G_int11[bp1][bp2][c][c] =   40; /*   C */
   this->G_int11[bp1][bp2][c][g] =   40; /*   G */
   this->G_int11[bp1][bp2][c][u] =   40; /*   U */
   this->G_int11[bp1][bp2][g][a] =   40; /* G A */
   this->G_int11[bp1][bp2][g][c] =   40; /*   C */
   this->G_int11[bp1][bp2][g][g] = -210; /*   G */
   this->G_int11[bp1][bp2][g][u] =   40; /*   U */
   this->G_int11[bp1][bp2][u][a] =   40; /* U A */
   this->G_int11[bp1][bp2][u][c] =   40; /*   C */
   this->G_int11[bp1][bp2][u][g] =   40; /*   G */
   this->G_int11[bp1][bp2][u][u] =  -70; /*   U */
   /*    GU */
   bp2 = this->bp_idx[g][u];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110; /*   U */
   /*    UG */
   bp2 = this->bp_idx[u][g];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110; /*   U */
   /*    AU */
   bp2 = this->bp_idx[a][u];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  100; /*   U */
   /*    UA */
   bp2 = this->bp_idx[u][a];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110; /*   U */

   /* GU */ 
   bp1 = this->bp_idx[g][u];
   /*    CG */
   bp2 = this->bp_idx[c][g];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110; /*   U */
   /*    GC */
   bp2 = this->bp_idx[g][c];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110; /*   U */
   /*    GU */
   bp2 = this->bp_idx[g][u];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170; /*   U */
   /*    UG */
   bp2 = this->bp_idx[u][g];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170; /*   U */
   /*    AU */
   bp2 = this->bp_idx[a][u];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170; /*   U */
   /*    UA */
   bp2 = this->bp_idx[u][a];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170; /*   U */
   /* UG */ 
   bp1 = this->bp_idx[u][g];
   /*    CG */
   bp2 = this->bp_idx[c][g];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110; /*   U */
   /*    GC */
   bp2 = this->bp_idx[g][c];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110; /*   U */
   /*    GU */
   bp2 = this->bp_idx[g][u];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170; /*   U */
   /*    UG */
   bp2 = this->bp_idx[u][g];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170; /*   U */
   /*    AU */
   bp2 = this->bp_idx[a][u];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170; /*   U */
   /*    UA */
   bp2 = this->bp_idx[u][a];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170; /*   U */
   /* AU */
   bp1 = this->bp_idx[a][u];
   /*    CG */
   bp2 = this->bp_idx[c][g];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110; /*   U */
   /*    GC */
   bp2 = this->bp_idx[g][c];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  100; /*   U */
   /*    GU */
   bp2 = this->bp_idx[g][u];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170; /*   U */
   /*    UG */
   bp2 = this->bp_idx[u][g];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170; /*   U */
   /*    AU */
   bp2 = this->bp_idx[a][u];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  120; /*   U */
   /*    UA */
   bp2 = this->bp_idx[u][a];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  150; /*   U */
   /* UA */
   bp1 = this->bp_idx[u][a];
   /*    CG */
   bp2 = this->bp_idx[c][g];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110; /*   U */
   /*    GC */
   bp2 = this->bp_idx[g][c];
   this->G_int11[bp1][bp2][a][a] =  110; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110; /*   U */
   /*    GU */
   bp2 = this->bp_idx[g][u];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170; /*   U */
   /*    UG */
   bp2 = this->bp_idx[u][g];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170; /*   U */
   /*    AU */
   bp2 = this->bp_idx[a][u];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  150; /*   U */
   /*    UA */
   bp2 = this->bp_idx[u][a];
   this->G_int11[bp1][bp2][a][a] =  170; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170; /*   G */
   this->G_int11[bp1][bp2][u][u] =  180; /*   U */

   return 0;
}

static int
allocate_init_non_gc_penalty_for_bp (int a, int u, int g, int c,
                                     NN_scores* this,
                                     const char* file, const int line)
{
   this->non_gc_penalty_for_bp = (int*) XOBJ_MALLOC (
                                        sizeof (int) * this->bp_allowed_size,
                                        file, line);
   if (this->non_gc_penalty_for_bp == NULL)
   {
      return 1;
   }

   this->non_gc_penalty_for_bp[(int)this->bp_idx[c][g]] = 0;
   this->non_gc_penalty_for_bp[(int)this->bp_idx[g][c]] = 0;
   this->non_gc_penalty_for_bp[(int)this->bp_idx[a][u]] = 50;
   this->non_gc_penalty_for_bp[(int)this->bp_idx[g][u]] = 50;
   this->non_gc_penalty_for_bp[(int)this->bp_idx[u][a]] = 50;
   this->non_gc_penalty_for_bp[(int)this->bp_idx[u][g]] = 50;

   return 0;
}

static int
allocate_init_bulge_loop (NN_scores* this, const char* file, const int line)
{

   this->G_bulge_loop_size = 31;
   this->G_bulge_loop = XOBJ_MALLOC (  this->G_bulge_loop_size
                                       * sizeof (this->G_bulge_loop[0]),
                                         file, line);
   if (this->G_bulge_loop == NULL)
   {
      return 1;        
   }

   this->G_bulge_loop[ 0] = INT_UNDEF; /* min. loop length is 1 */
   this->G_bulge_loop[ 1] = 380;
   this->G_bulge_loop[ 2] = 280;
   this->G_bulge_loop[ 3] = 320;
   this->G_bulge_loop[ 4] = 360;
   this->G_bulge_loop[ 5] = 400;
   this->G_bulge_loop[ 6] = 440;
   this->G_bulge_loop[ 7] = 459;
   this->G_bulge_loop[ 8] = 470;
   this->G_bulge_loop[ 9] = 480;
   this->G_bulge_loop[10] = 490;
   this->G_bulge_loop[11] = 500;
   this->G_bulge_loop[12] = 510;
   this->G_bulge_loop[13] = 519;
   this->G_bulge_loop[14] = 527;
   this->G_bulge_loop[15] = 534;
   this->G_bulge_loop[16] = 541;
   this->G_bulge_loop[17] = 548;
   this->G_bulge_loop[18] = 554;
   this->G_bulge_loop[19] = 560;
   this->G_bulge_loop[20] = 565;
   this->G_bulge_loop[21] = 571;
   this->G_bulge_loop[22] = 576;
   this->G_bulge_loop[23] = 580;
   this->G_bulge_loop[24] = 585;
   this->G_bulge_loop[25] = 589;
   this->G_bulge_loop[26] = 594;
   this->G_bulge_loop[27] = 598;
   this->G_bulge_loop[28] = 602;
   this->G_bulge_loop[29] = 605;
   this->G_bulge_loop[30] = 609;

   return 0;
}

static int
allocate_init_int21 (const int a, const int u, const int g, const int c,
                     const unsigned long no_of_b,
                     NN_scores* this,
                     const char* file, const int line)
{
   unsigned long bp1, bp2;

   /* allocate memory */
   this->G_int21 = (int*****) XOBJ_MALLOC_ND(sizeof (*****this->G_int21),
                                             D_INT21,
                                             file, line,
                                             this->bp_allowed_size,
                                             this->bp_allowed_size,
                                             no_of_b,
                                             no_of_b,
                                             no_of_b);
   if (this->G_int21 == NULL)
   {
      return 1;        
   }
   this->G_int21_size = this->bp_allowed_size
                      * this->bp_allowed_size
                      * no_of_b
                      * no_of_b
                      * no_of_b;

   /* CG */
   bp1 = this->bp_idx[c][g];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int21[bp1][bp2][a][a][a] = 240; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 220; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 400; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 210; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 170; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 160; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 400; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 100; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] =  60; /*     C */
   this->G_int21[bp1][bp2][a][g][g] =  40; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 400; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 400; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 230; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 220; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 220; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 220; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 220; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 400; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 250; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 190; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 220; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 170; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][a][g] =  80; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 400; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 400; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 400; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 400; /*     U */
   this->G_int21[bp1][bp2][g][g][a] =  80; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 220; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 400; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 400; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 400; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 400; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 400; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 220; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 130; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 400; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 400; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 170; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 120; /*     U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int21[bp1][bp2][a][a][a] = 230; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 220; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 110; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 400; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 210; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 170; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 160; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 400; /*     U */
   this->G_int21[bp1][bp2][a][g][a] =  80; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] =  60; /*     C */
   this->G_int21[bp1][bp2][a][g][g] =  40; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 400; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 400; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 230; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 220; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 220; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 220; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 220; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 400; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 250; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 190; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 220; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 170; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][a][g] =  80; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 400; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 400; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 400; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 400; /*     U */
   this->G_int21[bp1][bp2][g][g][a] =  80; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 220; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 400; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 400; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 400; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 400; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 400; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 220; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 150; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 400; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 400; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 170; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 120; /*     U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */

   /* GC */
   bp1 = this->bp_idx[g][c];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int21[bp1][bp2][a][a][a] = 250; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 220; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 210; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 400; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 210; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 170; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 160; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 400; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 120; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] =  60; /*     C */
   this->G_int21[bp1][bp2][a][g][g] =  40; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 400; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 400; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 230; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 220; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 220; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 220; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 220; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 400; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 250; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 190; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 220; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 170; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][a][g] =  80; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 400; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 400; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 400; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 400; /*     U */
   this->G_int21[bp1][bp2][g][g][a] =  80; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 220; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 400; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 400; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 400; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 400; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 400; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 220; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 120; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 400; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 400; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 170; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 120; /*     U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int21[bp1][bp2][a][a][a] = 240; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 220; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 400; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 210; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 170; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 160; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 400; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 100; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] =  60; /*     C */
   this->G_int21[bp1][bp2][a][g][g] =  40; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 400; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 400; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 230; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 220; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 220; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 220; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 220; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 400; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 250; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 190; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 220; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 170; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][a][g] =  80; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 400; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 400; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 400; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 400; /*     U */
   this->G_int21[bp1][bp2][g][g][a] =  80; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 220; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 400; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 400; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 400; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 400; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 400; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 400; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 220; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 130; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 400; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 400; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 400; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 170; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 400; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 120; /*     U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */

   /* GU */
   bp1 = this->bp_idx[g][u];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */

   /* UG */
   bp1 = this->bp_idx[u][g];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */

   /* AU */
   bp1 = this->bp_idx[a][u];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */

   /* UA */
   bp1 = this->bp_idx[u][a];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int21[bp1][bp2][a][a][a] = 320; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200; /*     U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int21[bp1][bp2][a][a][a] = 390; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270; /*     U */

   return 0;
}

static int
allocate_init_int22 (const int a, const int u, const int g, const int c,
                     const unsigned long no_of_b,
                     NN_scores* this,
                     const char* file, const int line)
{
   unsigned long bp1, bp2;

   /* allocate memory */
   this->G_int22 = (int******) XOBJ_MALLOC_ND(sizeof (******this->G_int22),
                                              D_INT22,
                                              file, line,
                                              this->bp_allowed_size,
                                              this->bp_allowed_size,
                                              no_of_b,
                                              no_of_b,
                                              no_of_b,
                                              no_of_b);
   if (this->G_int22 == NULL)
   {
      return 1;        
   }
   this->G_int22_size = this->bp_allowed_size
                      * this->bp_allowed_size
                      * no_of_b
                      * no_of_b
                      * no_of_b
                      * no_of_b;

   /* CG */
   bp1 = this->bp_idx[c][g];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  130; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  160; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  120; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   20; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   30; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   60; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  160; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   60; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  210; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  190; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  130; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   30; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   60; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  140; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  -40; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -110; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -60; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   40; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  140; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =   80; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  130; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  120; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  210; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  110; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  140; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  120; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   20; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  150; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  130; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  150; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  140; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  120; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  150; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  140; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   20; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  150; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  130; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -150; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  -20; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -40; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  120; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  150; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =    0; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  130; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  110; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   30; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   60; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   50; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   30; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   20; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =   90; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =    0; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -70; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   90; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   60; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  150; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  130; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  130; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -70; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =   90; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =    0; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =   90; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -60; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -260; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  130; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =   90; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -110; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  130; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  -40; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  120; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -150; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   30; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  130; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -60; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  130; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -110; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  150; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   60; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  140; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -40; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =   50; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  130; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -60; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  130; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -40; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -420; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  -50; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  140; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -40; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   50; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  110; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -260; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =   50; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =  -40; /*       U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int22[bp1][bp2][a][a][a][a] =   50; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =   60; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  110; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  -30; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   10; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] = -160; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  110; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] = -100; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  -60; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =   70; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =   50; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =   20; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   40; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   50; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  140; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   10; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -70; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -80; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  180; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  120; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  -50; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  -60; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  -60; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  150; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  130; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  100; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  130; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  120; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  -70; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =   70; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =   40; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  100; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  110; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  100; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  130; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  120; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =    0; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =   30; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   70; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =   70; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  100; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -190; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  -30; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -70; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  110; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  140; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  120; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] = -150; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  -20; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  -20; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  -10; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   20; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  -20; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =   90; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  -40; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =   90; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =    0; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] = -170; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   30; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   70; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  -10; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  110; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  100; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   20; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =   50; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =    0; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -50; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =   60; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =    0; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =   90; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -90; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] = -100; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -300; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  120; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  100; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] = -130; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  -60; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -240; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =   90; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =   60; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  120; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -160; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   30; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   40; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] = -160; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =   50; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  110; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -160; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =   30; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  120; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  -60; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   30; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   20; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] = -160; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =   10; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =   50; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -60; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -70; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -440; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] = -100; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  120; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   30; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  -10; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -410; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =   10; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =   40; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] = -100; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   60; /*       U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  200; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  180; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   80; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   80; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  190; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  230; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   80; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   60; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   40; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  160; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   90; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  190; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  130; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  180; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   70; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  180; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  210; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  160; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   70; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  120; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   80; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   80; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  110; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  100; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  100; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   80; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   60; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -20; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   90; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   70; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  170; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -20; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   30; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   40; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   40; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -160; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  140; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  130; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -60; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  180; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -80; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  100; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -40; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   70; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -310; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   60; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -210; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   10; /*       U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  200; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   60; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  130; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  260; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  260; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  100; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  130; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  140; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   40; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   20; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -40; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =    0; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  130; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  150; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  190; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  190; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   90; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -90; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   40; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   30; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   40; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  110; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  140; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  100; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   60; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =    0; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  100; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   80; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =    0; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   40; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =    0; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =    0; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -200; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -70; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  190; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -110; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   70; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   30; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   20; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  150; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -220; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  120; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30; /*       U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  200; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  180; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   80; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   80; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  190; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  230; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   80; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   60; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   40; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  160; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   90; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  190; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  130; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  180; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   70; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  180; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  210; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  160; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   70; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  120; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   80; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   80; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  110; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  100; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  100; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   80; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   60; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -20; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   90; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   70; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  170; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -20; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   30; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   40; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   40; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -160; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  140; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  130; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -60; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  180; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -80; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  100; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -40; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   70; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -310; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   60; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -210; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   10; /*       U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int22[bp1][bp2][a][a][a][a] =  200; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   60; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  130; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  260; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  260; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  100; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  130; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  140; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   40; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   20; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -40; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =    0; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  130; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  150; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  190; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  190; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   90; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -90; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   40; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   30; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   40; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  110; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  140; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  100; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   60; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =    0; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  100; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   80; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =    0; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   40; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =    0; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =    0; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -200; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -70; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  190; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -110; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   70; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   30; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   20; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  150; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -220; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  120; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30; /*       U */

/*yy*/
   /* GC */
   bp1 = this->bp_idx[g][c];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int22[bp1][bp2][a][a][a][a] =   50; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   40; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  130; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  100; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  -20; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   70; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =   60; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =    0; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] = -100; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  -10; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -160; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -60; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  -10; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =   60; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  140; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =   30; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  110; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  100; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  100; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  110; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  -40; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  110; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  120; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  150; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  150; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  130; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  130; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  120; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  120; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  120; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  -70; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  -60; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  120; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  150; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -160; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  -60; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -50; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  120; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  120; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  120; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =    0; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  100; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  100; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   30; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   30; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   30; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  -30; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  -50; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  -70; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =   90; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] = -150; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] = -170; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] = -130; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   10; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  -60; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =   70; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  -20; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   40; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  -10; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] = -160; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  -60; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  -90; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  -60; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] = -160; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -410; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =   40; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   30; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =   90; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -240; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =   50; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =   10; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =   70; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =    0; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -190; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  -20; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   20; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =   50; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =    0; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =   30; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  -10; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   20; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =   40; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =   20; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -80; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =   50; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] = -100; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] = -160; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -440; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] = -100; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  170; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   20; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =    0; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -300; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =   60; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =   10; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] = -100; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   60; /*       U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int22[bp1][bp2][a][a][a][a] =  150; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  120; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   10; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  120; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =   90; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  -50; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  -80; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] = -190; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  120; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =   90; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  180; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =   90; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =   80; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   10; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  -20; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] = -130; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  -70; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -200; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] = -130; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  100; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  -30; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] = -160; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  150; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =   20; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  120; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  100; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  100; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  100; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  -80; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =   20; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =   30; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =   90; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =   90; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  100; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  100; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  100; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  100; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =    0; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =    0; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  -10; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =   90; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =   90; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  150; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -190; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  -90; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -90; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  100; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  100; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  110; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] = -150; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  -50; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  -50; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   20; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   20; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   30; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  -50; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  -80; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =   90; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] = -150; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] = -260; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] = -150; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  -80; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] = -160; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =   20; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  -50; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  -80; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =   50; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] = -150; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] = -190; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  -90; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  -90; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  -60; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] = -190; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] = -100; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -450; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =   30; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  -50; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] = -150; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  -60; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -410; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =   30; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  -50; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =   80; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =    0; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -190; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   20; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  -80; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] = -190; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =   30; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =    0; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -200; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =   20; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =    0; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   20; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] = -100; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] = -190; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  -70; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  -10; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] = -130; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =   50; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] = -100; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] = -190; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -490; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  -90; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =    0; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   30; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] = -150; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -450; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  -50; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  -70; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =  -50; /*       U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  210; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   60; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   60; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  240; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   90; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   60; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  140; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   70; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -60; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   10; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  170; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   40; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =   70; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   60; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  170; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  160; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   60; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  120; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   50; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   50; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   40; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  140; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   80; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   80; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   80; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   10; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  -10; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  -90; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] = -110; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  -10; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =    0; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  -80; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =   80; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =   10; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   80; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   10; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] = -110; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  -10; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  -60; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -50; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   40; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -310; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =   80; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =   10; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  130; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -210; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =   80; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =   10; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  150; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  210; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =   60; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -130; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   70; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =   70; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -120; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  100; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =   60; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   60; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =   70; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -50; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   50; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =   60; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   50; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -250; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  -20; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =    0; /*       U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  210; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  140; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   30; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  210; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  270; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  110; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   80; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  120; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   10; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   30; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -100; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  160; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   30; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  100; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  230; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  140; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  140; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   80; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  190; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  190; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   80; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  120; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -90; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   10; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   10; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   30; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  130; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  140; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  100; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   10; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  -30; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] = -110; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -90; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   10; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   10; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  -60; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  110; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  110; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -90; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =   10; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  -50; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -90; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =    0; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -350; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =   80; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =   10; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   40; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -220; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  110; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  190; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =   30; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  150; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  210; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =   40; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -150; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   70; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   90; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =   90; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -100; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =   90; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   90; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =   90; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -90; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -390; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   10; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =   60; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   40; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -260; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30; /*       U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  210; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   60; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   60; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  240; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   90; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   60; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  140; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   70; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -60; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   10; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  170; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   40; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =   70; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   60; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  170; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  160; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   60; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  120; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   50; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   50; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   40; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  140; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   80; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   80; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   80; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   10; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  -10; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  -90; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] = -110; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  -10; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =    0; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  -80; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =   80; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =   10; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   80; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   10; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] = -110; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  -10; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  -60; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -50; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   40; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -310; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =   80; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =   10; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  130; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -210; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =   80; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =   10; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  150; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  210; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =   60; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -130; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   70; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =   70; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -120; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  100; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =   60; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   60; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =   70; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -50; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   50; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =   60; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   50; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -250; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  -20; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =    0; /*       U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int22[bp1][bp2][a][a][a][a] =  210; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  140; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   30; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  210; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  270; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  110; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   80; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  120; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   10; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   30; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -100; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  160; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   30; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  100; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  230; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  140; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  140; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   80; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  190; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  190; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   80; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  120; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -90; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   10; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   10; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   30; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  130; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  140; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  100; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   10; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  -30; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] = -110; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -90; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   10; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   10; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  -60; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  110; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  110; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -90; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =   10; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  -50; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -90; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =    0; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -350; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =   80; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =   10; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   40; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -220; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  110; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  190; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =   30; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  150; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  210; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =   40; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -150; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   70; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   90; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =   90; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -100; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =   90; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   90; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =   90; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -90; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -390; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   10; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =   60; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   40; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -260; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30; /*       U */

/*yy*/
   /* GU */
   bp1 = this->bp_idx[g][u];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  200; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   80; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   90; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  240; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  280; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  180; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  100; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   90; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  160; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   30; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -80; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  270; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   70; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  180; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  180; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  230; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   80; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  210; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  210; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  180; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   80; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -90; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =    0; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  180; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   60; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  150; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   90; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   90; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   90; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   80; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  130; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   70; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   50; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -20; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  110; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   90; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  180; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -20; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  130; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -10; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -40; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -210; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  180; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  140; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -60; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  180; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  230; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   60; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   40; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  130; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   80; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  100; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   40; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   40; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -310; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =    0; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  100; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  150; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -160; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =   90; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   60; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   10; /*       U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int22[bp1][bp2][a][a][a][a] =  210; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   60; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   10; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =    0; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] = -110; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  180; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   60; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  160; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  150; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =   70; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   70; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   60; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  160; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =    0; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -120; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =   40; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  -80; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  210; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  -10; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =   80; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  160; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  150; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =   60; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =   60; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =   60; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   50; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  140; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  140; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -130; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  -30; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  -90; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =   10; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =   10; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   90; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   80; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   80; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   90; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   60; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] = -110; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  -60; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   50; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   60; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   70; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   50; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -50; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  130; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -50; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -250; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  -10; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -210; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  170; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  150; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -60; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =   70; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   50; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   60; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  -20; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  140; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   40; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -30; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   50; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   10; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -310; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =   10; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =   80; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   50; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =    0; /*       U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  140; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  140; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   80; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   20; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  190; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  240; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  230; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =   10; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  140; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  140; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   30; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  140; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  190; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   30; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   80; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   90; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -110; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  190; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -10; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   90; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  230; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   20; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   90; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -210; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -110; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  140; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  110; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   60; /*       U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  160; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  250; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  250; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  160; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -20; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  180; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  250; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  250; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  250; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  150; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  190; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -30; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   70; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  170; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  110; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   90; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  100; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  210; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  150; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  130; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  250; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  210; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   50; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -150; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  250; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  220; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  230; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   70; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   70; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -120; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80; /*       U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  140; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  140; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   80; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   20; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  190; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  240; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  230; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =   10; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  140; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  140; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   30; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  140; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  190; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   30; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   80; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   90; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -110; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  190; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -10; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   90; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  230; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   20; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   90; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -210; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -110; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  140; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  110; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   60; /*       U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  160; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  250; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  250; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  160; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -20; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  180; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  250; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  250; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  250; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  150; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  190; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -30; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   70; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  170; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  110; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   90; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  100; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  210; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  150; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  130; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  250; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  210; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   50; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -150; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  250; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  220; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  230; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   70; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   70; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -120; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80; /*       U */

/*yy*/
   /* UG */
   bp1 = this->bp_idx[u][g];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  200; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  100; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  240; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  130; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  280; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  100; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  100; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   30; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -70; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   10; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  270; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  160; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  120; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  190; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  160; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  150; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   60; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  190; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  180; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  160; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  210; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   60; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -110; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   30; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =    0; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   40; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   70; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  120; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   90; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  100; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   90; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =    0; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   40; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  130; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  130; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   80; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  170; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =    0; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   40; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  130; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   70; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   10; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -220; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  160; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -70; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  190; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   20; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -40; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   40; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   30; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   30; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -200; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  120; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   20; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30; /*       U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int22[bp1][bp2][a][a][a][a] =  210; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   80; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   10; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   10; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  180; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   80; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  150; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =   90; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   70; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   70; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =    0; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -100; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  190; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =   40; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  -60; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =   10; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  210; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  190; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  170; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  140; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  -30; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  110; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =   80; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  140; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  140; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =   40; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =   90; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =   60; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   30; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  140; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -150; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  150; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] = -110; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =   40; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =   10; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   70; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  110; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   80; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  110; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   80; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   30; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -90; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   80; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   30; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  180; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  130; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   90; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -30; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   10; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  100; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   70; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -30; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -260; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  190; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   10; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -220; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  190; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  100; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  110; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -100; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  100; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  100; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   90; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =    0; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -390; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  -10; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   40; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -350; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =   30; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  110; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30; /*       U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  160; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  100; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   30; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  130; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  210; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  220; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  110; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  160; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -10; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  130; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  230; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  120; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  150; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  190; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  100; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   90; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  160; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  240; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   90; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  110; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -120; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  240; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  250; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  130; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  130; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  230; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -150; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  230; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  170; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80; /*       U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  270; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  230; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  230; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  170; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   80; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  150; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   70; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  170; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  270; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   70; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   70; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -160; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  230; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  280; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  280; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  230; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -290; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -160; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   90; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =  110; /*       U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  160; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  100; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   30; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  130; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  210; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  220; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  110; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  160; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -10; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  130; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  230; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  120; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  150; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  190; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  100; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   90; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  160; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  240; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   90; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  110; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -120; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  240; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  250; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  130; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  130; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  230; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -150; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  230; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  170; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80; /*       U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  270; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  230; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  230; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  170; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   80; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  150; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   70; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  170; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  270; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   70; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   70; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -160; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  230; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  280; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  280; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  230; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -290; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -160; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   90; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =  110; /*       U */

/*yy*/
   /* AU */
   bp1 = this->bp_idx[a][u];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  200; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   80; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   90; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  240; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  280; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  180; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  100; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   90; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  160; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   30; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -80; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  270; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   70; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  180; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  180; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  230; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   80; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  210; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  210; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  180; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   80; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -90; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =    0; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  180; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   60; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  150; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   90; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   90; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   90; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   80; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  130; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   70; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   50; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -20; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  110; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   90; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  180; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -20; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  130; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -10; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -40; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -210; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  180; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  140; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -60; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  180; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  230; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   60; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   40; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  130; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   80; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  100; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   40; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   40; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -310; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =    0; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  100; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  150; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -160; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =   90; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   60; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   10; /*       U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int22[bp1][bp2][a][a][a][a] =  210; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   60; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   10; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =    0; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] = -110; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  180; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   60; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  160; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  150; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =   70; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   70; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   60; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  160; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =    0; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -120; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =   40; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  -80; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  210; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  -10; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =   80; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  160; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  150; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =   60; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =   60; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =   60; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   50; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  140; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  140; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -130; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  -30; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  -90; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =   10; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =   10; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   90; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   80; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   80; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   90; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   60; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] = -110; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  -60; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   50; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   60; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   70; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   50; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -50; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  130; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -50; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -250; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  170; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  -10; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -210; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  170; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  150; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -60; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =   70; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   50; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   60; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  -20; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  140; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   40; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -30; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   50; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   10; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -310; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =   10; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =   80; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   50; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =    0; /*       U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  140; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  140; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   80; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   20; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  190; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  240; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  230; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =   10; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  140; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  140; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   30; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  140; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  190; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   30; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   80; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   90; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -110; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  190; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -10; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   90; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  230; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   20; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   90; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -210; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -110; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  140; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  110; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   60; /*       U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  160; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  250; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  250; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  160; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -20; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  180; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  250; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  250; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  250; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  150; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  190; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -30; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   70; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  170; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  110; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   90; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  100; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  210; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  150; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  130; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  250; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  210; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   50; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -150; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  250; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  220; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  230; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   70; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   70; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -120; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80; /*       U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  140; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  140; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  140; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   80; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   20; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  190; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  240; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  230; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =   10; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  140; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  140; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   30; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  140; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  190; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   30; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   80; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   90; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -110; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  190; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -10; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   90; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  230; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   20; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   90; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -210; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -110; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  140; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  110; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   60; /*       U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  160; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  250; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  250; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  160; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -20; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  180; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  250; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  250; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  250; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  150; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  190; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -30; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   70; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  170; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  160; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  110; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   90; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  100; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  210; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  150; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  130; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  250; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  230; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  210; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   50; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -150; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  250; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  220; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  230; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -10; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   70; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   70; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -120; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  160; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  130; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80; /*       U */

/*yy*/
   /* UA */
   bp1 = this->bp_idx[u][a];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  200; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  100; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  240; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  130; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  280; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  100; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  100; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   30; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -70; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   10; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  270; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   80; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  160; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  220; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  120; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  190; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  160; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  150; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   60; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  190; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  180; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  160; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  210; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   60; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -110; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   30; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =    0; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   40; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   70; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  120; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   90; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  100; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   90; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =    0; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   40; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  130; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  130; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   80; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  220; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  220; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  170; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =    0; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   40; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  130; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   70; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   10; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -220; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  160; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -70; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  190; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   20; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -40; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   40; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   30; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   30; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  210; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  120; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -200; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  120; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   20; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30; /*       U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int22[bp1][bp2][a][a][a][a] =  210; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   80; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   10; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   10; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  180; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   80; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  150; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =   90; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =   90; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   70; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   70; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  180; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =    0; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -100; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  190; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =   40; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  -60; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =   10; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  210; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  110; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  190; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  170; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  140; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  -30; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  110; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =   80; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  140; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  150; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  140; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =   40; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =   90; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =   60; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   30; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  140; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -150; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  150; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] = -110; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =   40; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =   10; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   70; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  110; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   80; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  110; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   80; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   30; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -90; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   80; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  120; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   30; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  180; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  180; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  130; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   90; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =   80; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   40; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -30; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   10; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  100; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   70; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -30; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -70; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -260; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  190; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   10; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -220; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  190; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  140; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  100; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  110; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  180; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -100; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  100; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  100; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   90; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -90; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =    0; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -390; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  -10; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  110; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   40; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -350; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =   30; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  110; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   10; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30; /*       U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  160; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  100; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   30; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  130; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  210; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  220; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  110; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  160; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -10; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  130; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  230; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  120; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  150; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  190; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  100; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   90; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  160; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  240; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   90; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  110; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -120; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  240; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  250; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  130; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  130; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  230; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -150; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  230; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  170; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80; /*       U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  270; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  230; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  230; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  170; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   80; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  150; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   70; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  170; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  270; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   70; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   70; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -160; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  230; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  280; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  280; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  230; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -290; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -160; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   90; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =  110; /*       U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  150; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  160; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  150; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  100; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   30; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  130; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  210; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  220; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  110; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  160; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -10; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  130; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   90; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  230; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  120; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  170; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  150; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  190; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  100; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   90; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  160; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  110; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  240; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   90; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  100; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  150; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  110; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -120; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  240; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -30; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  170; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   50; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  250; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  260; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -20; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  180; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   50; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  130; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  130; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  230; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -150; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  230; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  170; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80; /*       U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int22[bp1][bp2][a][a][a][a] =  280; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  170; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  270; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  280; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  170; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =    0; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   70; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  190; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  230; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  230; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  230; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  280; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  170; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   80; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  150; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  170; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   80; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   70; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  170; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  270; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   70; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  110; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  110; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  150; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   70; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -160; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -30; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  270; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  230; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  230; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  150; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  280; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =    0; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  280; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  280; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  190; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  270; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  230; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -290; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   90; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  220; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -160; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  220; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  190; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   90; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =  110; /*       U */

/*yy*/


   return 0;
}

static int
allocate_init_dangle5 (const int a, const int u, const int g, const int c, 
                       const unsigned long size,
                       NN_scores* this,
                       const char* file, const int line)
{
   this->G_dangle5 = (int**) XOBJ_MALLOC_2D (this->bp_allowed_size, size,
                                             sizeof (int), file, line);

   this->G_dangle5_size = this->bp_allowed_size * size;

   if (this->G_dangle5 == NULL)
   {
      return 1;
   }   

   /* CG */
   this->G_dangle5[(int)this->bp_idx[c][g]][a] = -50;
   this->G_dangle5[(int)this->bp_idx[c][g]][c] = -30;
   this->G_dangle5[(int)this->bp_idx[c][g]][g] = -20;
   this->G_dangle5[(int)this->bp_idx[c][g]][u] = -10;

   /* GC */
   this->G_dangle5[(int)this->bp_idx[g][c]][a] = -20;
   this->G_dangle5[(int)this->bp_idx[g][c]][c] = -30;
   this->G_dangle5[(int)this->bp_idx[g][c]][g] =   0;
   this->G_dangle5[(int)this->bp_idx[g][c]][u] =   0;

   /* GU */ 
   this->G_dangle5[(int)this->bp_idx[g][u]][a] = -30;
   this->G_dangle5[(int)this->bp_idx[g][u]][c] = -30;
   this->G_dangle5[(int)this->bp_idx[g][u]][g] = -40;
   this->G_dangle5[(int)this->bp_idx[g][u]][u] = -20;

   /* UG */
   this->G_dangle5[(int)this->bp_idx[u][g]][a] = -30;
   this->G_dangle5[(int)this->bp_idx[u][g]][c] = -10;
   this->G_dangle5[(int)this->bp_idx[u][g]][g] = -20;
   this->G_dangle5[(int)this->bp_idx[u][g]][u] = -20;

   /* AU */
   this->G_dangle5[(int)this->bp_idx[a][u]][a] = -30;
   this->G_dangle5[(int)this->bp_idx[a][u]][c] = -30;
   this->G_dangle5[(int)this->bp_idx[a][u]][g] = -40;
   this->G_dangle5[(int)this->bp_idx[a][u]][u] = -20;

   /* UA */
   this->G_dangle5[(int)this->bp_idx[u][a]][a] = -30;
   this->G_dangle5[(int)this->bp_idx[u][a]][c] = -10;
   this->G_dangle5[(int)this->bp_idx[u][a]][g] = -20;
   this->G_dangle5[(int)this->bp_idx[u][a]][u] = -20;

   return 0;
}

static int
allocate_init_dangle3 (const int a, const int u, const int g, const int c, 
                       const unsigned long size,
                       NN_scores* this,
                       const char* file, const int line)
{
   this->G_dangle3 = (int**) XOBJ_MALLOC_2D (this->bp_allowed_size, size,
                                             sizeof (int), file, line);

   this->G_dangle3_size = this->bp_allowed_size * size;

   if (this->G_dangle3 == NULL)
   {
      return 1;
   }   

   /* CG */
   this->G_dangle3[(int)this->bp_idx[c][g]][a] = -110;
   this->G_dangle3[(int)this->bp_idx[c][g]][c] =  -40;
   this->G_dangle3[(int)this->bp_idx[c][g]][g] = -130;
   this->G_dangle3[(int)this->bp_idx[c][g]][u] =  -60;

   /* GC */
   this->G_dangle3[(int)this->bp_idx[g][c]][a] = -170;
   this->G_dangle3[(int)this->bp_idx[g][c]][c] =  -80;
   this->G_dangle3[(int)this->bp_idx[g][c]][g] = -170;
   this->G_dangle3[(int)this->bp_idx[g][c]][u] = -120;

   /* GU */ 
   this->G_dangle3[(int)this->bp_idx[g][u]][a] = -70;
   this->G_dangle3[(int)this->bp_idx[g][u]][c] = -10;
   this->G_dangle3[(int)this->bp_idx[g][u]][g] = -70;
   this->G_dangle3[(int)this->bp_idx[g][u]][u] = -10;

   /* UG */
   this->G_dangle3[(int)this->bp_idx[u][g]][a] = -80;
   this->G_dangle3[(int)this->bp_idx[u][g]][c] = -50;
   this->G_dangle3[(int)this->bp_idx[u][g]][g] = -80;
   this->G_dangle3[(int)this->bp_idx[u][g]][u] = -60;

   /* AU */
   this->G_dangle3[(int)this->bp_idx[a][u]][a] = -70;
   this->G_dangle3[(int)this->bp_idx[a][u]][c] = -10;
   this->G_dangle3[(int)this->bp_idx[a][u]][g] = -70;
   this->G_dangle3[(int)this->bp_idx[a][u]][u] = -10;

   /* UA */
   this->G_dangle3[(int)this->bp_idx[u][a]][a] = -80;
   this->G_dangle3[(int)this->bp_idx[u][a]][c] = -50;
   this->G_dangle3[(int)this->bp_idx[u][a]][g] = -80;
   this->G_dangle3[(int)this->bp_idx[u][a]][u] = -60;

   return 0;
}

void
tetra_loop_swap_entries (unsigned long src, unsigned long dest,
                         NN_scores* this)
{
   unsigned long i;
   char tmp[D_TL];
   int G_tmp;
   char* ctmp;

   /* copy scores first */
   G_tmp = this->G_tetra_loop[dest];
   this->G_tetra_loop[dest] = this->G_tetra_loop[src];
   this->G_tetra_loop[src] = G_tmp;

   if ((src == 0)||(dest == 0))
   {
      for (i = 0; i < D_TL; i++)
      {
         tmp[i] = this->tetra_loop[dest][i];
      }   
      for (i = 0; i < D_TL; i++)
      {
         this->tetra_loop[dest][i] = this->tetra_loop[src][i];
      }
      for (i = 0; i < D_TL; i++)
      {
         this->tetra_loop[src][i] = tmp[i];
      }
   }
   else
   {
      /* temp. cpy dest */
      ctmp = this->tetra_loop[dest];
      
      /* cpy src to dest */
      this->tetra_loop[dest] = this->tetra_loop[src];

      /* restore dest in src */
      this->tetra_loop[src] = ctmp;
   }

}

static int
tetra_loop_cmp (unsigned long idx1, unsigned long idx2, NN_scores* this)
{
   unsigned long i;

   for (i = 0; i < D_TL; i++)
   {
      if (this->tetra_loop[idx1][i] != this->tetra_loop[idx2][i])
      {
         return this->tetra_loop[idx1][i] - this->tetra_loop[idx2][i];
      }
   }

   return 0;
}

static void
tetra_loop_qsort (const unsigned long left, const unsigned long right,
                  NN_scores* this)
{
   unsigned long pivot = left;
   unsigned long i;
   unsigned long r;

   if (left >= right)
   {
      return;
   }

   r = right - 1;

   /* mfprintf (stderr, "left: %lu right: %lu ", left, r); */

   /* sort according to pivot */
   i = left + 1;
   while (i <= r)   
   {
      if (tetra_loop_cmp (pivot, i, this) < 0)
      {
         tetra_loop_swap_entries (i, r, this);
         r--;
      }
      else
      {
         i++;
      }
   }

   tetra_loop_swap_entries (pivot, r, this);   

   /* mfprintf (stderr, "middle: %lu\n", r); */

   tetra_loop_qsort (left, r, this);
   tetra_loop_qsort (r+1, right, this);
}

static int
allocate_init_tetra_loop (const int a, const int u, const int g, const int c,
                          NN_scores* this,
                          const char* file, const int line)
{
   char* l;

   this->tetra_loop_size = 30;
   this->tetra_loop = (char**) XOBJ_MALLOC_2D (this->tetra_loop_size, D_TL,
                                               sizeof (char),
                                               file, line);
   if (this->tetra_loop == NULL)
   {
      return 1;
   }

   this->G_tetra_loop = (int*) XOBJ_MALLOC (
                                           sizeof (int) * this->tetra_loop_size,
                                           file, line);

   /* GGGGAC -300 */
   l = this->tetra_loop[0];
   l[0] = g; l[1] = g; l[2] = g; l[3] = g; l[4] = a; l[5] = c;
   this->G_tetra_loop[0] = -300;

   /* GGUGAC -300 */
   l = this->tetra_loop[1];
   l[0] = g; l[1] = g; l[2] = u; l[3] = g; l[4] = a; l[5] = c;
   this->G_tetra_loop[1] = -300;

   /* CGAAAG -300 */
   l = this->tetra_loop[2];
   l[0] = c; l[1] = g; l[2] = a; l[3] = a; l[4] = a; l[5] = g;
   this->G_tetra_loop[2] = -300;

   /* GGAGAC -300 */
   l = this->tetra_loop[3];
   l[0] = g; l[1] = g; l[2] = a; l[3] = g; l[4] = a; l[5] = c;
   this->G_tetra_loop[3] = -300;

   /* CGCAAG -300 */
   l = this->tetra_loop[4];
   l[0] = c; l[1] = g; l[2] = c; l[3] = a; l[4] = a; l[5] = g;
   this->G_tetra_loop[4] = -300;

   /* GGAAAC -300 */
   l = this->tetra_loop[5];
   l[0] = g; l[1] = g; l[2] = a; l[3] = a; l[4] = a; l[5] = c;
   this->G_tetra_loop[5] = -300;

   /* CGGAAG -300 */
   l = this->tetra_loop[6];
   l[0] = c; l[1] = g; l[2] = g; l[3] = a; l[4] = a; l[5] = g;
   this->G_tetra_loop[6] = -300;

   /* CUUCGG -300 */
   l = this->tetra_loop[7];
   l[0] = c; l[1] = u; l[2] = u; l[3] = c; l[4] = g; l[5] = g;
   this->G_tetra_loop[7] = -300;

   /* CGUGAG -300 */
   l = this->tetra_loop[8];
   l[0] = c; l[1] = g; l[2] = u; l[3] = g; l[4] = a; l[5] = g;
   this->G_tetra_loop[8] = -300;

   /* CGAAGG -250 */
   l = this->tetra_loop[9];
   l[0] = c; l[1] = g; l[2] = a; l[3] = a; l[4] = g; l[5] = g;
   this->G_tetra_loop[9] = -250;

   /* CUACGG -250 */
   l = this->tetra_loop[10];
   l[0] = c; l[1] = u; l[2] = a; l[3] = c; l[4] = g; l[5] = g;
   this->G_tetra_loop[10] = -250;

   /* GGCAAC -250 */
   l = this->tetra_loop[11];
   l[0] = g; l[1] = g; l[2] = c; l[3] = a; l[4] = a; l[5] = c;
   this->G_tetra_loop[11] = -250;

   /* CGCGAG -250 */
   l = this->tetra_loop[12];
   l[0] = c; l[1] = g; l[2] = c; l[3] = g; l[4] = a; l[5] = g;
   this->G_tetra_loop[12] = -250;

   /* UGAGAG -250 */
   l = this->tetra_loop[13];
   l[0] = u; l[1] = g; l[2] = a; l[3] = g; l[4] = a; l[5] = g;
   this->G_tetra_loop[13] = -250;

   /* CGAGAG -200 */
   l = this->tetra_loop[14];
   l[0] = c; l[1] = g; l[2] = a; l[3] = g; l[4] = a; l[5] = g;
   this->G_tetra_loop[14] = -200;

   /* AGAAAU -200 */
   l = this->tetra_loop[15];
   l[0] = a; l[1] = g; l[2] = a; l[3] = a; l[4] = a; l[5] = u;
   this->G_tetra_loop[15] = -200;

   /* CGUAAG -200 */
   l = this->tetra_loop[16];
   l[0] = c; l[1] = g; l[2] = u; l[3] = a; l[4] = a; l[5] = g;
   this->G_tetra_loop[16] = -200;

   /* CUAACG -200 */
   l = this->tetra_loop[17];
   l[0] = c; l[1] = u; l[2] = a; l[3] = a; l[4] = c; l[5] = g;
   this->G_tetra_loop[17] = -200;

   /* UGAAAG -200 */
   l = this->tetra_loop[18];
   l[0] = u; l[1] = g; l[2] = a; l[3] = a; l[4] = a; l[5] = g;
   this->G_tetra_loop[18] = -200;
 
   /* GGAAGC -150 */
   l = this->tetra_loop[19];
   l[0] = g; l[1] = g; l[2] = a; l[3] = a; l[4] = g; l[5] = c;
   this->G_tetra_loop[19] = -150;

   /* GGGAAC -150 */
   l = this->tetra_loop[20];
   l[0] = g; l[1] = g; l[2] = g; l[3] = a; l[4] = a; l[5] = c;
   this->G_tetra_loop[20] = -150;

   /* UGAAAA -150 */
   l = this->tetra_loop[21];
   l[0] = u; l[1] = g; l[2] = a; l[3] = a; l[4] = a; l[5] = a;
   this->G_tetra_loop[21] = -150;
 
   /* AGCAAU -150 */
   l = this->tetra_loop[22];
   l[0] = a; l[1] = g; l[2] = c; l[3] = a; l[4] = a; l[5] = u;
   this->G_tetra_loop[22] = -150;
 
   /* AGUAAU -150 */
   l = this->tetra_loop[23];
   l[0] = a; l[1] = g; l[2] = u; l[3] = a; l[4] = a; l[5] = u;
   this->G_tetra_loop[23] = -150;

   /* CGGGAG -150 */
   l = this->tetra_loop[24];
   l[0] = c; l[1] = g; l[2] = g; l[3] = g; l[4] = a; l[5] = g;
   this->G_tetra_loop[24] = -150;

   /* AGUGAU -150 */
   l = this->tetra_loop[25];
   l[0] = a; l[1] = g; l[2] = u; l[3] = g; l[4] = a; l[5] = u;
   this->G_tetra_loop[25] = -150;

   /* GGCGAC -150 */
   l = this->tetra_loop[26];
   l[0] = g; l[1] = g; l[2] = c; l[3] = g; l[4] = a; l[5] = c;
   this->G_tetra_loop[26] = -150;

   /* GGGAGC -150 */
   l = this->tetra_loop[27];
   l[0] = g; l[1] = g; l[2] = g; l[3] = a; l[4] = g; l[5] = c;
   this->G_tetra_loop[27] = -150;

   /* GUGAAC -150 */
   l = this->tetra_loop[28];
   l[0] = g; l[1] = u; l[2] = g; l[3] = a; l[4] = a; l[5] = c;
   this->G_tetra_loop[28] = -150;
 
   /* UGGAAA -150 */
   l = this->tetra_loop[29];
   l[0] = u; l[1] = g; l[2] = g; l[3] = a; l[4] = a; l[5] = a;
   this->G_tetra_loop[29] = -150;


   /*nn_scores_fprintf_tetra_loop(stderr, this, sigma);
     mfprintf (stderr, "\n\n");*/

   tetra_loop_qsort (0, this->tetra_loop_size, this);

   /*mfprintf (stderr, "\n\n");
   nn_scores_fprintf_tetra_loop(stderr, this, sigma);*/

   return 0;
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
nn_scores_new_init (Alphabet* sigma, const char* file, const int line)
{
   NN_scores* this;
   /*unsigned long i;*/
   char a, u, g, c;

   assert (sigma);

   this = nn_scores_new (file, line);

   if (this == NULL)
   {
      return NULL;
   }


   if (! alphabet_is_standard_rna (sigma))
   {
      nn_scores_delete (this);
      return NULL;
   }
   
   /* fetch nucleotide identifiers from alphabet (just for convenience) */
   a = alphabet_base_2_no('A', sigma);
   u = alphabet_base_2_no('U', sigma);
   g = alphabet_base_2_no('G', sigma);
   c = alphabet_base_2_no('C', sigma);
   
   /* create table of allowed base pairs */
   if (allocate_init_bp_allowed (a, u, g, c, this, file, line))
   {
      nn_scores_delete (this);
      return NULL;
   }

   /* create table of base pair indeces */
   if (allocate_init_bp_idx (alphabet_size (sigma), a, u, g, c, this,
                             file, line))
   {
      nn_scores_delete (this);
      return NULL;
   }
   
   /* create table of stacking energies */
   if (allocate_init_G_stack (a, u, g, c, this, file, line))
   {
      nn_scores_delete (this);
      return NULL;
   }

   /* create table for mismatch stacking energies */
   if(allocate_init_G_mm_stack_size (a, u, g, c, alphabet_size (sigma), this,
                                     file, line))
   {
      nn_scores_delete (this);
      return NULL;       
   }

   /* init hairpin loops */
   if (allocate_init_hairpin_loop (this, file, line))
   {
      nn_scores_delete (this);
      return NULL;
   }
   
   /* init hairpin closing base pair energies */
   if (allocate_init_mismatch_hairpin (a, u, g, c,
                                       alphabet_size (sigma),
                                       this,
                                       file, line))
   {
      nn_scores_delete (this);
      return NULL;
   }

   /* init bulge loops */
   if (allocate_init_bulge_loop (this, file, line))
   {
      nn_scores_delete (this);
      return NULL;
   }

   /* init penalties for non-GC closing base pairs of loops */
   if (allocate_init_non_gc_penalty_for_bp (a, u, g, c, this, file,line))
   {
      nn_scores_delete (this);
      return NULL;      
   }

   /* init tetra loop index and score values */
   if (allocate_init_tetra_loop (a, u, g, c, this, file, line))
   {
      nn_scores_delete (this);
      return NULL;      
   }

   if (allocate_init_dangle5 (a, u, g, c, alphabet_size (sigma), this,
                              file, line))
   {
      nn_scores_delete (this);
      return NULL; 
   }

   if (allocate_init_dangle3 (a, u, g, c, alphabet_size (sigma), this,
                              file, line))
   {
      nn_scores_delete (this);
      return NULL; 
   }

   if (allocate_init_internal_loop (this, file, line))
   {
      nn_scores_delete (this);
      return NULL;
   }

   if (allocate_init_int11 (a, u, g, c, alphabet_size (sigma), this,
                     file, line))
   {
      nn_scores_delete (this);
      return NULL; 
   }

   if (allocate_init_int21 (a, u, g, c, alphabet_size (sigma), this,
                     file, line))
   {
      nn_scores_delete (this);
      return NULL; 
   }

   if (allocate_init_int22 (a, u, g, c, alphabet_size (sigma), this,
                     file, line))
   {
      nn_scores_delete (this);
      return NULL; 
   }

   if (allocate_init_mismatch_interior (a, u, g, c,
                                        alphabet_size (sigma),
                                        this,
                                        file, line))
   {
      nn_scores_delete (this);
      return NULL;
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
     XFREE (this->G_hairpin_loop);
     XFREE_ND (D_MM_H, (void**) this->G_mismatch_hairpin);
     XFREE (this->G_bulge_loop);
     XFREE (this->non_gc_penalty_for_bp);
     XFREE_2D ((void**)this->tetra_loop);
     XFREE (this->G_tetra_loop);
     XFREE (this->G_internal_loop);
     XFREE_ND (D_INT11, (void**) this->G_int11);
     XFREE_ND (D_INT21, (void**) this->G_int21);
     XFREE_ND (D_INT22, (void**) this->G_int22);
     XFREE_ND (D_MM_I, (void**) this->G_mismatch_interior);
     XFREE_2D ((void**)this->G_dangle5);
     XFREE_2D ((void**)this->G_dangle3);
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

int
nn_scores_get_G_extloop_multiloop (const char* seq,
                                   const unsigned long unpaired,
                                   const unsigned long nstems,
                                   unsigned long (*stems)[No_Of_Strands],
                                   const unsigned long ndangle5,
                                   unsigned long (*dangle5)[No_Of_Dangles],
                                   const unsigned long ndangle3,
                                   unsigned long (*dangle3)[No_Of_Dangles],
                                   const bool is_multiloop,
                                   const NN_scores* scheme)
{
   unsigned long i;
   int G = 0;

   assert (scheme);
   assert (scheme->non_gc_penalty_for_bp);
   assert (scheme->bp_idx);
   assert (seq);
   assert (stems);

   /* penalty for non gc basepair initiating a stem */
   for (i = 0; i < nstems; i++)
   {
      G += scheme->non_gc_penalty_for_bp[
         (int)scheme->bp_idx[
            (int)seq[stems[i][P5_Strand]]
                            ]
                            [
            (int)seq[stems[i][P3_Strand]]
                            ]
                                        ];
   }

   /* 5' dangle */
   for (i = 0; i < ndangle5; i++)
   {

      G += scheme->G_dangle5[
         (int)scheme->bp_idx[
            (int)seq[dangle5[i][P5_Dangle]]
                            ]
                            [
            (int)seq[dangle5[i][P3_Dangle]]
                            ]
         ][(int)seq[dangle5[i][Ne_Dangle]]
                            ];
   }

   /* 3' dangle */
   for (i = 0; i < ndangle3; i++)
   {
      G += scheme->G_dangle3[
         (int)scheme->bp_idx[
            (int)seq[dangle3[i][P5_Dangle]]
                            ]
                            [
            (int)seq[dangle3[i][P3_Dangle]]
                            ]
         ][(int)seq[dangle3[i][Ne_Dangle]]
                            ];
   }

   /* linear multiloop energy */
   if (is_multiloop)
   {
      G += NN_ML_OFFSET + unpaired * NN_ML_UNPAIRED + nstems * NN_ML_STEMS;
   }

   return G;
}

/** @brief Return the stacking score for a set of paired bases.
 *
 * @params[in] i i component of the upstream pair of the stack (5' end).
 * @params[in] j j component of the upstream pair of the stack (pairs with i).
 * @params[in] jm1 j-1 component of the downstream pair (pairs with ip1).
 * @params[in] ip1 i+1 component of the downstream pair.
 * @params[in] scheme The scoring scheme.
 */
long
nn_scores_get_G_stack (const char i, const char j,
                       const char jm1, const char ip1,
                       const NN_scores* scheme)
{
   assert (scheme);
   assert (scheme->G_stack);
   assert (scheme->bp_idx);
   assert ((unsigned)   i < sqrtf (scheme->bp_idx_size));
   assert ((unsigned)   j < sqrtf (scheme->bp_idx_size));
   assert ((unsigned) ip1 < sqrtf (scheme->bp_idx_size));
   assert ((unsigned) jm1 < sqrtf (scheme->bp_idx_size));
   assert (  (unsigned) scheme->bp_idx[(int)i][(int)j] 
           < sqrtf (scheme->G_stack_size));
   assert (  (unsigned) scheme->bp_idx[(int)jm1][(int)ip1]
           < sqrtf (scheme->G_stack_size));
   
   return scheme->G_stack[(int) scheme->bp_idx[(int)i][(int)j]]
                         [(int) scheme->bp_idx[(int)jm1][(int)ip1]];
}

/** @brief Return the mismatch stacking score for a set of bases.
 *
 * @params[in] i i component of the upstream pair of the stack (5' end).
 * @params[in] j j component of the upstream pair of the stack (pairs with i).
 * @params[in] k position j-1.
 * @params[in] l position i+1.
 * @params[in] scheme The scoring scheme.
 */
long
nn_scores_get_G_mm_stack (const char i, const char j,
                          const char k, const char l,
                          const NN_scores* scheme)
{
   assert (scheme);
   assert (scheme->G_mm_stack);
   assert (scheme->bp_idx);
   assert ((unsigned) i < sqrtf (scheme->bp_idx_size));
   assert ((unsigned) j < sqrtf (scheme->bp_idx_size));
   assert ((unsigned) k < sqrtf (scheme->bp_idx_size));
   assert ((unsigned) l < sqrtf (scheme->bp_idx_size));
   assert (  (unsigned) scheme->bp_idx[(int)i][(int)j] 
           < scheme->bp_allowed_size);
   assert (  (unsigned) scheme->bp_idx[(int)k][(int)l]
             < (scheme->G_mm_stack_size / scheme->bp_allowed_size));
   
   return scheme->G_mm_stack[(int) scheme->bp_idx[(int)i][(int)j]]
                            [(int) scheme->bp_idx[(int)k][(int)l]];
}

static __inline__ int
tetra_loop_cmp_seq (const char* seq,
                    const unsigned long i,
                    unsigned long loop,
                    const NN_scores* this)
{
   unsigned long k;

   for (k = 0; k < D_TL; k++)
   {
      if (this->tetra_loop[loop][k] != seq[i + k])
      {
         return this->tetra_loop[loop][k] - seq[i + k];
      }
   }

   return 0;
}

/** @brief Return the bonus score for a tetra loop.
 *
 * @param[in] seq transformed RNA sequence.
 * @param[in] i i component of the closing base pair.
 * @param[in] j j component of the closing base pair.
 * @param[in] this scoring sceme.
 */
int
nn_scores_get_G_tetra_loop (const char* seq,
                            const unsigned long i,
                            const NN_scores* this)
{
   unsigned long l, r, m; /* search interval */
   int cmp;

   assert (seq);
   assert (this);

   l = 0;
   r = this->tetra_loop_size;

   m = 0;
   cmp = 1;
   while (l < r)
   {
      m = (l + r) / 2;

      /* mfprintf (stderr, "Search l:%lu r:%lu m:%lu\n", l, r, m); */

      cmp = tetra_loop_cmp_seq (seq, i, m, this);
      if (cmp < 0)
      {
         l = m + 1;
      }
      else
      {
         r = m;
      }

   }

   /* mfprintf (stderr, "Search l:%lu r:%lu m:%lu\n", l, r, m); */

   if ((l < this->tetra_loop_size)&&(tetra_loop_cmp_seq (seq, i, l, this) == 0))
   {
     /*  mfprintf (stderr, "FOUND\n"); */
      return this->G_tetra_loop[l];
   }

   return 0;
}

/** @brief Returns the score for a hairpin loop of certain size.
 *
 * @params[in] i 5' base of closing pair.
 * @params[in] j 3' base of closing pair.
 * @params[in] size Length of the loop (unpaired bases only).
 * @params[in] this Scoring scheme.
 */
int
nn_scores_get_G_hairpin_loop (const char* seq,
                              const unsigned long i,
                              const unsigned long j,
                              const unsigned long size,
                              const NN_scores* this)
{
   int G = 0;
   int bp;
   int bip1 = seq[i + 1];
   int bjm1 = seq[j - 1];

   assert (seq);
   assert (this);
   assert (j > 0);

   bp = (int)this->bp_idx[(int)seq[i]][(int)seq[j]];

   if (size < this->G_hairpin_loop_size)
   {
      G += this->G_hairpin_loop[size];
   }
   else
   {
      G += this->G_hairpin_loop[this->G_hairpin_loop_size - 1]
         + (NN_LXC37 * logf((float) size / (this->G_hairpin_loop_size - 1)));
   }

    /* mismatch penalty for the mismatch interior to the closing basepair of
       the hairpin. triloops get non-parameterised mismatch penalty */
   if (size == D_MM_H)
   {
      G += this->non_gc_penalty_for_bp[bp];
   }
   else
   {
      G += this->G_mismatch_hairpin[bp][bip1][bjm1];
   }
     /*  mfprintf (stderr, "G= %d\n", G); */
   /* tetraloop bonus */
   if (size == 4)
   {
      G += nn_scores_get_G_tetra_loop (seq, i, this);
   }
   /* mfprintf (stderr, "G= %d\n", G); */
   return G;
}

/** @brief Returns the score for a bulge loop.
 *
 * @params[in] i 5' base of closing pair.
 * @params[in] j 3' base of closing pair.
 * @params[in] size Length of the loop (unpaired bases only).
 * @params[in] this Scoring scheme.
 */
int
nn_scores_get_G_bulge_loop (const char* seq,
                            const unsigned long i1,
                            const unsigned long j1,
                            const unsigned long i2,
                            const unsigned long j2,
                            const unsigned long size,
                            const NN_scores* this)
{
   int G = 0;

   assert (seq);
   assert (this);
   assert (this->G_bulge_loop);
   assert (this->non_gc_penalty_for_bp);
   assert (i1 < j1);
   assert (i2 < j2);

   if (size < this->G_bulge_loop_size)
   {
      G += this->G_bulge_loop[size];
   }
   else
   {
      G += this->G_bulge_loop[this->G_bulge_loop_size - 1]
         + (int) (NN_LXC37 *
                  logf((float) size / (this->G_bulge_loop_size - 1)));
   }

   if (size == 1)
   {
      G += nn_scores_get_G_stack (seq[i1], seq[j1], seq[j2], seq[i2], this);
   }
   else
   {
      /* bulge loops larger than 1 get penalty term for non-gc closing
         basepairs */
      G += this->non_gc_penalty_for_bp[
         (int)this->bp_idx[(int)seq[i1]][(int)seq[j1]]];
      G += this->non_gc_penalty_for_bp[
         (int)this->bp_idx[(int)seq[j2]][(int)seq[i2]]];
   }

   return G;
}

int
nn_scores_get_G_internal_loop (const char* seq,
                               const unsigned long size1,
                               const unsigned long size2,
                               const unsigned long i1,
                               const unsigned long j1,
                               const unsigned long i2,
                               const unsigned long j2,
                               const NN_scores* this)
{
   int G = 0;
   int bp1, bp2;
   int bi1p, bi2m, bj2p, bj1m;  /* bi1p = seq[i1 + 1] */
   unsigned long size;

   assert (seq);
   assert (this);
   assert (this->bp_idx);
   assert (this->G_int11);
   assert (this->G_int21);
   assert (this->G_int22);
   assert (this->G_mismatch_interior);
   assert (i1 < j1);
   assert (i1 < i2);
   assert (i2 < j2);
   assert (j2 < j1);

   bp1 =  (int)this->bp_idx[(int)seq[i1]][(int)seq[j1]];
   bp2 =  (int)this->bp_idx[(int)seq[j2]][(int)seq[i2]];
   bi1p = (int)seq[i1 + 1];
   bi2m = (int)seq[i2 - 1];
   bj2p = (int)seq[j2 + 1];
   bj1m = (int)seq[j1 - 1];

   if ((size1 == 1) && (size2 == 1))
   {
      /* 1x1 internal loop */
      return this->G_int11[bp1][bp2][bi1p][bj2p];
   }
   else if ((size1 == 1) && (size2 == 2))
   {
      /* 1x2 internal loop */
      return this->G_int21[bp1][bp2][bi1p][bj2p][bj1m];
   }
   else if ((size1 == 2) && (size2 == 1))
   {
      /* 2x1 internal loop */
      /* note switched order of pt1 and pt2 compared to 1x2 loop */
      return this->G_int21[bp2][bp1][bj2p][bi1p][bi2m];
   }
   else if ((size1 == 2) && (size2 == 2))
   {
      /* 2x2 internal loop */
      return this->G_int22[bp1][bp2][bi1p][bi2m][bj2p][bj1m];
   }
   else
   {
      /* generic internal loop */
      size = size1 + size2;
      if (size < this->G_internal_loop_size)
      {
         G += this->G_internal_loop[size];
      }
      else
      {
         G += this->G_internal_loop[this->G_internal_loop_size - 1]
            + NN_LXC37 * logf ((float) size / (this->G_internal_loop_size - 1));
      }
      /* loop asymmetry contribution */
      G += (NN_NINIO_MAX < (labs (size1 - size2) * NN_NINIO_M) ? NN_NINIO_MAX
            : (labs (size1 - size2) * NN_NINIO_M));
      /* mismatch contribution */
      G += this->G_mismatch_interior[bp1][bi1p][bj1m];
      G += this->G_mismatch_interior[bp2][bj2p][bi2m];
   }

   return G;
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
   unsigned long pline_width = 2;
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
      /* row id */
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

/** @brief Print the hairpin loop energies of a scoring scheme to a stream.
 *
 * Prints hairpin loop energies to a stream. Form is "loop size: score".
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] scheme The scoring scheme.
 */
void
nn_scores_fprintf_G_hairpin_loop (FILE* stream, const NN_scores* scheme)
{
   unsigned long i;
   int rprec, rprec_idx, tmp, pline_width = 0;
   char* string;
   char* string_start;
   char* en_undef;

   assert (scheme != NULL);
   assert (scheme->G_hairpin_loop != NULL);

   /* dermine widest cell */
   for (i = 0; i < scheme->G_hairpin_loop_size; i++)
   {
      rprec = 0;
      tmp = scheme->G_hairpin_loop[i];
      if (tmp == INT_UNDEF)
      {
         tmp = 0;
      }

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
      
      if (rprec > pline_width) 
      {
         pline_width = rprec;
      }       
   }
   rprec = pline_width;
   rprec_idx = floor (log10 (scheme->G_hairpin_loop_size) + 1);

   /* add up components of a line */
   pline_width += 3; /*: \n*/
   pline_width += rprec_idx;

   /* allocate memory for the string and the undef symbol
      therefore we add "rprec + 1" which is the 2 (for null terminating the
      strings) */
   string = XMALLOC (sizeof (char) *
                     (pline_width * scheme->G_hairpin_loop_size) + 2 + rprec);
   if (string == NULL)
   {
      return;
   }

   en_undef = string;
   for (i = 0; i < (unsigned long) rprec; i++)
   {
      en_undef[i] = '-';
   }
   en_undef[i] = '\0';


   string_start = string + rprec + 1;
   string = string + rprec + 1;

   /* start printing */
   for (i = 0; i < scheme->G_hairpin_loop_size; i++)
   {
      /* store index */
      msprintf (string, "%*ld", rprec_idx, i);
      string+= rprec_idx;

      msprintf (string, ": ");
      string+= 2;

      if (scheme->G_hairpin_loop[i] == INT_UNDEF)
      {
         msprintf (string, "%s", en_undef);         
      }
      else
      {
         msprintf (string, "%*i", rprec, scheme->G_hairpin_loop[i]);
      }
      string+= rprec;

      string[0] = '\n';
      string++;
   }

   /* actually print */
   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   XFREE (en_undef);  /* en_undef is the start of the whole memory */
}

/** @brief Print the mismatch hairpin energies of a scoring scheme to a stream.
 *
 * Prints hairpin loop closing base pair energies to a stream.
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] scheme The scoring scheme.
 * @params[in] sigma The alphabet.
 */
void
nn_scores_fprintf_G_mismatch_hairpin (FILE* stream,
                                      const NN_scores* scheme,
                                      const Alphabet* sigma)
{
   unsigned long i, j, k;
   int tmp;
   int rprec;
   unsigned long pline_width = 2;
   unsigned long x;             /* 3' base of loop pair */
   unsigned long y;             /* 5' base of loop pair */
   unsigned long z;             /* stacked base pair */
   char* string;
   char* string_start;
   char* header;

   assert (scheme != NULL);
   assert (scheme->G_mismatch_hairpin != NULL);
   assert (scheme->bp_allowed != NULL);
   assert (sigma != NULL);

   x = y = alphabet_size (sigma);
   z = scheme->bp_allowed_size;

   /* dermine widest cell */
   for (i = 0; i < z; i++)      /* for all base pairs */
   {
      for (j = 0; j < x; j++)   /* for all 3' bases */
      {
         for (k = 0; k < y; k++) /* for all 5' bases */
         {
            rprec = 0;
            tmp = scheme->G_mismatch_hairpin[i][j][k];

            /* fetch '-' symbol */
            if (tmp < 0)
            {
               tmp *= (-1);
               rprec++;
            }

            /* get no. of digits */
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
   }
   rprec = pline_width;

   /* add up components of a line */
   pline_width += 3;            /*\s|\s*/
   pline_width *= x;
   pline_width += 1;            /* + N */
   pline_width += 1;            /* + \n */

   /* create header line for each table */
   string = (char*) XMALLOC (sizeof (char) * pline_width);
   if (string == NULL)
   {
      return;
   }
   header = string;

   if (rprec > 1)
   {
      rprec--;
   }

   for (i = 0; i < x; i++)
   {
      msprintf (string, "  | %*c", rprec, alphabet_no_2_base (i, sigma));
      string += rprec;
      string += 4;
   }
   string[0] = '\n';
   string++;
   string[0] = '\0';

   if ((rprec + 1) > 1)
   {
      rprec++;
   }

   /* now for the tables */
   /* add to pline_width: room for base pair + table head */
   pline_width *= y;
   pline_width += strlen (header); /* table header */
   pline_width += 4;               /* NN:\n */
   string = (char*) XMALLOC (sizeof (char) * ((pline_width * z) + 1));
   if (string == NULL)
   {
      return;
   }
   string_start = string;  
   
   for (i = 0; i < z; i++)
   {
      /* start new bp table */
      msprintf (string, "%c%c:\n", 
                alphabet_no_2_base (scheme->bp_allowed[i][0], sigma),
                alphabet_no_2_base (scheme->bp_allowed[i][1], sigma));
      string += 4;
      msprintf (string, "%s", header);
      string += strlen (header);

      for (j = 0; j < y; j++)
      {
         msprintf (string, "%c", alphabet_no_2_base (j, sigma));
         string++;
         for (k = 0; k < x; k++)
         {
            msprintf (string, " | %*i", rprec,
                      scheme->G_mismatch_hairpin[i][j][k]);
            string += rprec;
            string += 3;
         }
         string[0] = '\n';
         string++;
      }
   }

   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   XFREE (string_start);
   XFREE (header);
}

/** @brief Print the bulge loop energies of a scoring scheme to a stream.
 *
 * Prints bulge loop energies to a stream. Form is "loop size: score".
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] scheme The scoring scheme.
 */
void
nn_scores_fprintf_G_bulge_loop (FILE* stream, const NN_scores* scheme)
{
   unsigned long i;
   int rprec, rprec_idx, tmp, pline_width = 0;
   char* string;
   char* string_start;
   char* en_undef;

   assert (scheme != NULL);
   assert (scheme->G_bulge_loop != NULL);

   /* dermine widest cell */
   for (i = 0; i < scheme->G_bulge_loop_size; i++)
   {
      rprec = 0;
      tmp = scheme->G_bulge_loop[i];
      if (tmp == INT_UNDEF)
      {
         tmp = 0;
      }

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
      
      if (rprec > pline_width) 
      {
         pline_width = rprec;
      }       
   }
   rprec = pline_width;
   rprec_idx = floor (log10 (scheme->G_bulge_loop_size) + 1);

   /* add up components of a line */
   pline_width += 3; /*: \n*/
   pline_width += rprec_idx;

   /* allocate memory for the string and the undef symbol
      therefore we add "rprec + 1" which is the 2 (for null terminating the
      strings) */
   string = XMALLOC (sizeof (char) *
                     (pline_width * scheme->G_bulge_loop_size) + 2 + rprec);
   if (string == NULL)
   {
      return;
   }

   en_undef = string;
   for (i = 0; i < (unsigned long) rprec; i++)
   {
      en_undef[i] = '-';
   }
   en_undef[i] = '\0';


   string_start = string + rprec + 1;
   string = string + rprec + 1;

   /* start printing */
   for (i = 0; i < scheme->G_bulge_loop_size; i++)
   {
      /* store index */
      msprintf (string, "%*ld", rprec_idx, i);
      string+= rprec_idx;

      msprintf (string, ": ");
      string+= 2;

      if (scheme->G_bulge_loop[i] == INT_UNDEF)
      {
         msprintf (string, "%s", en_undef);         
      }
      else
      {
         msprintf (string, "%*i", rprec, scheme->G_bulge_loop[i]);
      }
      string+= rprec;

      string[0] = '\n';
      string++;
   }

   /* actually print */
   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   XFREE (en_undef);  /* en_undef is the start of the whole memory */
}


/** @brief Print the penalties for non-GC closing base pairs.
 *
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] scheme The scoring scheme.
 * @params[in] sigma The alphabet.
 */
void
nn_scores_fprintf_non_gc_penalty_for_bp(FILE* stream,
                                        const NN_scores* scheme,
                                        const Alphabet* sigma)
{
   unsigned long i;
   int tmp;
   int rprec;
   unsigned long pline_width = 2;
   char* string;
   char* string_start;

   assert (scheme != NULL);
   assert (scheme->non_gc_penalty_for_bp != NULL);
   assert (scheme->bp_allowed != NULL);
   assert (sigma != NULL);

   /* determine highest no. of digits */
   for (i = 0; i < scheme->bp_allowed_size; i++) /* for all base pairs */
   {
      rprec = 0;
      tmp = scheme->non_gc_penalty_for_bp[i];
      
      /* fetch '-' symbol */
      if (tmp < 0)
      {
         tmp *= (-1);
         rprec++;
      }

      /* get no. of digits */
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
   rprec = pline_width;

   /* add up components of a line */
   pline_width += 4;            /*NN:\s*/
   pline_width += 1;            /* + \n */

   string = (char*) XMALLOC (sizeof (char)
                             * ((pline_width * scheme->bp_allowed_size) + 1));
   if (string == NULL)
   {
      return;
   }
   string_start = string;  
   
   for (i = 0; i < scheme->bp_allowed_size; i++)
   {
      msprintf (string, "%c%c: ", 
                alphabet_no_2_base (scheme->bp_allowed[i][0], sigma),
                alphabet_no_2_base (scheme->bp_allowed[i][1], sigma));
      string += 4;

      msprintf (string, "%*i", rprec, scheme->non_gc_penalty_for_bp[i]);
      string += rprec;

      string[0] = '\n';
      string++;
   }

   string[0] = '\0';
   mfprintf (stream, "%s", string_start);  

   XFREE (string_start);
}

/** @brief Print the bonus scores for tetra loops.
 *
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] scheme The scoring scheme.
 * @params[in] sigma The alphabet.
 */
void
nn_scores_fprintf_tetra_loop(FILE* stream,
                             const NN_scores* scheme,
                             const Alphabet* sigma)
{
   unsigned long i, j;
   int tmp;
   int rprec;
   unsigned long pline_width = 2;
   char* string;
   char* string_start;

   assert (scheme != NULL);
   assert (scheme->tetra_loop != NULL);
   assert (scheme->G_tetra_loop != NULL);
   assert (sigma != NULL);

   /* find largest no. */
   for (i = 0; i < scheme->tetra_loop_size; i++)
   {
      rprec = 0;
      tmp = scheme->G_tetra_loop[i];
      
      /* fetch '-' symbol */
      if (tmp < 0)
      {
         tmp *= (-1);
         rprec++;
      }
      
      /* get no. of digits */
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
   rprec = pline_width;   

   /* add up components of a line */
   pline_width += 10;            /*N-NNNN-N:\s*/
   pline_width += 1;            /* + \n */

   /* allocate buffer */
   string = (char*) XMALLOC (sizeof (char)
                             * ((pline_width * scheme->tetra_loop_size) + 1));
   if (string == NULL)
   {
      return;
   }
   string_start = string;  

   /* print */
   for (i = 0; i < scheme->tetra_loop_size; i++)
   {
      msprintf (string, "%c-", alphabet_no_2_base (scheme->tetra_loop[i][0],
        sigma));
      /*msprintf (string, "%i-", scheme->tetra_loop[i][0]);*/
      string += 2;      
      for (j = 1; j < (D_TL - 1); j++)
      {
         msprintf (string, "%c", alphabet_no_2_base (scheme->tetra_loop[i][j],
           sigma));
         /*msprintf (string, "%i", scheme->tetra_loop[i][j]);*/
         string ++;
      }
      msprintf (string, "-%c", alphabet_no_2_base (scheme->tetra_loop[i][j],
        sigma));
      /*msprintf (string, "-%i", scheme->tetra_loop[i][j]);*/
      string += 2; 

      msprintf (string, ": %*i", rprec, scheme->G_tetra_loop[i]);
      string += rprec;
      string += 2;

      string[0] = '\n';
      string++;
   }

   string[0] = '\0';
   mfprintf (stream, "%s", string_start);  

   XFREE (string_start);
}

/** @brief Print 5' dangling end scores.
 *
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] scheme The scoring scheme.
 * @params[in] sigma The alphabet.
 */
void
nn_scores_fprintf_G_dangle5(FILE* stream,
                            const NN_scores* scheme,
                            const Alphabet* sigma)
{
   unsigned long i, j;
   int tmp;
   int rprec;
   unsigned long pline_width = 1;
   unsigned long columns;
   unsigned long rows;
   char* string;
   char* string_start;

   assert (scheme != NULL);
   assert (scheme->G_dangle5 != NULL);
   assert (scheme->bp_allowed != NULL);
   assert (sigma != NULL);

   rows = scheme->bp_allowed_size;
   columns = alphabet_size (sigma);

   /* dermine widest cell */
   for (i = 0; i < rows; i++)
   {
      for (j = 0; j < columns; j++)
      {
         rprec = 0;
         tmp = scheme->G_dangle5[i][j];
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
   pline_width *= columns;
   pline_width += 2;            /* + NN */
   pline_width += 1;            /* + \n */

   /* alloc memory for the string */
   string = XMALLOC (sizeof (char) * ((pline_width * (rows + 1)) + 1));
   if (string == NULL)
   {
      return;
   }

   string_start = string;
   
   /* print base pairs horizontally */
   msprintf (string, "  "); /*skip first "NN"*/
   string += 2;

   for (i = 0; i < columns; i++)
   {
      msprintf (string, " | %*c", rprec, alphabet_no_2_base (i, sigma));
      string += (rprec + 3);
   }

   string[0] = '\n';
   string++;

   /* print matrix */
   for (i = 0; i < rows; i++)
   {
      /* row id */
      msprintf (string, "%c%c",
                alphabet_no_2_base (scheme->bp_allowed[i][0], sigma),
                alphabet_no_2_base (scheme->bp_allowed[i][1], sigma));
      string += 2;
      
      for (j = 0; j < columns; j++)
      {
         msprintf (string, " | ");
         string += 3;
         
         msprintf (string, "%*d", rprec,
                   scheme->G_dangle5[i][j]);
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

/** @brief Print 3' dangling end scores.
 *
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] scheme The scoring scheme.
 * @params[in] sigma The alphabet.
 */
void
nn_scores_fprintf_G_dangle3(FILE* stream,
                            const NN_scores* scheme,
                            const Alphabet* sigma)
{
   unsigned long i, j;
   int tmp;
   int rprec;
   unsigned long pline_width = 1;
   unsigned long columns;
   unsigned long rows;
   char* string;
   char* string_start;

   assert (scheme != NULL);
   assert (scheme->G_dangle3 != NULL);
   assert (scheme->bp_allowed != NULL);
   assert (sigma != NULL);

   rows = scheme->bp_allowed_size;
   columns = alphabet_size (sigma);

   /* dermine widest cell */
   for (i = 0; i < rows; i++)
   {
      for (j = 0; j < columns; j++)
      {
         rprec = 0;
         tmp = scheme->G_dangle3[i][j];
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
   pline_width *= columns;
   pline_width += 2;            /* + NN */
   pline_width += 1;            /* + \n */

   /* alloc memory for the string */
   string = XMALLOC (sizeof (char) * ((pline_width * (rows + 1)) + 1));
   if (string == NULL)
   {
      return;
   }

   string_start = string;
   
   /* print base pairs horizontally */
   msprintf (string, "  "); /*skip first "NN"*/
   string += 2;

   for (i = 0; i < columns; i++)
   {
      msprintf (string, " | %*c", rprec, alphabet_no_2_base (i, sigma));
      string += (rprec + 3);
   }

   string[0] = '\n';
   string++;

   /* print matrix */
   for (i = 0; i < rows; i++)
   {
      /* row id */
      msprintf (string, "%c%c",
                alphabet_no_2_base (scheme->bp_allowed[i][0], sigma),
                alphabet_no_2_base (scheme->bp_allowed[i][1], sigma));
      string += 2;
      
      for (j = 0; j < columns; j++)
      {
         msprintf (string, " | ");
         string += 3;
         
         msprintf (string, "%*d", rprec,
                   scheme->G_dangle3[i][j]);
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

static int
get_ndigits_matrix (const unsigned long cols,
                    const unsigned long rows,
                    int** matrix)
{
   int tmp = 0;
   int cval;
   unsigned long i, j;

   /* find largest non-negative no. */
   for (i = 0; i < rows; i++)
   {
      for (j = 0; j < cols; j++)
      {
         if (matrix[i][j] == 0)
         {
            cval = 0;
         }
         else if (matrix[i][j] < 0)
         {
            cval = floor (log10 (matrix[i][j] * (-1)) + 2);
         }
         else
         {
            cval = floor (log10 (matrix[i][j]) + 1);
         }
         
         if (cval > tmp)
         {
            tmp = cval;
         }
      } 
   }

   if (tmp == 0)
   {
      return 1;
   }

   return tmp;
}

/** @brief Print the internal loop energies of a scoring scheme to a stream.
 *
 * Prints internal loop energies to a stream. Form is "loop size: score".
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] scheme The scoring scheme.
 */
void
nn_scores_fprintf_G_internal_loop (FILE* stream, const NN_scores* scheme)
{
   unsigned long i;
   int rprec, rprec_idx, tmp, pline_width = 0;
   char* string;
   char* string_start;
   char* en_undef;

   assert (scheme != NULL);
   assert (scheme->G_internal_loop != NULL);

   /* dermine widest cell */
   for (i = 0; i < scheme->G_internal_loop_size; i++)
   {
      rprec = 0;
      tmp = scheme->G_internal_loop[i];
      if (tmp == INT_UNDEF)
      {
         tmp = 0;
      }

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
      
      if (rprec > pline_width) 
      {
         pline_width = rprec;
      }       
   }
   rprec = pline_width;
   rprec_idx = floor (log10 (scheme->G_internal_loop_size) + 1);

   /* add up components of a line */
   pline_width += 3; /*: \n*/
   pline_width += rprec_idx;

   /* allocate memory for the string and the undef symbol
      therefore we add "rprec + 1" which is the 2 (for null terminating the
      strings) */
   string = XMALLOC (sizeof (char) *
                     (pline_width * scheme->G_internal_loop_size) + 2 + rprec);
   if (string == NULL)
   {
      return;
   }

   en_undef = string;
   for (i = 0; i < (unsigned long) rprec; i++)
   {
      en_undef[i] = '-';
   }
   en_undef[i] = '\0';


   string_start = string + rprec + 1;
   string = string + rprec + 1;

   /* start printing */
   for (i = 0; i < scheme->G_internal_loop_size; i++)
   {
      /* store index */
      msprintf (string, "%*ld", rprec_idx, i);
      string+= rprec_idx;

      msprintf (string, ": ");
      string+= 2;

      if (scheme->G_internal_loop[i] == INT_UNDEF)
      {
         msprintf (string, "%s", en_undef);         
      }
      else
      {
         msprintf (string, "%*i", rprec, scheme->G_internal_loop[i]);
      }
      string+= rprec;

      string[0] = '\n';
      string++;
   }

   /* actually print */
   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   XFREE (en_undef);  /* en_undef is the start of the whole memory */
}

/** @brief Print the parameters for 1x1 internal loops to a stream.
 *
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] scheme The scoring scheme.
 * @params[in] sigma The alphabet.
 */
void
nn_scores_fprintf_G_int11 (FILE* stream,
                           const NN_scores* scheme,
                           const Alphabet* sigma)
{
   int rprec = 0;
   int tmp;
   unsigned long i, j, k, l;
   unsigned long store_size = 0, no = 0, asize = 0, hsize = 0;
   char* string;
   char* header;
   char* string_start;

   assert (scheme != NULL);
   assert (scheme->G_int11 != NULL);
   assert (scheme->bp_allowed != NULL);
   assert (sigma != NULL);

   /* get largest no. of digits of all table entries */
   asize = alphabet_size (sigma);
   for (i = 0; i < scheme->bp_allowed_size; i++)
   {
      for (j = 0; j < scheme->bp_allowed_size; j++)
      {
         tmp = get_ndigits_matrix (asize, asize, scheme->G_int11[i][j]);
         if (tmp > rprec)
         {
            rprec = tmp;
         }
      }
   }

   /* calc. space needs */
   /* 1st D: bp "NN:\n" */
   no = scheme->bp_allowed_size;
   store_size += (4 * no);

   /* 2nd D: bp "  NN:\n" */
   no *= no;
   store_size += (6 * no);

   /* 3rd D: b: "       | A | C | G | U\n" */
   hsize = (asize * rprec) + 6 + (3 * asize); 
   store_size += (hsize * no);

   header = XMALLOC (sizeof (char) * (hsize + 1));
   if (header == NULL)
   {
      return;
   }
   string_start = header;

   msprintf (header, "     ");
   header += 5;
   for (i = 0; i < asize; i++)
   {
      msprintf (header, " | %*c", rprec, alphabet_no_2_base(i, sigma));
      header += rprec + 3;
   }
   msprintf (header, "\n");
   header++;
   header[0] = '\0';
   header = string_start;

   /* 4th D: b: "    N | ...\n" */
   store_size += (5 * asize * no);     /* "    N"  */
   store_size += (asize * no);         /* "\n" */
   no *= (asize * asize);
   store_size += (3 * no);             /* " | " */
   store_size += (rprec * no);

   /* alloc memory for the string */
   string = XMALLOC (sizeof (char) * (store_size + 1));
   if (string == NULL)
   {
      return;
   }
   string_start = string;

   /* print to string*/
   for (i = 0; i < scheme->bp_allowed_size; i++)
   {
      msprintf (string, "%c%c:\n",
                           alphabet_no_2_base(scheme->bp_allowed[i][0], sigma),
                           alphabet_no_2_base(scheme->bp_allowed[i][1], sigma));
      string += 4;
      
      for (j = 0; j < scheme->bp_allowed_size; j++)
      {
         msprintf (string, "  %c%c:\n",
                           alphabet_no_2_base(scheme->bp_allowed[j][0], sigma),
                           alphabet_no_2_base(scheme->bp_allowed[j][1], sigma));
         string += 6;
         
         msprintf (string, "%s", header);
         string += hsize;
         
         for (k = 0; k < asize; k++)
         {
            msprintf (string, "    %c", alphabet_no_2_base(k, sigma));
            string += 5;

            for (l = 0; l < asize; l++)
            {
               msprintf (string, " | %*i", rprec, scheme->G_int11[i][j][k][l]);
               string += 3 + rprec;
            }

            msprintf (string, "\n");
            string += 1;
         }
      }
   }   

   /* print matrix */
   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   XFREE (string_start);
   XFREE (header);
}

/** @brief Print the parameters for 2x1 internal loops to a stream.
 *
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] scheme The scoring scheme.
 * @params[in] sigma The alphabet.
 */
void
nn_scores_fprintf_G_int21 (FILE* stream,
                           const NN_scores* scheme,
                           const Alphabet* sigma)
{
   int rprec = 0;
   int tmp;
   unsigned long i, j, k, l, m;
   unsigned long store_size = 0, no = 0, asize = 0, hsize = 0;
   char* string;
   char* header;
   char* string_start;

   assert (scheme != NULL);
   assert (scheme->G_int21 != NULL);
   assert (scheme->bp_allowed != NULL);
   assert (sigma != NULL);

   /* get largest no. of digits of all table entries */
   asize = alphabet_size (sigma);
   for (i = 0; i < scheme->bp_allowed_size; i++)
   {
      for (j = 0; j < scheme->bp_allowed_size; j++)
      {
         for (k = 0; k < asize; k++)
         {
            tmp = get_ndigits_matrix (asize, asize, scheme->G_int21[i][j][k]);
            if (tmp > rprec)
            {
               rprec = tmp;
            }
         }
      }
   }

   /* calc. space needs */
   /* 1st D: bp "NN:\n" */
   no = scheme->bp_allowed_size;
   store_size += (4 * no);

   /* 2nd D: bp "  NN:\n" */
   no *= no;
   store_size += (6 * no);

   /* 3rd D: b "    N:\n" */
   no *= asize;
   store_size += (7 * no);

   /* 4th D: b: "        | A | C | G | U\n" */
   hsize = (asize * rprec) + 7 + (3 * asize);
   store_size += (hsize * no);

   header = XMALLOC (sizeof (char) * (hsize + 1));
   if (header == NULL)
   {
      return;
   }
   string_start = header;

   msprintf (header, "      ");
   header += 6;
   for (i = 0; i < asize; i++)
   {
      msprintf (header, " | %*c", rprec, alphabet_no_2_base(i, sigma));
      header += rprec + 3;
   }
   msprintf (header, "\n");
   header++;
   header[0] = '\0';
   header = string_start;

   /* 4th D: b: "     N | ...\n" */
   store_size += (6 * asize * no);     /* "     N"  */
   store_size += (asize * no);         /* "\n" */
   no *= (asize * asize);
   store_size += (3 * no);             /* " | " */
   store_size += (rprec * no);

   /* alloc memory for the string */
   string = XMALLOC (sizeof (char) * (store_size + 1));
   if (string == NULL)
   {
      return;
   }
   string_start = string;

   /* print to string*/
   for (i = 0; i < scheme->bp_allowed_size; i++)
   {
      msprintf (string, "%c%c:\n",
                alphabet_no_2_base(scheme->bp_allowed[i][0], sigma),
                alphabet_no_2_base(scheme->bp_allowed[i][1], sigma));
      string += 4;
      
      for (j = 0; j < scheme->bp_allowed_size; j++)
      {
         msprintf (string, "  %c%c:\n",
                           alphabet_no_2_base(scheme->bp_allowed[j][0], sigma),
                           alphabet_no_2_base(scheme->bp_allowed[j][1], sigma));
         string += 6;
         
         for (k = 0; k < asize; k++)
         {
            msprintf (string, "    %c:\n", alphabet_no_2_base(k, sigma));
            string += 7;     
            msprintf (string, "%s", header);
            string += hsize;
         
            for (l = 0; l < asize; l++)
            {
               msprintf (string, "     %c", alphabet_no_2_base(l, sigma));
               string += 6;
               
            for (m = 0; m < asize; m++)
            {
               msprintf (string, " | %*i",
                         rprec, scheme->G_int21[i][j][k][l][m]);
               string += 3 + rprec;
            }

            msprintf (string, "\n");
            string += 1;
            }
         }
      }
   }

   /* print matrix */
   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   XFREE (string_start);
   XFREE (header);
}

void
nn_scores_fprintf_G_int22 (FILE* stream,
                           const NN_scores* scheme,
                           const Alphabet* sigma)
{
   int rprec = 0;
   int tmp;
   unsigned long i, j, k, l, m, n;
   unsigned long store_size = 0, no = 0, asize = 0, hsize = 0;
   char* string;
   char* header;
   char* string_start;

   assert (scheme != NULL);
   assert (scheme->G_int22 != NULL);
   assert (scheme->bp_allowed != NULL);
   assert (sigma != NULL);

   /* get largest no. of digits of all table entries */
   asize = alphabet_size (sigma);
   for (i = 0; i < scheme->bp_allowed_size; i++)
   {
      for (j = 0; j < scheme->bp_allowed_size; j++)
      {
         for (k = 0; k < asize; k++)
         {
            for (l = 0; l < asize; l++)
            {
               tmp = get_ndigits_matrix (asize,
                                         asize,
                                         scheme->G_int22[i][j][k][l]);
               if (tmp > rprec)
               {
                  rprec = tmp;
               }
            }
         }
      }
   }

   /* calc. space needs */
   /* 1st D: bp "NN:\n" */
   no = scheme->bp_allowed_size;
   store_size += (4 * no);

   /* 2nd D: bp "  NN:\n" */
   no *= no;
   store_size += (6 * no);

   /* 3rd D: b "    N:\n" */
   no *= asize;
   store_size += (7 * no);

   /* 4th D: b "     N:\n" */
   no *= asize;
   store_size += (8 * no);

   /* 5th D: b: "        | A | C | G | U\n" */
   hsize = (asize * rprec) + 8 + (3 * asize);
   store_size += (hsize * no);

   header = XMALLOC (sizeof (char) * (hsize + 1));
   if (header == NULL)
   {
      return;
   }
   string_start = header;

   msprintf (header, "       ");
   header += 7;
   for (i = 0; i < asize; i++)
   {
      msprintf (header, " | %*c", rprec, alphabet_no_2_base(i, sigma));
      header += rprec + 3;
   }
   msprintf (header, "\n");
   header++;
   header[0] = '\0';
   header = string_start;

   /* 6th D: b: "      N | ...\n" */
   store_size += (7 * asize * no);     /* "     N"  */
   store_size += (asize * no);         /* "\n" */
   no *= (asize * asize);
   store_size += (3 * no);             /* " | " */
   store_size += (rprec * no);

   /* alloc memory for the string */
   string = XMALLOC (sizeof (char) * (store_size + 1));
   if (string == NULL)
   {
      return;
   }
   string_start = string;

   /* print to string*/
   for (i = 0; i < scheme->bp_allowed_size; i++)
   {
      msprintf (string, "%c%c:\n",
                alphabet_no_2_base(scheme->bp_allowed[i][0], sigma),
                alphabet_no_2_base(scheme->bp_allowed[i][1], sigma));
      string += 4;
      
      for (j = 0; j < scheme->bp_allowed_size; j++)
      {
         msprintf (string, "  %c%c:\n",
                           alphabet_no_2_base(scheme->bp_allowed[j][0], sigma),
                           alphabet_no_2_base(scheme->bp_allowed[j][1], sigma));
         string += 6;
         
         for (k = 0; k < asize; k++)
         {
            msprintf (string, "    %c:\n", alphabet_no_2_base(k, sigma));
            string += 7;

            for (l = 0; l < asize; l++)
            {
               msprintf (string, "     %c:\n", alphabet_no_2_base(l, sigma));
               string += 8;
               msprintf (string, "%s", header);
               string += hsize;
         
               for (m = 0; m < asize; m++)
               {
                  msprintf (string, "      %c", alphabet_no_2_base(m, sigma));
                  string += 7;
                  
                  for (n = 0; n < asize; n++)
                  {
                     msprintf (string, " | %*i",
                               rprec, scheme->G_int22[i][j][k][l][m][n]);
                     string += 3 + rprec;
                  }
                  
            msprintf (string, "\n");
            string += 1;
               }
            }
         }
      }
   }

/*    /\* print matrix *\/ */
   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   XFREE (string_start);
   XFREE (header);
}

/** @brief Print the mismatch hairpin energies of a scoring scheme to a stream.
 *
 * Prints hairpin loop closing base pair energies to a stream.
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] scheme The scoring scheme.
 * @params[in] sigma The alphabet.
 */
void
nn_scores_fprintf_G_mismatch_interior (FILE* stream,
                                       const NN_scores* scheme,
                                       const Alphabet* sigma)
{
   unsigned long i, j, k;
   int tmp;
   int rprec;
   unsigned long pline_width = 2;
   unsigned long x;             /* 3' base of loop pair */
   unsigned long y;             /* 5' base of loop pair */
   unsigned long z;             /* stacked base pair */
   char* string;
   char* string_start;
   char* header;

   assert (scheme != NULL);
   assert (scheme->G_mismatch_interior != NULL);
   assert (scheme->bp_allowed != NULL);
   assert (sigma != NULL);

   x = y = alphabet_size (sigma);
   z = scheme->bp_allowed_size;

   /* dermine widest cell */
   for (i = 0; i < z; i++)      /* for all base pairs */
   {
      for (j = 0; j < x; j++)   /* for all 3' bases */
      {
         for (k = 0; k < y; k++) /* for all 5' bases */
         {
            rprec = 0;
            tmp = scheme->G_mismatch_interior[i][j][k];

            /* fetch '-' symbol */
            if (tmp < 0)
            {
               tmp *= (-1);
               rprec++;
            }

            /* get no. of digits */
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
   }
   rprec = pline_width;

   /* add up components of a line */
   pline_width += 3;            /*\s|\s*/
   pline_width *= x;
   pline_width += 1;            /* + N */
   pline_width += 1;            /* + \n */

   /* create header line for each table */
   string = (char*) XMALLOC (sizeof (char) * pline_width);
   if (string == NULL)
   {
      return;
   }
   header = string;

   if (rprec > 1)
   {
      rprec--;
   }

   for (i = 0; i < x; i++)
   {
      msprintf (string, "  | %*c", rprec, alphabet_no_2_base (i, sigma));
      string += rprec;
      string += 4;
   }
   string[0] = '\n';
   string++;
   string[0] = '\0';

   if ((rprec + 1) > 1)
   {
      rprec++;
   }

   /* now for the tables */
   /* add to pline_width: room for base pair + table head */
   pline_width *= y;
   pline_width += strlen (header); /* table header */
   pline_width += 4;               /* NN:\n */
   string = (char*) XMALLOC (sizeof (char) * ((pline_width * z) + 1));
   if (string == NULL)
   {
      return;
   }
   string_start = string;  
   
   for (i = 0; i < z; i++)
   {
      /* start new bp table */
      msprintf (string, "%c%c:\n", 
                alphabet_no_2_base (scheme->bp_allowed[i][0], sigma),
                alphabet_no_2_base (scheme->bp_allowed[i][1], sigma));
      string += 4;
      msprintf (string, "%s", header);
      string += strlen (header);

      for (j = 0; j < y; j++)
      {
         msprintf (string, "%c", alphabet_no_2_base (j, sigma));
         string++;
         for (k = 0; k < x; k++)
         {
            msprintf (string, " | %*i", rprec,
                      scheme->G_mismatch_interior[i][j][k]);
            string += rprec;
            string += 3;
         }
         string[0] = '\n';
         string++;
      }
   }

   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   XFREE (string_start);
   XFREE (header);
}

/******************************   Miscellaneous   *****************************/

unsigned long
nn_scores_bp_2_idx (const char base1, const char base2, const NN_scores* scheme)
{
   assert (scheme);

   /*mprintf ("(%c,%c) = %lu\n", alphabet_no_2_base(base1, sigma),
                               alphabet_no_2_base(base2, sigma),
                               base1 * alphabet_size (sigma) + base2);*/

   return scheme->bp_idx[(int)base1][(int)base2];
}

bool
nn_scores_is_allowed_basepair (const char base1, const char base2,
                               void* this)
{
   NN_scores* scheme = (NN_scores*) this;

   assert (scheme);


   if ((unsigned long) scheme->bp_idx[(int)base1][(int)base2]
       < scheme->bp_allowed_size)
   {
      return true;
   }

   return false;
}
