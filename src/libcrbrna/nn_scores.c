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
#include <stdarg.h>
#include <math.h>
#include <libcrbbasic/crbbasic.h>
/*#include "alphabet.h"*/
#include "nn_scores.h"

/* no. of canonical base pairs + whobble GU */
#define NO_ALLOWED_BP 6
/* dimensions of the G_mismatch_hairpin table */
#define D_MM_H 3                /* minimal loop size */
#define D_TL   6                /* size of a tetraloop + closing bp */
#define NN_LXC37 107.856        /* ask marco about this */


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
      this->G_stack                 = NULL;
      this->G_stack_size            = 0;
      this->G_mm_stack              = NULL;
      this->G_mm_stack_size         = 0;
      this->G_hairpin_loop          = NULL;
      this->G_hairpin_loop_size     = 0;
      this->G_mismatch_hairpin      = NULL;
      this->G_mismatch_hairpin_size = 0;
      this->non_gc_penalty_for_bp   = NULL;
      this->tetra_loop              = NULL;
      this->G_tetra_loop            = NULL;
      this->tetra_loop_size         = 0;
      this->bp_idx                  = NULL;
      this->bp_allowed              = NULL;
      this->bp_allowed_size         = 0;
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

void
tetra_loop_swap_entries (unsigned long src, unsigned long dest,
                         unsigned long col,
                         NN_scores* this)
{
   unsigned long i;
   char tmp[D_TL];
   int G_tmp;
   char* ctmp;

   if (this->tetra_loop[dest][col] == this->tetra_loop[src][col])
   {
      return;
   }

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

static void
radix_sort_tetra_loop (unsigned long alpha_size, NN_scores* this, Alphabet* sigma)
{
   unsigned long left_border[UCHAR_MAX + 1];
   unsigned long i;
   char cd;                     /* current digit */
   unsigned long col = D_TL;
   unsigned long lb;

   for (col = D_TL; col > 0; col--)
   {
   memset (left_border, 0, (UCHAR_MAX + 1)*sizeof (unsigned long));

   /* find borders */
   for (i = 0; i < this->tetra_loop_size; i++)
   {
      left_border[this->tetra_loop[i][col - 1] + 1]++;
   }

   mfprintf (stderr, "lb of 0: %lu\n", left_border[0]);
   for (i = 1; i <= alpha_size; i++)
   {
      left_border[i] += left_border[i - 1];
      mfprintf (stderr, "lb of %lu: %lu\n", i, left_border[i]);
   }

   for (i = 0; i < this->tetra_loop_size; i++)
   {
      cd = this->tetra_loop[i][col - 1];

      /* find correct entry for current position */
      if ((int)cd == 0)
         lb = 0;
      else
         lb = left_border[(int)cd - 1];
      mfprintf (stderr, "Working on %2lu (%2i) (%2lu < %2lu < %2lu) ", i, cd,
                lb, i, left_border[(int)cd]);
      while ((lb > i) || (i+1 > left_border[(int)cd]))
      {
         mfprintf (stderr, "SWAP: %lu -> %lu\n", i, left_border[(int)cd]);
         tetra_loop_swap_entries (i, left_border[(int)cd], col - 1, this);
         nn_scores_fprintf_tetra_loop(stderr, this, sigma);
         left_border[(int)cd]++;
         
         cd = this->tetra_loop[i][col - 1];
         if ((int)cd == 0)
            lb = 0;
         else
            lb = left_border[(int)cd - 1];
         mfprintf (stderr, "Working on %2lu (%2i) (%2lu < %2lu < %2lu) ", i, cd,
                   left_border[(int)cd], i, left_border[cd + 1]);
      }
      mfprintf (stderr, "\n");

   }
   }
}

static int
allocate_init_tetra_loop (const int a, const int u, const int g, const int c,
                          unsigned long alpha_size,
                          NN_scores* this,
                          Alphabet* sigma,
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

   radix_sort_tetra_loop (alpha_size, this, sigma);

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

   /* init penalties for non-GC closing base pairs of loops */
   if (allocate_init_non_gc_penalty_for_bp (a, u, g, c, this, file,line))
   {
      nn_scores_delete (this);
      return NULL;      
   }

   /* init tetra loop index and score values */
   if (allocate_init_tetra_loop (a, u, g, c, alphabet_size (sigma), this, sigma,
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
     XFREE (this->non_gc_penalty_for_bp);
     XFREE_2D ((void**)this->tetra_loop);
     XFREE (this->G_tetra_loop);
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
   int bp = (int)this->bp_idx[(int)seq[i]][(int)seq[j]];
   int bip1 = seq[i + 1];
   int bjm1 = seq[j - 1];

   assert (seq);
   assert (this);
   assert (j > 0);

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
      G+= this->G_mismatch_hairpin[bp][bip1][bjm1];
   }
   
   /* tetraloop bonus */
   if (size == 4)
   {
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
      string += 2;      
      for (j = 1; j < (D_TL - 1); j++)
      {
         msprintf (string, "%c", alphabet_no_2_base (scheme->tetra_loop[i][j],
                                                     sigma));
         string ++;
      }
      msprintf (string, "-%c", alphabet_no_2_base (scheme->tetra_loop[i][j],
                                                  sigma));
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
