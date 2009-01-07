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
 *  ToDo:
 *         - change scheme->bp_idx_size to be non quadratic
 *           (get rid of sqrt testing, calc. quadratic size on demand)
 *         - functions: all seq. related input in 5' to 3' order,e.g.get_G_int12
 */


#include <config.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
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
#define NN_LXC37 107.856f       /* ask marco about this */
#define NN_ML_OFFSET 340
#define NN_ML_UNPAIRED 0
#define NN_ML_STEMS 40
#define NN_NINIO_M 50
#define NN_NINIO_MAX 300

struct NN_scores {
      /*c*/float** G_stack;                         /* stacking energies */
      unsigned long G_stack_size;
      /*c*/float** G_mm_stack;                   /* stacks with one mismatch */
      unsigned long G_mm_stack_size;
      /*c*/float* G_hairpin_loop;                   /* hairpin loops */
      unsigned long G_hairpin_loop_size;
      /*c*/float*** G_mismatch_hairpin;        /* hairpin loop closing bp */
      unsigned long G_mismatch_hairpin_size;
      /*c*/float* non_gc_penalty_for_bp;      /* penalty for closing non-GC */
      char** tetra_loop;               /* sorted list of possible loops */
      /*c*/float* G_tetra_loop;                     /* scores */
      unsigned long tetra_loop_size;
      /*c*/float* G_bulge_loop;                     /* bulge loops */
      unsigned long G_bulge_loop_size;
      /* internal loops */
      /*c*/float* G_internal_loop;                  /* generic loops */
      unsigned long G_internal_loop_size;
      /*c*/float**** G_int11;                       /* 1x1 loops */
      unsigned long G_int11_size;
      /*c*/float***** G_int21;                      /* 2x1 loops */
      unsigned long G_int21_size;
      /*c*/float****** G_int22;                      /* 2x2 loops */
      unsigned long G_int22_size;
      /*c*/float*** G_mismatch_interior;         /* interior loop closing bp */
      unsigned long G_mismatch_interior_size;
      /* internal loops */
      /*c*/float** G_dangle5;                  /* 5' dangling end, bp + base */
      unsigned long G_dangle5_size;
      /*c*/float** G_dangle3;                 /* 3' dangling end, bp + base */
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
   NN_scores* this = XOBJ_MALLOC(sizeof (*this), file, line);
   
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
                                               sizeof (**this->bp_allowed),
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
   char i;

   this->bp_idx = (char**) XOBJ_MALLOC_2D (size, size, sizeof (**this->bp_idx),
                                           file, line);
   if (this->bp_idx == NULL)
   {
      return 1;        
   }
   this->bp_idx_size = size * size;

   /* place indeces of allowed base pairs at the begining of the table */
   for (i = 0; (unsigned long) i < this->bp_allowed_size; i++)
   {
      this->bp_idx[(int) this->bp_allowed[(int) i][0]]
                  [(int) this->bp_allowed[(int) i][1]] = i;
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
allocate_init_G_stack (char a, char u, char g, char c,
                       float offset,
                       NN_scores* this,
                       const char* file, const int line)
{
   /*unsigned long i;*/

   this->G_stack_size = this->bp_allowed_size;
   
   /* allocate matrix */
   this->G_stack = (float**) XOBJ_MALLOC_2D (this->G_stack_size,
                                             this->G_stack_size,
                                             sizeof (**this->G_stack),
                                             file, line);
   this->G_stack_size *= this->G_stack_size;
   
   if (this->G_stack == NULL)
   {
      return 1;
   }

   /* BEGIN_STACKING_ENERGIES */
   /* CG CG */
   /* 5'- CG
          GC - 5' */
   this->G_stack[(int) this->bp_idx[(int)c][(int)g]]
                [(int) this->bp_idx[(int)c][(int)g]] = -240 - offset;

   /* CG GC */
   /* 5'- CC
          GG - 5' */
   this->G_stack[(int) this->bp_idx[(int)c][(int)g]]
                [(int) this->bp_idx[(int)g][(int)c]] = -330 - offset;

   /* CG GU */
   /* 5'- CU
          GG - 5' */
   this->G_stack[(int) this->bp_idx[(int)c][(int)g]]
                [(int) this->bp_idx[(int)g][(int)u]] = -210 - offset;

   /* CG UG */
   /* 5'- CG
          GU - 5' */
   this->G_stack[(int) this->bp_idx[(int)c][(int)g]]
                [(int) this->bp_idx[(int)u][(int)g]] = -140 - offset;

   /* CG AU */
   /* 5'- CU
          GA - 5' */
   this->G_stack[(int) this->bp_idx[(int)c][(int)g]]
                [(int) this->bp_idx[(int)a][(int)u]] = -210 - offset;

   /* CG UA */
   /* 5'- CA
          GU - 5' */
   this->G_stack[(int) this->bp_idx[(int)c][(int)g]]
                [(int) this->bp_idx[(int)u][(int)a]] = -210 - offset;

   /* GC CG */
   /* 5'- GG
          CC - 5' */
   this->G_stack[(int) this->bp_idx[(int)g][(int)c]]
                [(int) this->bp_idx[(int)c][(int)g]] = -330 - offset;

   /* GC GC */
   /* 5'- GC
          CG - 5' */
   this->G_stack[(int) this->bp_idx[(int)g][(int)c]]
                [(int) this->bp_idx[(int)g][(int)c]] = -340 - offset;

   /* GC GU */
   /* 5'- GU
          CG - 5' */
   this->G_stack[(int) this->bp_idx[(int)g][(int)c]]
                [(int) this->bp_idx[(int)g][(int)u]] = -250 - offset;

   /* GC UG */
   /* 5'- GG
          CU - 5' */
   this->G_stack[(int) this->bp_idx[(int)g][(int)c]]
                [(int) this->bp_idx[(int)u][(int)g]] = -150 - offset;

   /* GC AU */
   /* 5'- GU
          CA - 5' */
   this->G_stack[(int) this->bp_idx[(int)g][(int)c]]
                [(int) this->bp_idx[(int)a][(int)u]] = -220 - offset;

   /* GC UA */
   /* 5'- GA
          CU - 5' */
   this->G_stack[(int) this->bp_idx[(int)g][(int)c]]
                [(int) this->bp_idx[(int)u][(int)a]] = -240 - offset;

   /* GU CG */
   /* 5'- GG
          UC - 5' */
   this->G_stack[(int) this->bp_idx[(int)g][(int)u]]
                [(int) this->bp_idx[(int)c][(int)g]] = -210 - offset;

   /* GU GC */
   /* 5'- GC
          UG - 5' */
   this->G_stack[(int) this->bp_idx[(int)g][(int)u]]
                [(int) this->bp_idx[(int)g][(int)c]] = -250 - offset;

   /* GU GU */
   /* 5'- GU
          UG - 5' */
   this->G_stack[(int) this->bp_idx[(int)g][(int)u]]
                [(int) this->bp_idx[(int)g][(int)u]] =  130 - offset;

   /* GU UG */
   /* 5'- GG
          UU - 5' */
   this->G_stack[(int) this->bp_idx[(int)g][(int)u]]
                [(int) this->bp_idx[(int)u][(int)g]] =  -50 - offset;

   /* GU AU */
   /* 5'- GU
          UA - 5' */
   this->G_stack[(int) this->bp_idx[(int)g][(int)u]]
                [(int) this->bp_idx[(int)a][(int)u]] = -140 - offset;

   /* GU UA */
   /* 5'- GA
          UU - 5' */
   this->G_stack[(int) this->bp_idx[(int)g][(int)u]]
                [(int) this->bp_idx[(int)u][(int)a]] = -130 - offset;

   /* UG CG */
   /* 5'- UG
          GC - 5' */
   this->G_stack[(int) this->bp_idx[(int)u][(int)g]]
                [(int) this->bp_idx[(int)c][(int)g]] = -140 - offset;

   /* UG GC */
   /* 5'- UC
          GG - 5' */
   this->G_stack[(int) this->bp_idx[(int)u][(int)g]]
                [(int) this->bp_idx[(int)g][(int)c]] = -150 - offset;

   /* UG GU */
   /* 5'- UU
          GG - 5' */
   this->G_stack[(int) this->bp_idx[(int)u][(int)g]]
                [(int) this->bp_idx[(int)g][(int)u]] =  -50 - offset;

   /* UG UG */
   /* 5'- UG
          GU - 5' */
   this->G_stack[(int) this->bp_idx[(int)u][(int)g]]
                [(int) this->bp_idx[(int)u][(int)g]] =   30 - offset;

   /* UG AU */
   /* 5'- UU
          GA - 5' */
   this->G_stack[(int) this->bp_idx[(int)u][(int)g]]
                [(int) this->bp_idx[(int)a][(int)u]] =  -60 - offset;

   /* UG UA */
   /* 5'- UA
          GU - 5' */
   this->G_stack[(int) this->bp_idx[(int)u][(int)g]]
                [(int) this->bp_idx[(int)u][(int)a]] = -100 - offset;

   /* AU CG */
   /* 5'- AG
          UC - 5' */
   this->G_stack[(int) this->bp_idx[(int)a][(int)u]]
                [(int) this->bp_idx[(int)c][(int)g]] = -210 - offset;

   /* AU GC */
   /* 5'- AC
          UG - 5' */
   this->G_stack[(int) this->bp_idx[(int)a][(int)u]]
                [(int) this->bp_idx[(int)g][(int)c]] = -220 - offset;

   /* AU GU */
   /* 5'- AU
          UG - 5' */
   this->G_stack[(int) this->bp_idx[(int)a][(int)u]]
                [(int) this->bp_idx[(int)g][(int)u]] = -140 - offset;

   /* AU UG */
   /* 5'- AG
          UU - 5' */
   this->G_stack[(int) this->bp_idx[(int)a][(int)u]]
                [(int) this->bp_idx[(int)u][(int)g]] =  -60 - offset;

   /* AU AU */
   /* 5'- AU
          UA - 5' */
   this->G_stack[(int) this->bp_idx[(int)a][(int)u]]
                [(int) this->bp_idx[(int)a][(int)u]] = -110 - offset;

   /* AU UA */
   /* 5'- AA
          UU - 5' */
   this->G_stack[(int) this->bp_idx[(int)a][(int)u]]
                [(int) this->bp_idx[(int)u][(int)a]] =  -90 - offset;

   /* UA CG */
   /* 5'- UG
          AC - 5' */
   this->G_stack[(int) this->bp_idx[(int)u][(int)a]]
                [(int) this->bp_idx[(int)c][(int)g]] = -210 - offset;

   /* UA GC */
   /* 5'- UC
          AG - 5' */
   this->G_stack[(int) this->bp_idx[(int)u][(int)a]]
                [(int) this->bp_idx[(int)g][(int)c]] = -240 - offset;

   /* UA GU */
   /* 5'- UU
          AG - 5' */
   this->G_stack[(int) this->bp_idx[(int)u][(int)a]]
                [(int) this->bp_idx[(int)g][(int)u]] = -130 - offset;

   /* UA UG */
   /* 5'- UG
          AU - 5' */
   this->G_stack[(int) this->bp_idx[(int)u][(int)a]]
                [(int) this->bp_idx[(int)u][(int)g]] = -100 - offset;

   /* UA AU */
   /* 5'- UU
          AA - 5' */
   this->G_stack[(int) this->bp_idx[(int)u][(int)a]]
                [(int) this->bp_idx[(int)a][(int)u]] =  -90 - offset;

   /* UA UA */
   /* 5'- UA
          AU - 5' */
   this->G_stack[(int) this->bp_idx[(int)u][(int)a]]
                [(int) this->bp_idx[(int)u][(int)a]] = -130 - offset;

   /* END_STACKING_ENERGIES */

   return 0;
}

static int
allocate_init_G_mm_stack_size (int a, int u, int g, int c,
                               unsigned long size,
                               float offset,
                               NN_scores* this,
                               const char* file, const int line)
{
   unsigned long bp1;

   this->G_mm_stack_size = size * size;
   
   this->G_mm_stack = (float**) XOBJ_MALLOC_2D (this->bp_allowed_size,
                                                this->G_mm_stack_size,
                                                sizeof (**this->G_mm_stack),
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

   /* BEGIN_MISMATCH_STACK */
   /* stacks containing a mismatch */
   /* mi: param from mismatch_interior table */
   /* mh: param from mismatch_hairpin */

   bp1 = this->bp_idx[c][g];
   /* CG AA */
   /* mh: -150 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][a]] =  -75 - offset;
   /* CG AC */
   /* mh: -150 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][c]] =  -75 - offset;
   /* CG AG */
   /* mh: -140 mi: -110*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][g]] = -125 - offset;
   /* CG AU */
   /* mh: -180 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][u]] =  -90 - offset;

   /* CG CA */
   /* mh: -100 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][a]] =  -50 - offset;
   /* CG CC */
   /* mh: -90 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][c]] =  -45 - offset;
   /* CG CG */
   /* mh: -290 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][g]] = -145 - offset;
   /* CG CU */
   /* mh: -80 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][u]] =  -40 - offset;

   /* CG GA */
   /* mh: -220 mi: -110*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][a]] = -165 - offset;
   /* CG GC */
   /* mh: -200 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][c]] = -100 - offset;
   /* CG GG */
   /* mh: -160 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][g]] =  -80 - offset;
   /* CG GU */
   /* mh: -110 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][u]] =  -55 - offset;

   /* CG UA */
   /* mh: -170 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][a]] =  -85 - offset;
   /* CG UC */
   /* mh: -140 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][c]] =  -70 - offset;
   /* CG UG */
   /* mh: -180 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][g]] =  -90 - offset;
   /* CG UU */
   /* mh: -200 mi: -70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][u]] = -135 - offset;

   bp1 = this->bp_idx[g][c];
   /* GC AA */
   /* mh: -110 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][a]] =  -55 - offset;
   /* GC AC */
   /* mh: -150 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][c]] =  -75 - offset;
   /* GC AG */
   /* mh: -130 mi: -110*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][g]] = -120 - offset;
   /* GC AU */
   /* mh: -210 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][u]] = -105 - offset;

   /* GC CA */
   /* mh: -110 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][a]] =  -55 - offset;
   /* GC CC */
   /* mh: -70 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][c]] =  -35 - offset;
   /* GC CG */
   /* mh: -240 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][g]] = -120 - offset;
   /* GC CU */
   /* mh: -50 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][u]] =  -25 - offset;

   /* GC GA */
   /* mh: -240 mi: -110*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][a]] = -175 - offset;
   /* GC GC */
   /* mh: -290 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][c]] = -145 - offset;
   /* GC GG */
   /* mh: -140 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][g]] =  -70 - offset;
   /* GC GU */
   /* mh: -120 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][u]] =  -60 - offset;

   /* GC UA */
   /* mh: -190 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][a]] =  -95 - offset;
   /* GC UC */
   /* mh: -100 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][c]] =  -50 - offset;
   /* GC UG */
   /* mh: -220 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][g]] = -110 - offset;
   /* GC UU */
   /* mh: -150 mi: -70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][u]] = -110 - offset;

   bp1 = this->bp_idx[g][u];
   /* GU AA */
   /* mh: 20 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][a]] =   45 - offset;
   /* GU AC */
   /* mh: -50 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][c]] =   10 - offset;
   /* GU AG */
   /* mh: -30 mi: -40*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][g]] =  -35 - offset;
   /* GU AU */
   /* mh: -30 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][u]] =   20 - offset;

   /* GU CA */
   /* mh: -10 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][a]] =   30 - offset;
   /* GU CC */
   /* mh: -20 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][c]] =   25 - offset;
   /* GU CG */
   /* mh: -150 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][g]] =  -40 - offset;
   /* GU CU */
   /* mh: -20 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][u]] =   25 - offset;

   /* GU GA */
   /* mh: -90 mi: -40*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][a]] =  -65 - offset;
   /* GU GC */
   /* mh: -110 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][c]] =  -20 - offset;
   /* GU GG */
   /* mh: -30 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][g]] =   20 - offset;
   /* GU GU */
   /* mh: 0 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][u]] =   35 - offset;

   /* GU UA */
   /* mh: -30 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][a]] =   20 - offset;
   /* GU UC */
   /* mh: -30 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][c]] =   20 - offset;
   /* GU UG */
   /* mh: -40 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][g]] =   15 - offset;
   /* GU UU */
   /* mh: -110 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][u]] =  -55 - offset;

   bp1 = this->bp_idx[u][g];
   /* UG AA */
   /* mh: -50 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][a]] =   10 - offset;
   /* UG AC */
   /* mh: -30 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][c]] =   20 - offset;
   /* UG AG */
   /* mh: -60 mi: -40*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][g]] =  -50 - offset;
   /* UG AU */
   /* mh: -50 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][u]] =   10 - offset;

   /* UG CA */
   /* mh: -20 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][a]] =   25 - offset;
   /* UG CC */
   /* mh: -10 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][c]] =   30 - offset;
   /* UG CG */
   /* mh: -170 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][g]] =  -50 - offset;
   /* UG CU */
   /* mh: 0 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][u]] =   35 - offset;

   /* UG GA */
   /* mh: -80 mi: -40*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][a]] =  -60 - offset;
   /* UG GC */
   /* mh: -120 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][c]] =  -25 - offset;
   /* UG GG */
   /* mh: -30 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][g]] =   20 - offset;
   /* UG GU */
   /* mh: -70 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][u]] =    0 - offset;

   /* UG UA */
   /* mh: -60 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][a]] =    5 - offset;
   /* UG UC */
   /* mh: -10 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][c]] =   30 - offset;
   /* UG UG */
   /* mh: -60 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][g]] =    5 - offset;
   /* UG UU */
   /* mh: -80 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][u]] =  -40 - offset;

   bp1 = this->bp_idx[a][u];
   /* AU AA */
   /* mh: -30 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][a]] =   20 - offset;
   /* AU AC */
   /* mh: -50 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][c]] =   10 - offset;
   /* AU AG */
   /* mh: -30 mi: -40*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][g]] =  -35 - offset;
   /* AU AU */
   /* mh: -30 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][u]] =   20 - offset;

   /* AU CA */
   /* mh: -10 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][a]] =   30 - offset;
   /* AU CC */
   /* mh: -20 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][c]] =   25 - offset;
   /* AU CG */
   /* mh: -150 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][g]] =  -40 - offset;
   /* AU CU */
   /* mh: -20 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][u]] =   25 - offset;

   /* AU GA */
   /* mh: -110 mi: -40*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][a]] =  -75 - offset;
   /* AU GC */
   /* mh: -120 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][c]] =  -25 - offset;
   /* AU GG */
   /* mh: -20 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][g]] =   25 - offset;
   /* AU GU */
   /* mh: 20 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][u]] =   45 - offset;

   /* AU UA */
   /* mh: -30 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][a]] =   20 - offset;
   /* AU UC */
   /* mh: -30 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][c]] =   20 - offset;
   /* AU UG */
   /* mh: -60 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][g]] =    5 - offset;
   /* AU UU */
   /* mh: -110 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][u]] =  -55 - offset;

   bp1 = this->bp_idx[u][a];
   /* UA AA */
   /* mh: -50 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][a]] =   10 - offset;
   /* UA AC */
   /* mh: -30 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][c]] =   20 - offset;
   /* UA AG */
   /* mh: -60 mi: -40*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][g]] =  -50 - offset;
   /* UA AU */
   /* mh: -50 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[a][u]] =   10 - offset;

   /* UA CA */
   /* mh: -20 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][a]] =   25 - offset;
   /* UA CC */
   /* mh: -10 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][c]] =   30 - offset;
   /* UA CG */
   /* mh: -120 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][g]] =  -25 - offset;
   /* UA CU */
   /* mh: 0 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[c][u]] =   35 - offset;

   /* UA GA */
   /* mh: -140 mi: -40*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][a]] =  -90 - offset;
   /* UA GC */
   /* mh: -120 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][c]] =  -25 - offset;
   /* UA GG */
   /* mh: -70 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][g]] =    0 - offset;
   /* UA GU */
   /* mh: -20 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[g][u]] =   25 - offset;

   /* UA UA */
   /* mh: -30 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][a]] =   20 - offset;
   /* UA UC */
   /* mh: -10 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][c]] =   30 - offset;
   /* UA UG */
   /* mh: -50 mi: 70*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][g]] =   10 - offset;
   /* UA UU */
   /* mh: -80 mi: 0*/
   this->G_mm_stack[bp1][(int) this->bp_idx[u][u]] =  -40 - offset;

   /* END_MISMATCH_STACK */

   return 0;
}

static int
allocate_init_hairpin_loop (float offset, NN_scores* this,
                            const char* file, const int line)
{

   this->G_hairpin_loop_size = 31;
   this->G_hairpin_loop = (float*) XOBJ_MALLOC (  this->G_hairpin_loop_size
                                            * sizeof (this->G_hairpin_loop[0]),
                                                  file, line);
   if (this->G_hairpin_loop == NULL)
   {
      return 1;        
   }

   /* BEGIN_HAIRPINS */
   this->G_hairpin_loop[ 0] = FLOAT_UNDEF;
   this->G_hairpin_loop[ 1] = FLOAT_UNDEF;
   this->G_hairpin_loop[ 2] = FLOAT_UNDEF;
   this->G_hairpin_loop[ 3] = 570 - offset;
   this->G_hairpin_loop[ 4] = 560 - offset;
   this->G_hairpin_loop[ 5] = 560 - offset;
   this->G_hairpin_loop[ 6] = 540 - offset;
   this->G_hairpin_loop[ 7] = 590 - offset;
   this->G_hairpin_loop[ 8] = 560 - offset;
   this->G_hairpin_loop[ 9] = 640 - offset;
   this->G_hairpin_loop[10] = 650 - offset;
   this->G_hairpin_loop[11] = 660 - offset;
   this->G_hairpin_loop[12] = 670 - offset;
   this->G_hairpin_loop[13] = 678 - offset;
   this->G_hairpin_loop[14] = 686 - offset;
   this->G_hairpin_loop[15] = 694 - offset;
   this->G_hairpin_loop[16] = 701 - offset;
   this->G_hairpin_loop[17] = 707 - offset;
   this->G_hairpin_loop[18] = 713 - offset;
   this->G_hairpin_loop[19] = 719 - offset;
   this->G_hairpin_loop[20] = 725 - offset;
   this->G_hairpin_loop[21] = 730 - offset;
   this->G_hairpin_loop[22] = 735 - offset;
   this->G_hairpin_loop[23] = 740 - offset;
   this->G_hairpin_loop[24] = 744 - offset;
   this->G_hairpin_loop[25] = 749 - offset;
   this->G_hairpin_loop[26] = 753 - offset;
   this->G_hairpin_loop[27] = 757 - offset;
   this->G_hairpin_loop[28] = 761 - offset;
   this->G_hairpin_loop[29] = 765 - offset;
   this->G_hairpin_loop[30] = 769 - offset;
   /* END_HAIRPINS */

   return 0;
}

static int
allocate_init_mismatch_hairpin (int a, int u, int g, int c,
                                const unsigned long no_of_b,
                                float offset,
                                NN_scores* this,
                                const char* file, const int line)
{   
   unsigned long bp1;

   /* allocate memory */
   this->G_mismatch_hairpin
      = (float***) XOBJ_MALLOC_ND(sizeof (***this->G_mismatch_hairpin),
                                  D_MM_H,
                                  file, line,
                                  this->bp_allowed_size, no_of_b, no_of_b);
   if (this->G_mismatch_hairpin == NULL)
   {
      return 1;        
   }
   this->G_mismatch_hairpin_size = this->bp_allowed_size * no_of_b * no_of_b;

   /* BEGIN_MISMATCH_HAIRPIN */
   /* CG */
   bp1 = this->bp_idx[c][g];
   this->G_mismatch_hairpin[bp1][a][a] = -150 - offset; /* AA */
   this->G_mismatch_hairpin[bp1][a][c] = -150 - offset; /* AC */
   this->G_mismatch_hairpin[bp1][a][g] = -140 - offset; /* AG */
   this->G_mismatch_hairpin[bp1][a][u] = -180 - offset; /* AU */

   this->G_mismatch_hairpin[bp1][c][a] = -100 - offset; /* CA */
   this->G_mismatch_hairpin[bp1][c][c] =  -90 - offset; /* CC */
   this->G_mismatch_hairpin[bp1][c][g] = -290 - offset; /* CG */
   this->G_mismatch_hairpin[bp1][c][u] =  -80 - offset; /* CU */

   this->G_mismatch_hairpin[bp1][g][a] = -220 - offset; /* GA */
   this->G_mismatch_hairpin[bp1][g][c] = -200 - offset; /* GC */
   this->G_mismatch_hairpin[bp1][g][g] = -160 - offset; /* GG */
   this->G_mismatch_hairpin[bp1][g][u] = -110 - offset; /* GU */

   this->G_mismatch_hairpin[bp1][u][a] = -170 - offset; /* UA */
   this->G_mismatch_hairpin[bp1][u][c] = -140 - offset; /* UC */
   this->G_mismatch_hairpin[bp1][u][g] = -180 - offset; /* UG */
   this->G_mismatch_hairpin[bp1][u][u] = -200 - offset; /* UU */

   /* GC */
   bp1 = this->bp_idx[g][c];
   this->G_mismatch_hairpin[bp1][a][a] = -110 - offset; /* AA */
   this->G_mismatch_hairpin[bp1][a][c] = -150 - offset; /* AC */
   this->G_mismatch_hairpin[bp1][a][g] = -130 - offset; /* AG */
   this->G_mismatch_hairpin[bp1][a][u] = -210 - offset; /* AU */

   this->G_mismatch_hairpin[bp1][c][a] = -110 - offset; /* CA */
   this->G_mismatch_hairpin[bp1][c][c] =  -70 - offset; /* CC */
   this->G_mismatch_hairpin[bp1][c][g] = -240 - offset; /* CG */
   this->G_mismatch_hairpin[bp1][c][u] =  -50 - offset; /* CU */

   this->G_mismatch_hairpin[bp1][g][a] = -240 - offset; /* GA */
   this->G_mismatch_hairpin[bp1][g][c] = -290 - offset; /* GC */
   this->G_mismatch_hairpin[bp1][g][g] = -140 - offset; /* GG */
   this->G_mismatch_hairpin[bp1][g][u] = -120 - offset; /* GU */

   this->G_mismatch_hairpin[bp1][u][a] = -190 - offset; /* UA */
   this->G_mismatch_hairpin[bp1][u][c] = -100 - offset; /* UC */
   this->G_mismatch_hairpin[bp1][u][g] = -220 - offset; /* UG */
   this->G_mismatch_hairpin[bp1][u][u] = -150 - offset; /* UU */

   /* GU */
   bp1 = this->bp_idx[g][u];
   this->G_mismatch_hairpin[bp1][a][a] =   20 - offset; /* AA */
   this->G_mismatch_hairpin[bp1][a][c] =  -50 - offset; /* AC */
   this->G_mismatch_hairpin[bp1][a][g] =  -30 - offset; /* AG */
   this->G_mismatch_hairpin[bp1][a][u] =  -30 - offset; /* AU */

   this->G_mismatch_hairpin[bp1][c][a] =  -10 - offset; /* CA */
   this->G_mismatch_hairpin[bp1][c][c] =  -20 - offset; /* CC */
   this->G_mismatch_hairpin[bp1][c][g] = -150 - offset; /* CG */
   this->G_mismatch_hairpin[bp1][c][u] =  -20 - offset; /* CU */

   this->G_mismatch_hairpin[bp1][g][a] =  -90 - offset; /* GA */
   this->G_mismatch_hairpin[bp1][g][c] = -110 - offset; /* GC */
   this->G_mismatch_hairpin[bp1][g][g] =  -30 - offset; /* GG */
   this->G_mismatch_hairpin[bp1][g][u] =    0 - offset; /* GU */

   this->G_mismatch_hairpin[bp1][u][a] =  -30 - offset; /* UA */
   this->G_mismatch_hairpin[bp1][u][c] =  -30 - offset; /* UC */
   this->G_mismatch_hairpin[bp1][u][g] =  -40 - offset; /* UG */
   this->G_mismatch_hairpin[bp1][u][u] = -110 - offset; /* UU */

   /* UG */
   bp1 = this->bp_idx[u][g];
   this->G_mismatch_hairpin[bp1][a][a] =  -50 - offset; /* AA */
   this->G_mismatch_hairpin[bp1][a][c] =  -30 - offset; /* AC */
   this->G_mismatch_hairpin[bp1][a][g] =  -60 - offset; /* AG */
   this->G_mismatch_hairpin[bp1][a][u] =  -50 - offset; /* AU */

   this->G_mismatch_hairpin[bp1][c][a] =  -20 - offset; /* CA */
   this->G_mismatch_hairpin[bp1][c][c] =  -10 - offset; /* CC */
   this->G_mismatch_hairpin[bp1][c][g] = -170 - offset; /* CG */
   this->G_mismatch_hairpin[bp1][c][u] =    0 - offset; /* CU */

   this->G_mismatch_hairpin[bp1][g][a] =  -80 - offset; /* GA */
   this->G_mismatch_hairpin[bp1][g][c] = -120 - offset; /* GC */
   this->G_mismatch_hairpin[bp1][g][g] =  -30 - offset; /* GG */
   this->G_mismatch_hairpin[bp1][g][u] =  -70 - offset; /* GU */

   this->G_mismatch_hairpin[bp1][u][a] =  -60 - offset; /* UA */
   this->G_mismatch_hairpin[bp1][u][c] =  -10 - offset; /* UC */
   this->G_mismatch_hairpin[bp1][u][g] =  -60 - offset; /* UG */
   this->G_mismatch_hairpin[bp1][u][u] =  -80 - offset; /* UU */

   /* AU */
   bp1 = this->bp_idx[a][u];
   this->G_mismatch_hairpin[bp1][a][a] =  -30 - offset; /* AA */
   this->G_mismatch_hairpin[bp1][a][c] =  -50 - offset; /* AC */
   this->G_mismatch_hairpin[bp1][a][g] =  -30 - offset; /* AG */
   this->G_mismatch_hairpin[bp1][a][u] =  -30 - offset; /* AU */

   this->G_mismatch_hairpin[bp1][c][a] =  -10 - offset; /* CA */
   this->G_mismatch_hairpin[bp1][c][c] =  -20 - offset; /* CC */
   this->G_mismatch_hairpin[bp1][c][g] = -150 - offset; /* CG */
   this->G_mismatch_hairpin[bp1][c][u] =  -20 - offset; /* CU */

   this->G_mismatch_hairpin[bp1][g][a] = -110 - offset; /* GA */
   this->G_mismatch_hairpin[bp1][g][c] = -120 - offset; /* GC */
   this->G_mismatch_hairpin[bp1][g][g] =  -20 - offset; /* GG */
   this->G_mismatch_hairpin[bp1][g][u] =   20 - offset; /* GU */

   this->G_mismatch_hairpin[bp1][u][a] =  -30 - offset; /* UA */
   this->G_mismatch_hairpin[bp1][u][c] =  -30 - offset; /* UC */
   this->G_mismatch_hairpin[bp1][u][g] =  -60 - offset; /* UG */
   this->G_mismatch_hairpin[bp1][u][u] = -110 - offset; /* UU */

   /* UA */
   bp1 = this->bp_idx[u][a];
   this->G_mismatch_hairpin[bp1][a][a] =  -50 - offset; /* AA */
   this->G_mismatch_hairpin[bp1][a][c] =  -30 - offset; /* AC */
   this->G_mismatch_hairpin[bp1][a][g] =  -60 - offset; /* AG */
   this->G_mismatch_hairpin[bp1][a][u] =  -50 - offset; /* AU */

   this->G_mismatch_hairpin[bp1][c][a] =  -20 - offset; /* CA */
   this->G_mismatch_hairpin[bp1][c][c] =  -10 - offset; /* CC */
   this->G_mismatch_hairpin[bp1][c][g] = -120 - offset; /* CG */
   this->G_mismatch_hairpin[bp1][c][u] =    0 - offset; /* CU */

   this->G_mismatch_hairpin[bp1][g][a] = -140 - offset; /* GA */
   this->G_mismatch_hairpin[bp1][g][c] = -120 - offset; /* GC */
   this->G_mismatch_hairpin[bp1][g][g] =  -70 - offset; /* GG */
   this->G_mismatch_hairpin[bp1][g][u] =  -20 - offset; /* GU */

   this->G_mismatch_hairpin[bp1][u][a] =  -30 - offset; /* UA */
   this->G_mismatch_hairpin[bp1][u][c] =  -10 - offset; /* UC */
   this->G_mismatch_hairpin[bp1][u][g] =  -50 - offset; /* UG */
   this->G_mismatch_hairpin[bp1][u][u] =  -80 - offset; /* UU */

   /* END_MISMATCH_HAIRPIN */

   return 0;
}

static int
allocate_init_mismatch_interior (int a, int u, int g, int c,
                                 const unsigned long no_of_b,
                                 float offset,
                                 NN_scores* this,
                                 const char* file, const int line)
{
   unsigned long bp1;

   /* allocate memory */
   this->G_mismatch_interior
      = (float***) XOBJ_MALLOC_ND(sizeof (***this->G_mismatch_interior),
                                  D_MM_I,
                                  file, line,
                                  this->bp_allowed_size, no_of_b, no_of_b);
   if (this->G_mismatch_interior == NULL)
   {
      return 1;        
   }
   this->G_mismatch_interior_size = this->bp_allowed_size * no_of_b * no_of_b;

   /* BEGIN_MISMATCH_INTERIOR */
   /* CG */
   bp1 = this->bp_idx[c][g];
   this->G_mismatch_interior[bp1][a][a] =    0 - offset; /* AA */
   this->G_mismatch_interior[bp1][a][c] =    0 - offset; /* AC */
   this->G_mismatch_interior[bp1][a][g] = -110 - offset; /* AG */
   this->G_mismatch_interior[bp1][a][u] =    0 - offset; /* AU */

   this->G_mismatch_interior[bp1][c][a] =    0 - offset; /* CA */
   this->G_mismatch_interior[bp1][c][c] =    0 - offset; /* CC */
   this->G_mismatch_interior[bp1][c][g] =    0 - offset; /* CG */
   this->G_mismatch_interior[bp1][c][u] =    0 - offset; /* CU */

   this->G_mismatch_interior[bp1][g][a] = -110 - offset; /* GA */
   this->G_mismatch_interior[bp1][g][c] =    0 - offset; /* GC */
   this->G_mismatch_interior[bp1][g][g] =    0 - offset; /* GG */
   this->G_mismatch_interior[bp1][g][u] =    0 - offset; /* GU */

   this->G_mismatch_interior[bp1][u][a] =    0 - offset; /* UA */
   this->G_mismatch_interior[bp1][u][c] =    0 - offset; /* UC */
   this->G_mismatch_interior[bp1][u][g] =    0 - offset; /* UG */
   this->G_mismatch_interior[bp1][u][u] =  -70 - offset; /* UU */

   /* GC */
   bp1 = this->bp_idx[g][c];
   this->G_mismatch_interior[bp1][a][a] =    0 - offset; /* AA */
   this->G_mismatch_interior[bp1][a][c] =    0 - offset; /* AC */
   this->G_mismatch_interior[bp1][a][g] = -110 - offset; /* AG */
   this->G_mismatch_interior[bp1][a][u] =    0 - offset; /* AU */

   this->G_mismatch_interior[bp1][c][a] =    0 - offset; /* CA */
   this->G_mismatch_interior[bp1][c][c] =    0 - offset; /* CC */
   this->G_mismatch_interior[bp1][c][g] =    0 - offset; /* CG */
   this->G_mismatch_interior[bp1][c][u] =    0 - offset; /* CU */

   this->G_mismatch_interior[bp1][g][a] = -110 - offset; /* GA */
   this->G_mismatch_interior[bp1][g][c] =    0 - offset; /* GC */
   this->G_mismatch_interior[bp1][g][g] =    0 - offset; /* GG */
   this->G_mismatch_interior[bp1][g][u] =    0 - offset; /* GU */

   this->G_mismatch_interior[bp1][u][a] =    0 - offset; /* UA */
   this->G_mismatch_interior[bp1][u][c] =    0 - offset; /* UC */
   this->G_mismatch_interior[bp1][u][g] =    0 - offset; /* UG */
   this->G_mismatch_interior[bp1][u][u] =  -70 - offset; /* UU */

   /* GU */
   bp1 = this->bp_idx[g][u];
   this->G_mismatch_interior[bp1][a][a] =   70 - offset; /* AA */
   this->G_mismatch_interior[bp1][a][c] =   70 - offset; /* AC */
   this->G_mismatch_interior[bp1][a][g] =  -40 - offset; /* AG */
   this->G_mismatch_interior[bp1][a][u] =   70 - offset; /* AU */

   this->G_mismatch_interior[bp1][c][a] =   70 - offset; /* CA */
   this->G_mismatch_interior[bp1][c][c] =   70 - offset; /* CC */
   this->G_mismatch_interior[bp1][c][g] =   70 - offset; /* CG */
   this->G_mismatch_interior[bp1][c][u] =   70 - offset; /* CU */

   this->G_mismatch_interior[bp1][g][a] =  -40 - offset; /* GA */
   this->G_mismatch_interior[bp1][g][c] =   70 - offset; /* GC */
   this->G_mismatch_interior[bp1][g][g] =   70 - offset; /* GG */
   this->G_mismatch_interior[bp1][g][u] =   70 - offset; /* GU */

   this->G_mismatch_interior[bp1][u][a] =   70 - offset; /* UA */
   this->G_mismatch_interior[bp1][u][c] =   70 - offset; /* UC */
   this->G_mismatch_interior[bp1][u][g] =   70 - offset; /* UG */
   this->G_mismatch_interior[bp1][u][u] =    0 - offset; /* UU */

   /* UG */
   bp1 = this->bp_idx[u][g];
   this->G_mismatch_interior[bp1][a][a] =   70 - offset; /* AA */
   this->G_mismatch_interior[bp1][a][c] =   70 - offset; /* AC */
   this->G_mismatch_interior[bp1][a][g] =  -40 - offset; /* AG */
   this->G_mismatch_interior[bp1][a][u] =   70 - offset; /* AU */

   this->G_mismatch_interior[bp1][c][a] =   70 - offset; /* CA */
   this->G_mismatch_interior[bp1][c][c] =   70 - offset; /* CC */
   this->G_mismatch_interior[bp1][c][g] =   70 - offset; /* CG */
   this->G_mismatch_interior[bp1][c][u] =   70 - offset; /* CU */

   this->G_mismatch_interior[bp1][g][a] =  -40 - offset; /* GA */
   this->G_mismatch_interior[bp1][g][c] =   70 - offset; /* GC */
   this->G_mismatch_interior[bp1][g][g] =   70 - offset; /* GG */
   this->G_mismatch_interior[bp1][g][u] =   70 - offset; /* GU */

   this->G_mismatch_interior[bp1][u][a] =   70 - offset; /* UA */
   this->G_mismatch_interior[bp1][u][c] =   70 - offset; /* UC */
   this->G_mismatch_interior[bp1][u][g] =   70 - offset; /* UG */
   this->G_mismatch_interior[bp1][u][u] =    0 - offset; /* UU */

   /* AU */
   bp1 = this->bp_idx[a][u];
   this->G_mismatch_interior[bp1][a][a] =   70 - offset; /* AA */
   this->G_mismatch_interior[bp1][a][c] =   70 - offset; /* AC */
   this->G_mismatch_interior[bp1][a][g] =  -40 - offset; /* AG */
   this->G_mismatch_interior[bp1][a][u] =   70 - offset; /* AU */

   this->G_mismatch_interior[bp1][c][a] =   70 - offset; /* CA */
   this->G_mismatch_interior[bp1][c][c] =   70 - offset; /* CC */
   this->G_mismatch_interior[bp1][c][g] =   70 - offset; /* CG */
   this->G_mismatch_interior[bp1][c][u] =   70 - offset; /* CU */

   this->G_mismatch_interior[bp1][g][a] =  -40 - offset; /* GA */
   this->G_mismatch_interior[bp1][g][c] =   70 - offset; /* GC */
   this->G_mismatch_interior[bp1][g][g] =   70 - offset; /* GG */
   this->G_mismatch_interior[bp1][g][u] =   70 - offset; /* GU */

   this->G_mismatch_interior[bp1][u][a] =   70 - offset; /* UA */
   this->G_mismatch_interior[bp1][u][c] =   70 - offset; /* UC */
   this->G_mismatch_interior[bp1][u][g] =   70 - offset; /* UG */
   this->G_mismatch_interior[bp1][u][u] =    0 - offset; /* UU */

   /* UA */
   bp1 = this->bp_idx[u][a];
   this->G_mismatch_interior[bp1][a][a] =   70 - offset; /* AA */
   this->G_mismatch_interior[bp1][a][c] =   70 - offset; /* AC */
   this->G_mismatch_interior[bp1][a][g] =  -40 - offset; /* AG */
   this->G_mismatch_interior[bp1][a][u] =   70 - offset; /* AU */

   this->G_mismatch_interior[bp1][c][a] =   70 - offset; /* CA */
   this->G_mismatch_interior[bp1][c][c] =   70 - offset; /* CC */
   this->G_mismatch_interior[bp1][c][g] =   70 - offset; /* CG */
   this->G_mismatch_interior[bp1][c][u] =   70 - offset; /* CU */

   this->G_mismatch_interior[bp1][g][a] =  -40 - offset; /* GA */
   this->G_mismatch_interior[bp1][g][c] =   70 - offset; /* GC */
   this->G_mismatch_interior[bp1][g][g] =   70 - offset; /* GG */
   this->G_mismatch_interior[bp1][g][u] =   70 - offset; /* GU */

   this->G_mismatch_interior[bp1][u][a] =   70 - offset; /* UA */
   this->G_mismatch_interior[bp1][u][c] =   70 - offset; /* UC */
   this->G_mismatch_interior[bp1][u][g] =   70 - offset; /* UG */
   this->G_mismatch_interior[bp1][u][u] =    0 - offset; /* UU */

   /* END_MISMATCH_INTERIOR */

   return 0;
}

static int
allocate_init_internal_loop (float offset, NN_scores* this,
                             const char* file, const int line)
{

   this->G_internal_loop_size = 31;
   this->G_internal_loop = XOBJ_MALLOC (  this->G_internal_loop_size
                                          * sizeof (this->G_internal_loop[0]),
                                          file, line);
   if (this->G_internal_loop == NULL)
   {
      return 1;
   }

   /* BEGIN_INTERNALS */
   this->G_internal_loop[ 0] = FLOAT_UNDEF;
   this->G_internal_loop[ 1] = FLOAT_UNDEF;
   this->G_internal_loop[ 2] = 410 - offset;
   this->G_internal_loop[ 3] = 510 - offset;
   this->G_internal_loop[ 4] = 170 - offset;
   this->G_internal_loop[ 5] = 180 - offset;
   this->G_internal_loop[ 6] = 200 - offset;
   this->G_internal_loop[ 7] = 220 - offset;
   this->G_internal_loop[ 8] = 230 - offset;
   this->G_internal_loop[ 9] = 240 - offset;
   this->G_internal_loop[10] = 250 - offset;
   this->G_internal_loop[11] = 260 - offset;
   this->G_internal_loop[12] = 270 - offset;
   this->G_internal_loop[13] = 278 - offset;
   this->G_internal_loop[14] = 286 - offset;
   this->G_internal_loop[15] = 294 - offset;
   this->G_internal_loop[16] = 301 - offset;
   this->G_internal_loop[17] = 307 - offset;
   this->G_internal_loop[18] = 313 - offset;
   this->G_internal_loop[19] = 319 - offset;
   this->G_internal_loop[20] = 325 - offset;
   this->G_internal_loop[21] = 330 - offset;
   this->G_internal_loop[22] = 335 - offset;
   this->G_internal_loop[23] = 340 - offset;
   this->G_internal_loop[24] = 345 - offset;
   this->G_internal_loop[25] = 349 - offset;
   this->G_internal_loop[26] = 353 - offset;
   this->G_internal_loop[27] = 357 - offset;
   this->G_internal_loop[28] = 361 - offset;
   this->G_internal_loop[29] = 365 - offset;
   this->G_internal_loop[30] = 369 - offset;
   /* END_INTERNALS */

   return 0;
}

static int
allocate_init_int11 (const int a, const int u, const int g, const int c,
                     const unsigned long no_of_b,
                     float offset,
                     NN_scores* this,
                     const char* file, const int line)
{
   unsigned long bp1, bp2;

   /* allocate memory */
   this->G_int11 = (float****) XOBJ_MALLOC_ND(sizeof (****this->G_int11),
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

   /* BEGIN_INT11_ENERGIES */
   /* CG */
   bp1 = this->bp_idx[c][g];
   /*    CG */
   bp2 = this->bp_idx[c][g];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =   40 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =   40 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =   40 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =   40 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =   40 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =   40 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =   40 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =   40 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =   40 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -140 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =   40 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =   40 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =   40 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =   40 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =   40 - offset; /*   U */
   /*    GC */
   bp2 = this->bp_idx[g][c];
   this->G_int11[bp1][bp2][a][a] =   40 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  -40 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =   40 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =   40 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =   30 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =   50 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =   40 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =   50 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  -10 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =   40 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -170 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =   40 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =   40 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =    0 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =   40 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  -30 - offset; /*   U */
   /*    GU */
   bp2 = this->bp_idx[g][u];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110 - offset; /*   U */
   /*    UG */
   bp2 = this->bp_idx[u][g];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110 - offset; /*   U */
   /*    AU */
   bp2 = this->bp_idx[a][u];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110 - offset; /*   U */
   /*    UA */
   bp2 = this->bp_idx[u][a];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110 - offset; /*   U */

   /* GC */
   bp1 = this->bp_idx[g][c];
   /*    CG */
   bp2 = this->bp_idx[c][g];
   this->G_int11[bp1][bp2][a][a] =   40 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =   30 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  -10 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =   40 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  -40 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =   50 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =   40 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =    0 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =   40 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =   40 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -170 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =   40 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =   40 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =   50 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =   40 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  -30 - offset; /*   U */
   /*    GC */
   bp2 = this->bp_idx[g][c];
   this->G_int11[bp1][bp2][a][a] =   80 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =   40 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =   40 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =   40 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =   40 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =   40 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =   40 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =   40 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =   40 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =   40 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -210 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =   40 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =   40 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =   40 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =   40 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  -70 - offset; /*   U */
   /*    GU */
   bp2 = this->bp_idx[g][u];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110 - offset; /*   U */
   /*    UG */
   bp2 = this->bp_idx[u][g];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110 - offset; /*   U */
   /*    AU */
   bp2 = this->bp_idx[a][u];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  100 - offset; /*   U */
   /*    UA */
   bp2 = this->bp_idx[u][a];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110 - offset; /*   U */

   /* GU */
   bp1 = this->bp_idx[g][u];
   /*    CG */
   bp2 = this->bp_idx[c][g];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110 - offset; /*   U */
   /*    GC */
   bp2 = this->bp_idx[g][c];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110 - offset; /*   U */
   /*    GU */
   bp2 = this->bp_idx[g][u];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170 - offset; /*   U */
   /*    UG */
   bp2 = this->bp_idx[u][g];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170 - offset; /*   U */
   /*    AU */
   bp2 = this->bp_idx[a][u];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170 - offset; /*   U */
   /*    UA */
   bp2 = this->bp_idx[u][a];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170 - offset; /*   U */

   /* UG */
   bp1 = this->bp_idx[u][g];
   /*    CG */
   bp2 = this->bp_idx[c][g];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110 - offset; /*   U */
   /*    GC */
   bp2 = this->bp_idx[g][c];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110 - offset; /*   U */
   /*    GU */
   bp2 = this->bp_idx[g][u];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170 - offset; /*   U */
   /*    UG */
   bp2 = this->bp_idx[u][g];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170 - offset; /*   U */
   /*    AU */
   bp2 = this->bp_idx[a][u];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170 - offset; /*   U */
   /*    UA */
   bp2 = this->bp_idx[u][a];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170 - offset; /*   U */

   /* AU */
   bp1 = this->bp_idx[a][u];
   /*    CG */
   bp2 = this->bp_idx[c][g];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110 - offset; /*   U */
   /*    GC */
   bp2 = this->bp_idx[g][c];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  100 - offset; /*   U */
   /*    GU */
   bp2 = this->bp_idx[g][u];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170 - offset; /*   U */
   /*    UG */
   bp2 = this->bp_idx[u][g];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170 - offset; /*   U */
   /*    AU */
   bp2 = this->bp_idx[a][u];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  120 - offset; /*   U */
   /*    UA */
   bp2 = this->bp_idx[u][a];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  150 - offset; /*   U */

   /* UA */
   bp1 = this->bp_idx[u][a];
   /*    CG */
   bp2 = this->bp_idx[c][g];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110 - offset; /*   U */
   /*    GC */
   bp2 = this->bp_idx[g][c];
   this->G_int11[bp1][bp2][a][a] =  110 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  110 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  110 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] = -100 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  110 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  110 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  110 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  110 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  110 - offset; /*   U */
   /*    GU */
   bp2 = this->bp_idx[g][u];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170 - offset; /*   U */
   /*    UG */
   bp2 = this->bp_idx[u][g];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  170 - offset; /*   U */
   /*    AU */
   bp2 = this->bp_idx[a][u];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  150 - offset; /*   U */
   /*    UA */
   bp2 = this->bp_idx[u][a];
   this->G_int11[bp1][bp2][a][a] =  170 - offset; /* A A */
   this->G_int11[bp1][bp2][a][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][a][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][a][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][c][a] =  170 - offset; /* C A */
   this->G_int11[bp1][bp2][c][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][c][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][c][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][g][a] =  170 - offset; /* G A */
   this->G_int11[bp1][bp2][g][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][g][g] =  -40 - offset; /*   G */
   this->G_int11[bp1][bp2][g][u] =  170 - offset; /*   U */
   this->G_int11[bp1][bp2][u][a] =  170 - offset; /* U A */
   this->G_int11[bp1][bp2][u][c] =  170 - offset; /*   C */
   this->G_int11[bp1][bp2][u][g] =  170 - offset; /*   G */
   this->G_int11[bp1][bp2][u][u] =  180 - offset; /*   U */

   /* END_INT11_ENERGIES */

   return 0;
}

static int
allocate_init_non_gc_penalty_for_bp (int a, int u, int g, int c,
                                     float offset,
                                     NN_scores* this,
                                     const char* file, const int line)
{
   this->non_gc_penalty_for_bp = /*(float*)*/ XOBJ_MALLOC (
      sizeof (*this->non_gc_penalty_for_bp)
      * this->bp_allowed_size,
      file, line);
   if (this->non_gc_penalty_for_bp == NULL)
   {
      return 1;
   }

   this->non_gc_penalty_for_bp[(int)this->bp_idx[a][u]] = 50 - offset;
   this->non_gc_penalty_for_bp[(int)this->bp_idx[c][g]] =  0 - offset;
   this->non_gc_penalty_for_bp[(int)this->bp_idx[g][c]] =  0 - offset;
   this->non_gc_penalty_for_bp[(int)this->bp_idx[g][u]] = 50 - offset;
   this->non_gc_penalty_for_bp[(int)this->bp_idx[u][a]] = 50 - offset;
   this->non_gc_penalty_for_bp[(int)this->bp_idx[u][g]] = 50 - offset;

   return 0;
}

static int
allocate_init_bulge_loop (float offset, NN_scores* this,
                          const char* file, const int line)
{

   this->G_bulge_loop_size = 31;
   this->G_bulge_loop = XOBJ_MALLOC (  this->G_bulge_loop_size
                                       * sizeof (this->G_bulge_loop[0]),
                                         file, line);
   if (this->G_bulge_loop == NULL)
   {
      return 1;        
   }

   /* BEGIN_BULGES */
   this->G_bulge_loop[ 0] = FLOAT_UNDEF;
   this->G_bulge_loop[ 1] = 380 - offset;
   this->G_bulge_loop[ 2] = 280 - offset;
   this->G_bulge_loop[ 3] = 320 - offset;
   this->G_bulge_loop[ 4] = 360 - offset;
   this->G_bulge_loop[ 5] = 400 - offset;
   this->G_bulge_loop[ 6] = 440 - offset;
   this->G_bulge_loop[ 7] = 459 - offset;
   this->G_bulge_loop[ 8] = 470 - offset;
   this->G_bulge_loop[ 9] = 480 - offset;
   this->G_bulge_loop[10] = 490 - offset;
   this->G_bulge_loop[11] = 500 - offset;
   this->G_bulge_loop[12] = 510 - offset;
   this->G_bulge_loop[13] = 519 - offset;
   this->G_bulge_loop[14] = 527 - offset;
   this->G_bulge_loop[15] = 534 - offset;
   this->G_bulge_loop[16] = 541 - offset;
   this->G_bulge_loop[17] = 548 - offset;
   this->G_bulge_loop[18] = 554 - offset;
   this->G_bulge_loop[19] = 560 - offset;
   this->G_bulge_loop[20] = 565 - offset;
   this->G_bulge_loop[21] = 571 - offset;
   this->G_bulge_loop[22] = 576 - offset;
   this->G_bulge_loop[23] = 580 - offset;
   this->G_bulge_loop[24] = 585 - offset;
   this->G_bulge_loop[25] = 589 - offset;
   this->G_bulge_loop[26] = 594 - offset;
   this->G_bulge_loop[27] = 598 - offset;
   this->G_bulge_loop[28] = 602 - offset;
   this->G_bulge_loop[29] = 605 - offset;
   this->G_bulge_loop[30] = 609 - offset;
   /* END_BULGES */

   return 0;
}

static int
allocate_init_int21 (const int a, const int u, const int g, const int c,
                     const unsigned long no_of_b,
                     float offset,
                     NN_scores* this,
                     const char* file, const int line)
{
   unsigned long bp1, bp2;

   /* allocate memory */
   this->G_int21 = (float*****) XOBJ_MALLOC_ND(sizeof (*****this->G_int21),
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

   /* BEGIN_INT21_ENERGIES */
   /* CG */
   bp1 = this->bp_idx[c][g];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int21[bp1][bp2][a][a][a] = 240 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 220 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 210 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 170 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 100 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] =  60 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] =  40 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 230 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 220 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 220 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 220 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 220 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 400 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 250 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 190 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 220 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 170 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] =  80 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 400 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] =  80 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 220 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 400 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 400 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 220 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 130 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 400 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 170 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 120 - offset; /*     U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int21[bp1][bp2][a][a][a] = 230 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 220 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 110 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 210 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 170 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] =  80 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] =  60 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] =  40 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 230 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 220 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 220 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 220 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 220 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 400 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 250 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 190 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 220 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 170 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] =  80 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 400 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] =  80 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 220 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 400 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 400 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 220 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 150 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 400 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 170 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 120 - offset; /*     U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */

   /* GC */
   bp1 = this->bp_idx[g][c];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int21[bp1][bp2][a][a][a] = 250 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 220 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 210 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 210 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 170 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 120 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] =  60 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] =  40 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 230 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 220 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 220 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 220 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 220 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 400 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 250 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 190 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 220 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 170 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] =  80 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 400 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] =  80 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 220 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 400 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 400 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 220 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 120 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 400 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 170 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 120 - offset; /*     U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int21[bp1][bp2][a][a][a] = 240 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 220 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 210 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 170 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 100 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] =  60 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] =  40 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 230 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 220 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 220 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 220 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 220 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 400 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 250 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 190 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 220 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 170 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] =  80 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 400 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] =  80 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 220 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 400 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 400 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 220 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 130 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 400 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 400 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 170 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 400 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 120 - offset; /*     U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */

   /* GU */
   bp1 = this->bp_idx[g][u];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */

   /* UG */
   bp1 = this->bp_idx[u][g];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */

   /* AU */
   bp1 = this->bp_idx[a][u];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */

   /* UA */
   bp1 = this->bp_idx[u][a];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int21[bp1][bp2][a][a][a] = 320 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 290 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 240 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 180 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 140 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 120 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 310 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 300 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 330 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 330 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 270 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 300 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 250 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 160 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 160 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 300 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 480 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 480 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 300 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 210 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 480 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 480 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 480 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 480 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 250 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 480 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 200 - offset; /*     U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int21[bp1][bp2][a][a][a] = 390 - offset; /* A A A */
   this->G_int21[bp1][bp2][a][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][a][a][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][c][a] = 360 - offset; /*   C A */
   this->G_int21[bp1][bp2][a][c][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][a][c][g] = 310 - offset; /*     G */
   this->G_int21[bp1][bp2][a][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][g][a] = 250 - offset; /*   G A */
   this->G_int21[bp1][bp2][a][g][c] = 210 - offset; /*     C */
   this->G_int21[bp1][bp2][a][g][g] = 190 - offset; /*     G */
   this->G_int21[bp1][bp2][a][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][a][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][a][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][a][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][a][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][a][a] = 380 - offset; /* C A A */
   this->G_int21[bp1][bp2][c][a][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][c][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][a][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][c][a] = 370 - offset; /*   C A */
   this->G_int21[bp1][bp2][c][c][c] = 400 - offset; /*     C */
   this->G_int21[bp1][bp2][c][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][c][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][c][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][c][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][c][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][c][u][a] = 400 - offset; /*   U A */
   this->G_int21[bp1][bp2][c][u][c] = 340 - offset; /*     C */
   this->G_int21[bp1][bp2][c][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][c][u][u] = 370 - offset; /*     U */
   this->G_int21[bp1][bp2][g][a][a] = 320 - offset; /* G A A */
   this->G_int21[bp1][bp2][g][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][a][g] = 230 - offset; /*     G */
   this->G_int21[bp1][bp2][g][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][g][c][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][c][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][g][a] = 230 - offset; /*   G A */
   this->G_int21[bp1][bp2][g][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][g][g] = 370 - offset; /*     G */
   this->G_int21[bp1][bp2][g][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][g][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][g][u][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][g][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][g][u][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][a][a] = 550 - offset; /* U A A */
   this->G_int21[bp1][bp2][u][a][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][a][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][a][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][c][a] = 550 - offset; /*   C A */
   this->G_int21[bp1][bp2][u][c][c] = 370 - offset; /*     C */
   this->G_int21[bp1][bp2][u][c][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][c][u] = 280 - offset; /*     U */
   this->G_int21[bp1][bp2][u][g][a] = 550 - offset; /*   G A */
   this->G_int21[bp1][bp2][u][g][c] = 550 - offset; /*     C */
   this->G_int21[bp1][bp2][u][g][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][g][u] = 550 - offset; /*     U */
   this->G_int21[bp1][bp2][u][u][a] = 550 - offset; /*   U A */
   this->G_int21[bp1][bp2][u][u][c] = 320 - offset; /*     C */
   this->G_int21[bp1][bp2][u][u][g] = 550 - offset; /*     G */
   this->G_int21[bp1][bp2][u][u][u] = 270 - offset; /*     U */

   /* END_INT21_ENERGIES */

   return 0;
}

static int
allocate_init_int22 (const int a, const int u, const int g, const int c,
                     const unsigned long no_of_b,
                     float offset,
                     NN_scores* this,
                     const char* file, const int line)
{
   unsigned long bp1, bp2;

   /* allocate memory */
   this->G_int22 = (float******) XOBJ_MALLOC_ND(sizeof (******this->G_int22),
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

   /* BEGIN_INT22_ENERGIES */
   /* CG */
   bp1 = this->bp_idx[c][g];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  130 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  120 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   20 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   30 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  160 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   60 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  210 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  190 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   30 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  -40 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -60 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   40 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  140 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  120 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  110 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   20 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  150 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  140 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  120 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   20 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -150 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  -20 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -40 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  120 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =    0 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   30 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   50 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   30 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   20 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =    0 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -70 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   60 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  150 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  130 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -70 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =    0 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -60 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -260 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  130 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  130 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  -40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -150 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -60 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   60 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =   50 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -60 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -40 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -420 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  -50 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   50 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -260 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =   50 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =  -40 - offset; /*       U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int22[bp1][bp2][a][a][a][a] =   50 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  110 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  -30 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   10 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] = -160 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  110 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] = -100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  -60 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =   70 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =   50 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =   20 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   40 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   50 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   10 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -70 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -80 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  180 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  -50 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  -60 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  -60 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  150 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  130 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  100 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  -70 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  100 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  100 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =    0 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   70 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -190 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  -30 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -70 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  110 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] = -150 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  -20 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  -20 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  -10 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   20 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  -20 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  -40 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =    0 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] = -170 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   70 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  -10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  110 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   20 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =    0 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -50 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =   60 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =    0 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -90 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] = -100 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -300 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  120 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] = -130 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  -60 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -240 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =   90 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =   60 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -160 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   40 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] = -160 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =   50 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -160 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  -60 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   20 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] = -160 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =   50 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -60 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -70 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -440 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] = -100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  -10 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -410 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =   40 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] = -100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   60 - offset; /*       U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  200 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  180 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   80 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  190 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  230 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   80 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   60 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  160 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  190 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   70 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  180 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  160 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   70 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  120 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   80 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  100 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   80 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   60 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -20 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   90 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  170 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -20 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   30 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   40 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  140 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -60 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  180 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -80 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -310 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   60 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -210 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   10 - offset; /*       U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  200 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   60 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  260 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  260 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  100 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  140 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   20 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -40 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  150 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  190 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  190 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   90 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -90 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   40 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   40 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  110 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  100 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   60 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =    0 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  100 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =    0 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   40 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =    0 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -70 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   20 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -220 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30 - offset; /*       U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  200 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  180 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   80 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  190 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  230 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   80 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   60 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  160 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  190 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   70 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  180 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  160 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   70 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  120 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   80 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  100 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   80 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   60 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -20 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   90 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  170 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -20 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   30 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   40 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  140 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -60 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  180 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -80 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -310 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   60 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -210 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   10 - offset; /*       U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int22[bp1][bp2][a][a][a][a] =  200 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   60 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  260 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  260 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  100 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  140 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   20 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -40 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  150 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  190 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  190 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   90 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -90 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   40 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   40 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  110 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  100 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   60 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =    0 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  100 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =    0 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   40 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =    0 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -70 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   20 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -220 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30 - offset; /*       U */

   /* GC */
   bp1 = this->bp_idx[g][c];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int22[bp1][bp2][a][a][a][a] =   50 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  130 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  -20 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =   60 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =    0 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] = -100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  -10 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -60 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  -10 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =   60 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  140 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  110 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  100 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  -40 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  150 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  130 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  120 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  -70 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  -60 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -160 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  -60 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -50 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  120 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =    0 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   30 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  -30 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  -50 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  -70 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] = -150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] = -170 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] = -130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   10 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  -60 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =   70 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  -20 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   40 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  -10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] = -160 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  -60 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  -90 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  -60 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] = -160 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -410 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =   40 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   30 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -240 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =   50 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -190 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  -20 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   20 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =   50 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =    0 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  -10 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   20 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =   20 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -80 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =   50 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] = -100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] = -160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -440 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] = -100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   20 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -300 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =   60 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =   10 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] = -100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   60 - offset; /*       U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int22[bp1][bp2][a][a][a][a] =  150 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  120 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  -50 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  -80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] = -190 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  120 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  180 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =   80 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   10 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  -20 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] = -130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  -70 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] = -130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  -30 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] = -160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  150 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =   20 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  120 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  100 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  -80 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =   20 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =   90 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  100 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =    0 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =    0 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  -10 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -190 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  -90 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -90 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  100 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] = -150 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  -50 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  -50 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   20 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   20 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  -50 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  -80 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] = -150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] = -260 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] = -150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  -80 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] = -160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =   20 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  -50 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  -80 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] = -150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] = -190 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  -90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  -90 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  -60 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] = -190 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] = -100 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -450 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =   30 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  -50 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] = -150 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  -60 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -410 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =   30 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  -50 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -190 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   20 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  -80 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] = -190 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =   20 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   20 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] = -100 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] = -190 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  -70 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  -10 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] = -130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =   50 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] = -100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] = -190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -490 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  -90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] = -150 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -450 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  -50 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  -70 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =  -50 - offset; /*       U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  210 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   60 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  240 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   90 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  140 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   70 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -60 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  170 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   40 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   60 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  170 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  160 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   60 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  120 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   50 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   50 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   40 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   80 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   10 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  -10 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  -90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] = -110 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  -10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =    0 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  -80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =   80 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   80 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] = -110 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  -10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  -60 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -50 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -310 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =   80 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -210 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =   80 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -120 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -50 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   50 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   50 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -250 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  -20 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =    0 - offset; /*       U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  210 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  210 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  270 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  110 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   30 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  160 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  230 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  140 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   80 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  190 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  190 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   80 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  120 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -90 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   10 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   30 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  100 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   10 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  -30 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] = -110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -90 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   10 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  -60 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  110 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  110 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -90 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  -50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -90 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -350 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =   80 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   40 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  110 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =   40 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -150 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -390 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   40 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -260 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30 - offset; /*       U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  210 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   60 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  240 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   90 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  140 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   70 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -60 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  170 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   40 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   60 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  170 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  160 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   60 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  120 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   50 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   50 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   40 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   80 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   10 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  -10 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  -90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] = -110 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  -10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =    0 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  -80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =   80 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   80 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] = -110 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  -10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  -60 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -50 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -310 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =   80 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -210 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =   80 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -120 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -50 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   50 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   50 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -250 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  -20 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =    0 - offset; /*       U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int22[bp1][bp2][a][a][a][a] =  210 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  210 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  270 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  110 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   30 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  160 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  230 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  140 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   80 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  190 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  190 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   80 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  120 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -90 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   10 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   30 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  100 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   10 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  -30 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] = -110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -90 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   10 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  -60 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  110 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  110 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -90 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  -50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -90 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -350 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =   80 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   40 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  110 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =   40 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] = -150 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -390 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   40 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -260 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30 - offset; /*       U */

   /* GU */
   bp1 = this->bp_idx[g][u];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  200 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  240 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  280 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  100 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   30 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  270 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  180 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   80 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  210 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  180 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   80 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -90 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =    0 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  180 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   60 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   90 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   80 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   70 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   50 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -20 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  110 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  180 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -20 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -10 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -40 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -210 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  180 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  140 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -60 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  180 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   60 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -310 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =    0 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -160 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   60 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   10 - offset; /*       U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int22[bp1][bp2][a][a][a][a] =  210 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   60 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   10 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] = -110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  180 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   60 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  150 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   70 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =    0 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -120 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =   40 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  -80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  210 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  -10 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  160 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =   60 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =   60 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   50 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -130 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  -30 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  -90 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =   10 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   90 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   90 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   60 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] = -110 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  -60 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   50 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   60 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   70 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   50 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -50 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -50 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -250 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  -10 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -210 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  170 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -60 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  -20 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -30 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   10 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -310 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =    0 - offset; /*       U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   20 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  240 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =   10 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  140 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   30 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  140 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   30 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   80 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   90 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  190 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   20 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -210 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   60 - offset; /*       U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  250 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  250 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -20 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  250 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  250 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  250 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  150 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  190 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -30 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  170 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  110 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  150 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  250 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   50 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  250 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -120 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80 - offset; /*       U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   20 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  240 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =   10 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  140 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   30 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  140 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   30 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   80 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   90 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  190 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   20 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -210 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   60 - offset; /*       U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  250 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  250 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -20 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  250 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  250 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  250 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  150 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  190 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -30 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  170 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  110 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  150 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  250 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   50 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  250 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -120 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80 - offset; /*       U */

   /* UG */
   bp1 = this->bp_idx[u][g];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  200 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  240 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  280 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  100 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   30 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -70 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  270 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  160 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  150 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   60 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  190 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  180 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  160 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   60 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -110 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =    0 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   40 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   70 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  100 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   90 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =    0 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  130 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =    0 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   70 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   10 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  160 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -70 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  190 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   20 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   20 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30 - offset; /*       U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int22[bp1][bp2][a][a][a][a] =  210 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   10 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   10 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  180 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  150 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   70 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =    0 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =   40 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  -60 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  210 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  170 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  140 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  -30 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  140 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  140 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =   40 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =   60 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   30 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -150 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  150 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] = -110 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =   40 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   70 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  110 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   80 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -90 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   80 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  180 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   90 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -30 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   70 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -30 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -260 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   10 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  190 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =    0 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -390 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  -10 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   40 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -350 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30 - offset; /*       U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  220 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  110 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  160 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -10 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  120 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  150 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  160 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  240 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  110 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  240 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -150 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80 - offset; /*       U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  230 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  230 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  170 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   80 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  150 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   70 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  170 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  270 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   70 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   70 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  230 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -290 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -160 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =  110 - offset; /*       U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  220 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  110 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  160 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -10 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  120 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  150 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  160 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  240 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  110 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  240 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -150 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80 - offset; /*       U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  230 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  230 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  170 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   80 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  150 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   70 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  170 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  270 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   70 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   70 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  230 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -290 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -160 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =  110 - offset; /*       U */

   /* AU */
   bp1 = this->bp_idx[a][u];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  200 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  240 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  280 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  100 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   30 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  270 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  180 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   80 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  210 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  180 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   80 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -90 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =    0 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  180 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   60 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   90 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   80 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   70 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   50 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -20 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  110 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  180 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -20 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -10 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -40 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -210 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  180 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  140 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -60 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  180 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   60 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -310 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =    0 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -160 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   60 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   10 - offset; /*       U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int22[bp1][bp2][a][a][a][a] =  210 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   60 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   10 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] = -110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  180 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   60 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  150 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   70 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =    0 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -120 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =   40 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  -80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  210 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  190 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  -10 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  160 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =   60 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =   60 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   50 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -130 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  -30 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  -90 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =   10 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   90 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =   90 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   60 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] = -110 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  -60 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   50 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   60 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   70 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   50 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -50 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -50 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -250 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  170 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  -10 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -210 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  170 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -60 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   60 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  -20 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  -30 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   10 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -310 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =    0 - offset; /*       U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   20 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  240 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =   10 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  140 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   30 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  140 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   30 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   80 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   90 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  190 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   20 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -210 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   60 - offset; /*       U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  250 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  250 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -20 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  250 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  250 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  250 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  150 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  190 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -30 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  170 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  110 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  150 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  250 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   50 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  250 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -120 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80 - offset; /*       U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   20 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  240 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =   10 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  140 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   30 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  140 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   30 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   80 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   90 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  190 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -10 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   20 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -210 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  140 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   60 - offset; /*       U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  250 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  250 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -20 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  250 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  250 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  250 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  250 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  150 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  190 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -30 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  170 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  150 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  110 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  150 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  250 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   50 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  250 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  310 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -120 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  160 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80 - offset; /*       U */

   /* UA */
   bp1 = this->bp_idx[u][a];
   /*   CG */
   bp2 = this->bp_idx[c][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  200 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  240 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  280 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  100 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =   30 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =  -70 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  270 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =   80 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  220 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  160 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  150 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =   60 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  190 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  180 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  160 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   60 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -110 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =    0 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  160 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   40 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   70 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  100 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   90 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =    0 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  130 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  220 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  220 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =    0 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   70 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   10 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  160 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -70 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  190 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   20 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   40 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -350 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   20 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30 - offset; /*       U */
   /*   GC */
   bp2 = this->bp_idx[g][c];
   this->G_int22[bp1][bp2][a][a][a][a] =  210 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =   10 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =   10 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  180 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  150 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =   70 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =   70 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  180 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =    0 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] = -100 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =   40 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  -60 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  210 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  170 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  140 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  -30 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  140 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  140 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =   40 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =   60 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =   30 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  210 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] = -150 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  150 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] = -110 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =   40 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =   10 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =   70 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  110 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =   80 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =  -90 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =   80 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  120 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  180 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  180 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =   90 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =   80 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =   40 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =  -30 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =   70 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  -30 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =  -70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -260 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =   10 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] = -220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  190 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  180 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] = -100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =   90 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  -90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =    0 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  240 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -390 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  -10 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =   40 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -350 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =   30 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  110 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   10 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   30 - offset; /*       U */
   /*   GU */
   bp2 = this->bp_idx[g][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  220 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  110 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  160 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -10 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  120 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  150 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  160 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  240 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  110 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  240 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -150 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80 - offset; /*       U */
   /*   UG */
   bp2 = this->bp_idx[u][g];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  230 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  230 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  170 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   80 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  150 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   70 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  170 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  270 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   70 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   70 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  230 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -290 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -160 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =  110 - offset; /*       U */
   /*   AU */
   bp2 = this->bp_idx[a][u];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  250 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  150 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  260 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  310 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  150 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  150 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  130 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =   30 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  230 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  210 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  220 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  110 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  160 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -10 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   90 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  120 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  140 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  150 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  190 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =  100 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   50 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  160 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  110 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  240 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   50 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  100 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  140 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =  110 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -20 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  240 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  170 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  260 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =  -20 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  180 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  130 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -250 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =  130 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -150 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =   80 - offset; /*       U */
   /*   UA */
   bp2 = this->bp_idx[u][a];
   this->G_int22[bp1][bp2][a][a][a][a] =  280 - offset; /* A A A A */
   this->G_int22[bp1][bp2][a][a][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][a][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][c][g] =  130 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][g][a] =  170 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][a][g][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][a][a] =  280 - offset; /*   C A A */
   this->G_int22[bp1][bp2][a][c][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][a][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][c][a] =  340 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][c][c][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][c][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][c][u][a] =  340 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][c][u][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][a][c][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][a][c][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][a][a] =  170 - offset; /*   G A A */
   this->G_int22[bp1][bp2][a][g][a][c] =  170 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][g][a] =  210 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][g][g][c] =  210 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][g][u][a] =  100 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][g][u][c] =    0 - offset; /*       C */
   this->G_int22[bp1][bp2][a][g][u][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][a][g][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][a][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][c][a] =  310 - offset; /*     C A */
   this->G_int22[bp1][bp2][a][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][g][a] =  220 - offset; /*     G A */
   this->G_int22[bp1][bp2][a][u][g][c] =  120 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][a][u][u][a] =  290 - offset; /*     U A */
   this->G_int22[bp1][bp2][a][u][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][a][u][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][a][u][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][a][a] =  230 - offset; /* C A A A */
   this->G_int22[bp1][bp2][c][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][a][u] =  310 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][c][a] =  190 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][a][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][g][a] =  130 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][a][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][a][a] =  230 - offset; /*   C A A */
   this->G_int22[bp1][bp2][c][c][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][a][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][c][a] =  230 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][c][c][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][c][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][c][u][a] =  230 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][c][u][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][c][c][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][c][u][u] =  250 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][a][a] =  130 - offset; /*   G A A */
   this->G_int22[bp1][bp2][c][g][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][a][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][g][a] =  170 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][g][u] =  240 - offset; /*       U */
   this->G_int22[bp1][bp2][c][g][u][a] =  -50 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][c][g][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][g][u][u] =   70 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][c][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][c][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][g][a] =   80 - offset; /*     G A */
   this->G_int22[bp1][bp2][c][u][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][g][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][c][u][u][a] =  150 - offset; /*     U A */
   this->G_int22[bp1][bp2][c][u][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][c][u][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][c][u][u][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][a][a] =  170 - offset; /* G A A A */
   this->G_int22[bp1][bp2][g][a][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][a][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][c][a] =  130 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][a][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][c][g] =  170 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][c][u] =   80 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][g][a] =   70 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][a][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][g][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][a][a] =  170 - offset; /*   C A A */
   this->G_int22[bp1][bp2][g][c][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][a][g] =  210 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][a][u] =  120 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][c][a] =  270 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][c][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][c][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][c][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][c][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][c][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][c][u][g] =  270 - offset; /*       G */
   this->G_int22[bp1][bp2][g][c][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][a][a] =   70 - offset; /*   G A A */
   this->G_int22[bp1][bp2][g][g][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][a][g] =  110 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][g][a] =  110 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][g][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][g][g] =  150 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][g][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][g][u][a] =   70 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][g][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][g][u][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][g][g][u][u] = -160 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][g][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][c][a] =  240 - offset; /*     C A */
   this->G_int22[bp1][bp2][g][u][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][c][g] =  240 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][g][u][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][g][g] =  160 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][g][u] =  -30 - offset; /*       U */
   this->G_int22[bp1][bp2][g][u][u][a] =  270 - offset; /*     U A */
   this->G_int22[bp1][bp2][g][u][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][g][u][u][g] =  230 - offset; /*       G */
   this->G_int22[bp1][bp2][g][u][u][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][a][a] =  200 - offset; /* U A A A */
   this->G_int22[bp1][bp2][u][a][a][c] =  340 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][a][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][a][u] =  290 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][a][c][c] =  230 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][c][g] =  -50 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][c][u] =  150 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][a][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][g][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][g][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][a][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][a][u][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][a][u][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][a][u][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][a][a] =  200 - offset; /*   C A A */
   this->G_int22[bp1][bp2][u][c][a][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][a][g] =    0 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][a][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][c][c][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][c][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][c][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][c][g][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][g][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][g][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][c][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][c][u][c] =  280 - offset; /*       C */
   this->G_int22[bp1][bp2][u][c][u][g] =  100 - offset; /*       G */
   this->G_int22[bp1][bp2][u][c][u][u] =  190 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][a][a] =  200 - offset; /*   G A A */
   this->G_int22[bp1][bp2][u][g][a][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][a][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][a][u] =  270 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][g][c][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][c][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][c][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][g][g][c] =  270 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][g][g] =   30 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][g][u] =  230 - offset; /*       U */
   this->G_int22[bp1][bp2][u][g][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][g][u][c] =  100 - offset; /*       C */
   this->G_int22[bp1][bp2][u][g][u][g] = -290 - offset; /*       G */
   this->G_int22[bp1][bp2][u][g][u][u] =   90 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][a][a] =  200 - offset; /*   U A A */
   this->G_int22[bp1][bp2][u][u][a][c] =  200 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][a][g] =  200 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][a][u] =  200 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][c][a] =  200 - offset; /*     C A */
   this->G_int22[bp1][bp2][u][u][c][c] =  250 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][c][g] =   70 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][c][u] =  160 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][g][a] =  200 - offset; /*     G A */
   this->G_int22[bp1][bp2][u][u][g][c] =  220 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][g][g] = -160 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][g][u] =  220 - offset; /*       U */
   this->G_int22[bp1][bp2][u][u][u][a] =  200 - offset; /*     U A */
   this->G_int22[bp1][bp2][u][u][u][c] =  190 - offset; /*       C */
   this->G_int22[bp1][bp2][u][u][u][g] =   90 - offset; /*       G */
   this->G_int22[bp1][bp2][u][u][u][u] =  110 - offset; /*       U */

   /* END_INT22_ENERGIES */

   return 0;
}

static int
allocate_init_dangle5 (const int a, const int u, const int g, const int c, 
                       const unsigned long size,
                       float offset,
                       NN_scores* this,
                       const char* file, const int line)
{
   this->G_dangle5 = (float**) XOBJ_MALLOC_2D (this->bp_allowed_size, size,
                                             sizeof (**this->G_dangle5),
                                             file, line);

   this->G_dangle5_size = this->bp_allowed_size * size;

   if (this->G_dangle5 == NULL)
   {
      return 1;
   }   

   /* BEGIN_DANGLE_5 */
   /* CG */
   this->G_dangle5[(int)this->bp_idx[c][g]][a] = -50 - offset;
   this->G_dangle5[(int)this->bp_idx[c][g]][c] = -30 - offset;
   this->G_dangle5[(int)this->bp_idx[c][g]][g] = -20 - offset;
   this->G_dangle5[(int)this->bp_idx[c][g]][u] = -10 - offset;

   /* GC */
   this->G_dangle5[(int)this->bp_idx[g][c]][a] = -20 - offset;
   this->G_dangle5[(int)this->bp_idx[g][c]][c] = -30 - offset;
   this->G_dangle5[(int)this->bp_idx[g][c]][g] =   0 - offset;
   this->G_dangle5[(int)this->bp_idx[g][c]][u] =   0 - offset;

   /* GU */
   this->G_dangle5[(int)this->bp_idx[g][u]][a] = -30 - offset;
   this->G_dangle5[(int)this->bp_idx[g][u]][c] = -30 - offset;
   this->G_dangle5[(int)this->bp_idx[g][u]][g] = -40 - offset;
   this->G_dangle5[(int)this->bp_idx[g][u]][u] = -20 - offset;

   /* UG */
   this->G_dangle5[(int)this->bp_idx[u][g]][a] = -30 - offset;
   this->G_dangle5[(int)this->bp_idx[u][g]][c] = -10 - offset;
   this->G_dangle5[(int)this->bp_idx[u][g]][g] = -20 - offset;
   this->G_dangle5[(int)this->bp_idx[u][g]][u] = -20 - offset;

   /* AU */
   this->G_dangle5[(int)this->bp_idx[a][u]][a] = -30 - offset;
   this->G_dangle5[(int)this->bp_idx[a][u]][c] = -30 - offset;
   this->G_dangle5[(int)this->bp_idx[a][u]][g] = -40 - offset;
   this->G_dangle5[(int)this->bp_idx[a][u]][u] = -20 - offset;

   /* UA */
   this->G_dangle5[(int)this->bp_idx[u][a]][a] = -30 - offset;
   this->G_dangle5[(int)this->bp_idx[u][a]][c] = -10 - offset;
   this->G_dangle5[(int)this->bp_idx[u][a]][g] = -20 - offset;
   this->G_dangle5[(int)this->bp_idx[u][a]][u] = -20 - offset;

   /* END_DANGLE_5 */

   return 0;
}

static int
allocate_init_dangle3 (const int a, const int u, const int g, const int c, 
                       const unsigned long size,
                       float offset,
                       NN_scores* this,
                       const char* file, const int line)
{
   this->G_dangle3 = (float**) XOBJ_MALLOC_2D (this->bp_allowed_size, size,
                                               sizeof (**this->G_dangle3),
                                               file, line);

   this->G_dangle3_size = this->bp_allowed_size * size;

   if (this->G_dangle3 == NULL)
   {
      return 1;
   }   

   /* BEGIN_DANGLE_3 */
   /* CG */
   this->G_dangle3[(int)this->bp_idx[c][g]][a] = -110 - offset;
   this->G_dangle3[(int)this->bp_idx[c][g]][c] =  -40 - offset;
   this->G_dangle3[(int)this->bp_idx[c][g]][g] = -130 - offset;
   this->G_dangle3[(int)this->bp_idx[c][g]][u] =  -60 - offset;

   /* GC */
   this->G_dangle3[(int)this->bp_idx[g][c]][a] = -170 - offset;
   this->G_dangle3[(int)this->bp_idx[g][c]][c] =  -80 - offset;
   this->G_dangle3[(int)this->bp_idx[g][c]][g] = -170 - offset;
   this->G_dangle3[(int)this->bp_idx[g][c]][u] = -120 - offset;

   /* GU */
   this->G_dangle3[(int)this->bp_idx[g][u]][a] =  -70 - offset;
   this->G_dangle3[(int)this->bp_idx[g][u]][c] =  -10 - offset;
   this->G_dangle3[(int)this->bp_idx[g][u]][g] =  -70 - offset;
   this->G_dangle3[(int)this->bp_idx[g][u]][u] =  -10 - offset;

   /* UG */
   this->G_dangle3[(int)this->bp_idx[u][g]][a] =  -80 - offset;
   this->G_dangle3[(int)this->bp_idx[u][g]][c] =  -50 - offset;
   this->G_dangle3[(int)this->bp_idx[u][g]][g] =  -80 - offset;
   this->G_dangle3[(int)this->bp_idx[u][g]][u] =  -60 - offset;

   /* AU */
   this->G_dangle3[(int)this->bp_idx[a][u]][a] =  -70 - offset;
   this->G_dangle3[(int)this->bp_idx[a][u]][c] =  -10 - offset;
   this->G_dangle3[(int)this->bp_idx[a][u]][g] =  -70 - offset;
   this->G_dangle3[(int)this->bp_idx[a][u]][u] =  -10 - offset;

   /* UA */
   this->G_dangle3[(int)this->bp_idx[u][a]][a] =  -80 - offset;
   this->G_dangle3[(int)this->bp_idx[u][a]][c] =  -50 - offset;
   this->G_dangle3[(int)this->bp_idx[u][a]][g] =  -80 - offset;
   this->G_dangle3[(int)this->bp_idx[u][a]][u] =  -60 - offset;

   /* END_DANGLE_3 */

   return 0;
}

static void
tetra_loop_swap_entries (unsigned long src, unsigned long dest,
                         NN_scores* this)
{
   unsigned long i;
   char tmp[D_TL];
   float G_tmp;
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
allocate_init_tetra_loop (const char a,
                          const char u,
                          const char g,
                          const char c,
                          float offset,
                          NN_scores* this,
                          const char* file, const int line)
{
   char* l;

   this->tetra_loop_size = 30;
   this->tetra_loop = (char**) XOBJ_MALLOC_2D (this->tetra_loop_size,(D_TL + 1),
                                               sizeof (**this->tetra_loop),
                                               file, line);
   if (this->tetra_loop == NULL)
   {
      return 1;
   }

   this->G_tetra_loop = (float*) XOBJ_MALLOC (sizeof (*this->G_tetra_loop)
                                            * this->tetra_loop_size,
                                            file, line);

   /* BEGIN_TETRA_LOOPS */
   /* AGAAAU -200 */
   l = this->tetra_loop[ 0];
   l[0] = a; l[1] = g; l[2] = a; l[3] = a; l[4] = a; l[5] = u; l[6] = '\0';
   this->G_tetra_loop[ 0] = -200 - offset;

   /* AGCAAU -150 */
   l = this->tetra_loop[ 1];
   l[0] = a; l[1] = g; l[2] = c; l[3] = a; l[4] = a; l[5] = u; l[6] = '\0';
   this->G_tetra_loop[ 1] = -150 - offset;

   /* AGUAAU -150 */
   l = this->tetra_loop[ 2];
   l[0] = a; l[1] = g; l[2] = u; l[3] = a; l[4] = a; l[5] = u; l[6] = '\0';
   this->G_tetra_loop[ 2] = -150 - offset;

   /* AGUGAU -150 */
   l = this->tetra_loop[ 3];
   l[0] = a; l[1] = g; l[2] = u; l[3] = g; l[4] = a; l[5] = u; l[6] = '\0';
   this->G_tetra_loop[ 3] = -150 - offset;

   /* CGAAAG -300 */
   l = this->tetra_loop[ 4];
   l[0] = c; l[1] = g; l[2] = a; l[3] = a; l[4] = a; l[5] = g; l[6] = '\0';
   this->G_tetra_loop[ 4] = -300 - offset;

   /* CGAAGG -250 */
   l = this->tetra_loop[ 5];
   l[0] = c; l[1] = g; l[2] = a; l[3] = a; l[4] = g; l[5] = g; l[6] = '\0';
   this->G_tetra_loop[ 5] = -250 - offset;

   /* CGAGAG -200 */
   l = this->tetra_loop[ 6];
   l[0] = c; l[1] = g; l[2] = a; l[3] = g; l[4] = a; l[5] = g; l[6] = '\0';
   this->G_tetra_loop[ 6] = -200 - offset;

   /* CGCAAG -300 */
   l = this->tetra_loop[ 7];
   l[0] = c; l[1] = g; l[2] = c; l[3] = a; l[4] = a; l[5] = g; l[6] = '\0';
   this->G_tetra_loop[ 7] = -300 - offset;

   /* CGCGAG -250 */
   l = this->tetra_loop[ 8];
   l[0] = c; l[1] = g; l[2] = c; l[3] = g; l[4] = a; l[5] = g; l[6] = '\0';
   this->G_tetra_loop[ 8] = -250 - offset;

   /* CGGAAG -300 */
   l = this->tetra_loop[ 9];
   l[0] = c; l[1] = g; l[2] = g; l[3] = a; l[4] = a; l[5] = g; l[6] = '\0';
   this->G_tetra_loop[ 9] = -300 - offset;

   /* CGGGAG -150 */
   l = this->tetra_loop[10];
   l[0] = c; l[1] = g; l[2] = g; l[3] = g; l[4] = a; l[5] = g; l[6] = '\0';
   this->G_tetra_loop[10] = -150 - offset;

   /* CGUAAG -200 */
   l = this->tetra_loop[11];
   l[0] = c; l[1] = g; l[2] = u; l[3] = a; l[4] = a; l[5] = g; l[6] = '\0';
   this->G_tetra_loop[11] = -200 - offset;

   /* CGUGAG -300 */
   l = this->tetra_loop[12];
   l[0] = c; l[1] = g; l[2] = u; l[3] = g; l[4] = a; l[5] = g; l[6] = '\0';
   this->G_tetra_loop[12] = -300 - offset;

   /* CUAACG -200 */
   l = this->tetra_loop[13];
   l[0] = c; l[1] = u; l[2] = a; l[3] = a; l[4] = c; l[5] = g; l[6] = '\0';
   this->G_tetra_loop[13] = -200 - offset;

   /* CUACGG -250 */
   l = this->tetra_loop[14];
   l[0] = c; l[1] = u; l[2] = a; l[3] = c; l[4] = g; l[5] = g; l[6] = '\0';
   this->G_tetra_loop[14] = -250 - offset;

   /* CUUCGG -300 */
   l = this->tetra_loop[15];
   l[0] = c; l[1] = u; l[2] = u; l[3] = c; l[4] = g; l[5] = g; l[6] = '\0';
   this->G_tetra_loop[15] = -300 - offset;

   /* GGAAAC -300 */
   l = this->tetra_loop[16];
   l[0] = g; l[1] = g; l[2] = a; l[3] = a; l[4] = a; l[5] = c; l[6] = '\0';
   this->G_tetra_loop[16] = -300 - offset;

   /* GGAAGC -150 */
   l = this->tetra_loop[17];
   l[0] = g; l[1] = g; l[2] = a; l[3] = a; l[4] = g; l[5] = c; l[6] = '\0';
   this->G_tetra_loop[17] = -150 - offset;

   /* GGAGAC -300 */
   l = this->tetra_loop[18];
   l[0] = g; l[1] = g; l[2] = a; l[3] = g; l[4] = a; l[5] = c; l[6] = '\0';
   this->G_tetra_loop[18] = -300 - offset;

   /* GGCAAC -250 */
   l = this->tetra_loop[19];
   l[0] = g; l[1] = g; l[2] = c; l[3] = a; l[4] = a; l[5] = c; l[6] = '\0';
   this->G_tetra_loop[19] = -250 - offset;

   /* GGCGAC -150 */
   l = this->tetra_loop[20];
   l[0] = g; l[1] = g; l[2] = c; l[3] = g; l[4] = a; l[5] = c; l[6] = '\0';
   this->G_tetra_loop[20] = -150 - offset;

   /* GGGAAC -150 */
   l = this->tetra_loop[21];
   l[0] = g; l[1] = g; l[2] = g; l[3] = a; l[4] = a; l[5] = c; l[6] = '\0';
   this->G_tetra_loop[21] = -150 - offset;

   /* GGGAGC -150 */
   l = this->tetra_loop[22];
   l[0] = g; l[1] = g; l[2] = g; l[3] = a; l[4] = g; l[5] = c; l[6] = '\0';
   this->G_tetra_loop[22] = -150 - offset;

   /* GGGGAC -300 */
   l = this->tetra_loop[23];
   l[0] = g; l[1] = g; l[2] = g; l[3] = g; l[4] = a; l[5] = c; l[6] = '\0';
   this->G_tetra_loop[23] = -300 - offset;

   /* GGUGAC -300 */
   l = this->tetra_loop[24];
   l[0] = g; l[1] = g; l[2] = u; l[3] = g; l[4] = a; l[5] = c; l[6] = '\0';
   this->G_tetra_loop[24] = -300 - offset;

   /* GUGAAC -150 */
   l = this->tetra_loop[25];
   l[0] = g; l[1] = u; l[2] = g; l[3] = a; l[4] = a; l[5] = c; l[6] = '\0';
   this->G_tetra_loop[25] = -150 - offset;

   /* UGAAAA -150 */
   l = this->tetra_loop[26];
   l[0] = u; l[1] = g; l[2] = a; l[3] = a; l[4] = a; l[5] = a; l[6] = '\0';
   this->G_tetra_loop[26] = -150 - offset;

   /* UGAAAG -200 */
   l = this->tetra_loop[27];
   l[0] = u; l[1] = g; l[2] = a; l[3] = a; l[4] = a; l[5] = g; l[6] = '\0';
   this->G_tetra_loop[27] = -200 - offset;

   /* UGAGAG -250 */
   l = this->tetra_loop[28];
   l[0] = u; l[1] = g; l[2] = a; l[3] = g; l[4] = a; l[5] = g; l[6] = '\0';
   this->G_tetra_loop[28] = -250 - offset;

   /* UGGAAA -150 */
   l = this->tetra_loop[29];
   l[0] = u; l[1] = g; l[2] = g; l[3] = a; l[4] = a; l[5] = a; l[6] = '\0';
   this->G_tetra_loop[29] = -150 - offset;

   /* END_TETRA_LOOPS */

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
 * The parameters can be modified by subtracting an offset to be defined. For
 * verbatim parameters just pass 0 as @c offset.\n
 * Returns @c NULL on error.
 *
 * @param[in] offset Value to subtract from parameters.
 * @param[in] sigma alphabet.
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
NN_scores*
nn_scores_new_init (float offset, Alphabet* sigma,
                    const char* file, const int line)
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
   if (allocate_init_G_stack (a, u, g, c, offset, this, file, line))
   {
      nn_scores_delete (this);
      return NULL;
   }

   /* create table for mismatch stacking energies */
   if(allocate_init_G_mm_stack_size (a, u, g, c, alphabet_size (sigma),
                                     offset, this,
                                     file, line))
   {
      nn_scores_delete (this);
      return NULL;       
   }

   /* init hairpin loops */
   if (allocate_init_hairpin_loop (offset, this, file, line))
   {
      nn_scores_delete (this);
      return NULL;
   }
   
   /* init hairpin closing base pair energies */
   if (allocate_init_mismatch_hairpin (a, u, g, c,
                                       alphabet_size (sigma),
                                       offset,
                                       this,
                                       file, line))
   {
      nn_scores_delete (this);
      return NULL;
   }

   /* init bulge loops */
   if (allocate_init_bulge_loop (offset, this, file, line))
   {
      nn_scores_delete (this);
      return NULL;
   }

   /* init penalties for non-GC closing base pairs of loops */
   if (allocate_init_non_gc_penalty_for_bp (a, u, g, c, offset, this,
                                            file, line))
   {
      nn_scores_delete (this);
      return NULL;      
   }

   /* init tetra loop index and score values */
   if (allocate_init_tetra_loop (a, u, g, c, offset, this, file, line))
   {
      nn_scores_delete (this);
      return NULL;      
   }

   if (allocate_init_dangle5 (a, u, g, c, alphabet_size (sigma), offset, this,
                              file, line))
   {
      nn_scores_delete (this);
      return NULL; 
   }

   if (allocate_init_dangle3 (a, u, g, c, alphabet_size (sigma), offset, this,
                              file, line))
   {
      nn_scores_delete (this);
      return NULL; 
   }

   if (allocate_init_internal_loop (offset, this, file, line))
   {
      nn_scores_delete (this);
      return NULL;
   }

   if (allocate_init_int11 (a, u, g, c, alphabet_size (sigma), offset, this,
                     file, line))
   {
      nn_scores_delete (this);
      return NULL; 
   }

   if (allocate_init_int21 (a, u, g, c, alphabet_size (sigma), offset, this,
                     file, line))
   {
      nn_scores_delete (this);
      return NULL; 
   }

   if (allocate_init_int22 (a, u, g, c, alphabet_size (sigma), offset, this,
                     file, line))
   {
      nn_scores_delete (this);
      return NULL; 
   }

   if (allocate_init_mismatch_interior (a, u, g, c,
                                        alphabet_size (sigma),
                                        offset,
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


/********************************   Altering   ********************************/
/** @brief Change the parameters in a scoring scheme by "random" values.
 *
 * This function could be used to add thermal noise to the Nearest Neighbour
 * scoring scheme. This means adding tiny bits of random numbers to each
 * parameter in the tables. This is neccessary at least for designing sequences,
 * if the scoring scheme contains equal values for different base
 * combinations.\n
 * To assure reproductibility, the initial seed for the random number generator
 * is a parameter @c seedval to be defined.
 *
 * @params[in] alpha_size size of the alphabet the scheme belongs to.
 * @params[in] seedval seed for the random number generator.
 * @params[in/ out] this scoring scheme to be changed.
 */
void
nn_scores_add_thermal_noise (unsigned long alpha_size,
                             long int seedval,
                             NN_scores* this)
{
   unsigned long i, j, k, l, m, n;
   float rval;

   assert (this);
   assert (this->non_gc_penalty_for_bp);
   assert (this->G_dangle5);
   assert (this->G_stack);

   /* init random number generator */
   srand48 (seedval);

   for (i = 0; i < this->bp_allowed_size; i++)
   {
      /* non_gc_penalty_for_bp */
      rval = (float) drand48 ();
      this->non_gc_penalty_for_bp[i] += (rval - 0.5f) / 100;

      /* G_dangle5, G_dangle3 */
      for (j = 0; j < alpha_size; j++)
      {
         rval = (float) drand48 ();
         this->G_dangle5[i][j] += (rval - 0.5f) / 100;

         rval = (float) drand48 ();
         this->G_dangle3[i][j] += (rval - 0.5f) / 100;
      }

      /* G_stack */
      for (j = 0; j < this->bp_allowed_size; j++)
      {
         rval = (float) drand48 ();
         this->G_stack[i][j] += (rval - 0.5f) / 100;
      }

      /* G_mm_stack */
      for (j = 0; j < this->bp_idx_size; j++)
      {
         rval = (float) drand48 ();
         this->G_mm_stack[i][j] += (rval - 0.5f) / 100;         
      }

      /* G_int11, G_int21 G_int22 */
      for (j = 0; j < this->bp_allowed_size; j++)
      {
         for (k = 0; k < alpha_size; k++)
         {
            for (l = 0; l < alpha_size; l++)
            {
               rval = (float) drand48 ();
               this->G_int11[i][j][k][l] += (rval - 0.5f) / 100;

               for (m = 0; m < alpha_size; m++)
               {
                  rval = (float) drand48 ();
                  this->G_int21[i][j][k][l][m] += (rval - 0.5f) / 100;

                  for (n = 0; n < alpha_size; n++)
                  {
                     rval = (float) drand48 ();
                     this->G_int22[i][j][k][l][m][n] += (rval - 0.5f) / 100;
                  }
               }
            }
         }
      }

      /* G_mismatch_interior, G_mismatch_hairpin */
      for (j = 0; j < alpha_size; j++)
      {
         for (k = 0; k < alpha_size; k++)
         {
            rval = (float) drand48 ();
            this->G_mismatch_interior[i][j][k] += (rval - 0.5f) / 100;

            rval = (float) drand48 ();
            this->G_mismatch_hairpin[i][j][k] += (rval - 0.5f) / 100;
         }
      }
   }

   /* G_tetra_loop */
   for (i = 0; i < this->tetra_loop_size; i++)
   {
      rval = (float) drand48 ();
      this->G_tetra_loop[i] += (rval - 0.5f) / 100;   
   }

   /* G_hairpin_loop */
   for (i = 0; i < this->G_hairpin_loop_size; i++)
   {
      if (this->G_hairpin_loop[i] < FLOAT_UNDEF)
      {
         rval = (float) drand48 ();
         this->G_hairpin_loop[i] += (rval - 0.5f) / 100;  
      } 
   }

   /* G_bulge_loop */
   for (i = 0; i < this->G_bulge_loop_size; i++)
   {
      if (this->G_bulge_loop[i] < FLOAT_UNDEF)
      {
         rval = (float) drand48 ();
         this->G_bulge_loop[i] += (rval - 0.5f) / 100;
      }   
   }

   /* G_internal_loop */
   for (i = 0; i < this->G_internal_loop_size; i++)
   {
      if (this->G_internal_loop[i] < FLOAT_UNDEF)
      {
         rval = (float) drand48 ();
         this->G_internal_loop[i] += (rval - 0.5f) / 100;
      }
   }
}


/*********************************   Access   *********************************/

/** @brief Return size of a tetra loop.
 *
 * Obviously this is 4. But we have to use this function for charma.
 *
 * @params[in] this The scoring scheme.
 */
unsigned long
nn_scores_get_size_tetra_loop (const NN_scores* this)
{
   assert (this);

   return D_TL - 2;             /* remove closing bp */
}

/** @brief Return size of a tetra loop including its closing base pair.
 *
 * Obviously this is 6. But we have to use this function for charma.
 *
 * @params[in] this The scoring scheme.
 */
unsigned long
nn_scores_get_size_tetra_loop_full (const NN_scores* this)
{
   assert (this);

   return D_TL;
}

/** @brief Return no. of tetra loop parameters.
 *
 * Just the no. of parameterised tetra loops.
 *
 * @param[in] this The scoring scheme.
 */
unsigned long
nn_scores_get_no_of_tetra_loops (const NN_scores* this)
{
   assert (this);

   return this->tetra_loop_size;
}

/** Get the sequence of a certain tetra loop.
 *
 * Returns a pointer to the sequence of a tetra loop from the parameter table.
 * Since you got no copy of the sequence, you must NOT free it after use! Be
 * sure to access only indeces smaller than
 * @c nn_scores_get_no_of_tetra_loops(). Please note that the sequence are
 * stored in transformed form.
 *
 * @param[in] i Index of loop to fetch from parameter table.
 * @param[in] this The scoring scheme.
 */
const char*
nn_scores_get_tetra_loop (const unsigned long i, const NN_scores* this)
{
   assert (this);
   assert (i < this->tetra_loop_size);

   return this->tetra_loop[i];
}

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

/** @brief Get penalty for non-GC closing base pairs 
 *
 * @param[in] i 5' pairing partner.
 * @param[in] j 3' pairing partner.
 * @param[in] this Scoring scheme
 */
float
nn_scores_get_G_non_gc_penalty_for_bp (const int i, const int j,
                                       const NN_scores* this)
{
   assert (this);
   assert (this->non_gc_penalty_for_bp);
   assert (this->bp_idx);
   assert (i < sqrtf ((float) this->bp_idx_size));
   assert (j < sqrtf ((float) this->bp_idx_size));

   return this->non_gc_penalty_for_bp[(int)this->bp_idx[i][j]];
}

float
nn_scores_get_G_dangle5 (const int i, const int j, const int im1,
                         const NN_scores* this)
{
   assert (this);
   assert (this->G_dangle5);
   assert (this->bp_idx);
   assert (i < sqrtf ((float) this->bp_idx_size));
   assert (j < sqrtf ((float) this->bp_idx_size));
   assert (im1 < sqrtf ((float) this->bp_idx_size));

   return this->G_dangle5[(int)this->bp_idx[i][j]][im1];
}

float
nn_scores_get_G_dangle3 (const int i, const int j, const int jp1,
                         const NN_scores* this)
{
   assert (this);
   assert (this->G_dangle3);
   assert (this->bp_idx);
   assert (i < sqrtf ((float) this->bp_idx_size));
   assert (j < sqrtf ((float) this->bp_idx_size));
   assert (jp1 < sqrtf ((float) this->bp_idx_size));

   return this->G_dangle3[(int)this->bp_idx[i][j]][jp1];
}

float
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
   float G = 0;

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
float
nn_scores_get_G_stack (const int i, const int j,
                       const int jm1, const int ip1,
                       const NN_scores* scheme)
{
   assert (scheme);
   assert (scheme->G_stack);
   assert (scheme->bp_idx);
   assert (i < sqrtf ((float) scheme->bp_idx_size));
   assert (j < sqrtf ((float) scheme->bp_idx_size));
   assert (ip1 < sqrtf ((float) scheme->bp_idx_size));
   assert (jm1 < sqrtf ((float) scheme->bp_idx_size));
   assert (scheme->bp_idx[i][j] 
           < sqrtf ((float) scheme->G_stack_size));
   assert (  scheme->bp_idx[jm1][ip1]
             < sqrtf ((float) scheme->G_stack_size));
   
   return scheme->G_stack[(int) scheme->bp_idx[i][j]]
                         [(int) scheme->bp_idx[jm1][ip1]];
}

/** @brief Return the mismatch stacking score for a set of bases.
 *
 * @params[in] i i component of the upstream pair of the stack (5' end).
 * @params[in] j j component of the upstream pair of the stack (pairs with i).
 * @params[in] k position j-1.
 * @params[in] l position i+1.
 * @params[in] scheme The scoring scheme.
 */
float
nn_scores_get_G_mm_stack (const char i, const char j,
                          const char k, const char l,
                          const NN_scores* scheme)
{
   assert (scheme);
   assert (scheme->G_mm_stack);
   assert (scheme->bp_idx);
   assert (i < sqrtf ((float) scheme->bp_idx_size));
   assert (j < sqrtf ((float) scheme->bp_idx_size));
   assert (k < sqrtf ((float) scheme->bp_idx_size));
   assert (l < sqrtf ((float) scheme->bp_idx_size));
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
 * @param[in] this scoring sceme.
 */
float
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

/** @brief Returns the penalty for the closing basepair of a hairpin.
 *
 * @param[in] i 5' base of closing pair.
 * @param[in] j 3' base of closing pair.
 * @param[in] ip1 i+1 base.
 * @param[in] jm1 j-1 base.
 * @param[in] size Loop size.
 * @param[in] this Scores.
 */
float
nn_scores_get_G_hairpin_mismatch (const int i,
                                  const int j,
                                  const int ip1,
                                  const int jm1,
                                  const unsigned long size,
                                  const NN_scores* this)
{

   /* mismatch penalty for the mismatch interior to the closing basepair of
      the hairpin. triloops get non-parameterised mismatch penalty */
   if (size == D_MM_H)
   {
      return this->non_gc_penalty_for_bp[(int)this->bp_idx[i][j]];
   }
   return this->G_mismatch_hairpin[(int)this->bp_idx[i][j]][ip1][jm1];
}

/** @brief Returns the score for a hairpin loop of certain size.
 *
 * @params[in] i 5' base of closing pair.
 * @params[in] j 3' base of closing pair.
 * @params[in] size Length of the loop (unpaired bases only).
 * @params[in] this Scoring scheme.
 */
float
nn_scores_get_G_hairpin_loop (const char* seq,
                              const unsigned long i,
                              const unsigned long j,
                              const unsigned long size,
                              const NN_scores* this)
{
   float G = 0;

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
         + /*(int)*/ (NN_LXC37 
                  * logf((float) size / (this->G_hairpin_loop_size - 1)));
   }

   G += nn_scores_get_G_hairpin_mismatch (seq[i], seq[j],
                                          seq[i + 1], seq[j - 1],
                                          size,
                                          this);

   /* tetraloop bonus */
   if (size == nn_scores_get_size_tetra_loop(this))
   {
      G += nn_scores_get_G_tetra_loop (seq, i, this);
   }
   /* mfprintf (stderr, "G= %d\n", G); */
   return G;
}

/** @brief Energy for both closing base pairs of a bulge loop.
 *
 * @param[in] i1 Partner i of the first pair.
 * @param[in] j1 Partner j of the first pair.
 * @param[in] i2 Partner i of the 2nd pair.
 * @param[in] j2 Partner j of the 2nd pair.
 * @param[in] size Size of the loop.
 * @param[in] this Scoring scheme.
 */
float
nn_scores_get_G_bulge_stack (const int bi1, const int bj1,
                             const int bj2, const int bi2,
                             const unsigned long size,
                             const NN_scores* this)
{
   float G = 0;
 
   assert (this);
   assert (this->non_gc_penalty_for_bp);

   if (size == 1)
   {
      G += nn_scores_get_G_stack (bi1, bj1, bj2, bi2, this);
   }
   else
   {
      /* bulge loops larger than 1 get penalty term for non-gc closing
         basepairs */
      G += this->non_gc_penalty_for_bp[(int)this->bp_idx[bi1][bj1]];
      G += this->non_gc_penalty_for_bp[(int)this->bp_idx[bj2][bi2]];
   }

   return G;
}

/** @brief Returns the score for a bulge loop.
 *
 * @params[in] i 5' base of closing pair.
 * @params[in] j 3' base of closing pair.
 * @params[in] size Length of the loop (unpaired bases only).
 * @params[in] this Scoring scheme.
 */
float
nn_scores_get_G_bulge_loop (const int bi1, const int bj1,
                            const int bi2, const int bj2,
                            const unsigned long size,
                            const NN_scores* this)
{
   float G = 0;

   assert (this);
   assert (this->G_bulge_loop);
   assert (this->non_gc_penalty_for_bp);
   assert (this->bp_idx);
   assert (bi1 < sqrtf ((float) this->bp_idx_size));
   assert (bj1 < sqrtf ((float) this->bp_idx_size));
   assert (bi2 < sqrtf ((float) this->bp_idx_size));
   assert (bj2 < sqrtf ((float) this->bp_idx_size));
   assert (  (unsigned) this->bp_idx[bi1][bj1] 
           < this->bp_allowed_size);
   assert (  (unsigned) this->bp_idx[bj2][bi2] 
           < this->bp_allowed_size);

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

   G += nn_scores_get_G_bulge_stack (bi1, bj1, bj2, bi2, size, this);

   return G;
}

/* note: everything given 5' -> 3' direction */
float
nn_scores_get_G_internal_2x2_loop (const int bi1,
                                   const int bj1,
                                   const int bi1p1,
                                   const int bi2m1,
                                   const int bj2,
                                   const int bi2,
                                   const int bj2p1,
                                   const int bj1m1,
                                   const NN_scores* this)
{
   assert (this);
   assert (this->bp_idx);
   assert (this->G_int22);
   assert (bi1 < sqrtf ((float) this->bp_idx_size));
   assert (bj1 < sqrtf ((float) this->bp_idx_size));
   assert (bi2 < sqrtf ((float) this->bp_idx_size));
   assert (bj2 < sqrtf ((float) this->bp_idx_size));
   assert ((unsigned)(this->bp_idx[bi1][bj1]
                      * this->bp_idx[bj2][bi2] 
                      * bi1p1
                      * bi2m1
                      * bj2p1
                      * bj1m1)
                      < this->G_int22_size);

   return this->G_int22[(int)this->bp_idx[bi1][bj1]] /* bp 1 */
                       [(int)this->bp_idx[bj2][bi2]] /* bp 2 */
                       [bi1p1][bi2m1][bj2p1][bj1m1];   /* unpaired bases */
}

float
nn_scores_get_G_internal_1x2_loop (const int bi1,
                                   const int bj1,
                                   const int bi1p1,
                                   const int bj2p1,
                                   const int bj1m1,
                                   const int bj2,
                                   const int bi2,
                                   const NN_scores* this)
{
   assert (this);
   assert (this->bp_idx);
   assert (this->G_int21);
   assert (bi1 < sqrtf ((float) this->bp_idx_size));
   assert (bj1 < sqrtf ((float) this->bp_idx_size));
   assert (bi2 < sqrtf ((float) this->bp_idx_size));
   assert (bj2 < sqrtf ((float) this->bp_idx_size));
   assert ((unsigned)(this->bp_idx[bi1][bj1]
                      * this->bp_idx[bj2][bi2]
                      * bi1p1
                      * bj2p1
                      * bj1m1)
           < this->G_int21_size);

   return this->G_int21[(int)this->bp_idx[bi1][bj1]] /* bp 1 */
                       [(int)this->bp_idx[bj2][bi2]] /* bp 2 */
                       [bi1p1][bj2p1][bj1m1];         /* unpaired bases */
}

float
nn_scores_get_G_internal_1x1_loop (const int bi1,
                                   const int bj1,
                                   const int bi1p1,
                                   const int bj1m1,
                                   const int bi2,
                                   const int bj2,
                                   const NN_scores* this)
{
   assert (this);
   assert (this->bp_idx);
   assert (this->G_int11);
   assert (bi1 < sqrtf ((float) this->bp_idx_size));
   assert (bj1 < sqrtf ((float) this->bp_idx_size));
   assert (bi2 < sqrtf ((float) this->bp_idx_size));
   assert (bj2 < sqrtf ((float) this->bp_idx_size));
   assert ((unsigned)(this->bp_idx[bi1][bj1] * this->bp_idx[bj2][bi2] * bi1p1 * bj1m1)
           < this->G_int11_size);

   return this->G_int11[(int)this->bp_idx[bi1][bj1]] /* bp 1 */
                       [(int)this->bp_idx[bj2][bi2]] /* bp 2 */
                       [bi1p1][bj1m1];               /* unpaired bases */
}

float
nn_scores_get_G_mismatch_interior (const int i,
                                   const int j,
                                   const int ip,
                                   const int jm,
                                   const NN_scores* this)
{
   assert (this);
   assert (this->bp_idx);
   assert (this->G_mismatch_interior);
   assert (i < sqrtf ((float) this->bp_idx_size));
   assert (j < sqrtf ((float) this->bp_idx_size));
   assert ((unsigned)(this->bp_idx[i][j] * ip * jm)
           < this->G_mismatch_interior_size);

   return this->G_mismatch_interior[(int)this->bp_idx[i][j]][ip][jm];
}

float
nn_scores_get_G_internal_loop (const char* seq,
                               const unsigned long size1,
                               const unsigned long size2,
                               const unsigned long pi1,
                               const unsigned long pj1,
                               const unsigned long pi2,
                               const unsigned long pj2,
                               const NN_scores* this)
{
   float G = 0;
   int bp1, bp2;
   int bi1p, bi2m, bj2p, bj1m;  /* bi1p = seq[pi1 + 1] */
   unsigned long size;

   assert (seq);
   assert (this);
   assert (this->bp_idx);
   assert (this->G_int11);
   assert (this->G_int21);
   assert (this->G_int22);
   assert (this->G_mismatch_interior);
   assert (pi1 < pj1);
   assert (pi1 < pi2);
   assert (pi2 < pj2);
   assert (pj2 < pj1);

   bp1 =  (int)this->bp_idx[(int)seq[pi1]][(int)seq[pj1]];
   bp2 =  (int)this->bp_idx[(int)seq[pj2]][(int)seq[pi2]];
   bi1p = (int)seq[pi1 + 1];
   bi2m = (int)seq[pi2 - 1];
   bj2p = (int)seq[pj2 + 1];
   bj1m = (int)seq[pj1 - 1];

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
      /* note switched order of bp1 and bp2 compared to 1x2 loop */
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
            + (int) (NN_LXC37 
                     * logf ((float) size / (this->G_internal_loop_size - 1)));
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
   string = XMALLOC (sizeof (*string) * ((scheme->bp_allowed_size * 3) + 1));
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
   rprec = floorf (log10f ((float) rprec) + 1);
   
   /* calculate linewidth */
   pline_width = rprec;
   pline_width += 3;            /* \s|\s */
   pline_width *= alpha_size;
   pline_width += 2;            /* N\n */

   string = XMALLOC (sizeof (*string) *
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
      msprintf (string, " | %*c", rprec, alphabet_no_2_base ((char) i, sigma));
      string += (rprec + 3);      
   }
   string[0] = '\n';
   string++;

   /* print index lines */
   for (i = 0; i < alpha_size; i++)
   {
      msprintf (string, "%c", alphabet_no_2_base ((char) i, sigma));
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
   float tmp;
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
            rprec += floorf (log10f ((float) tmp) + 1);
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
   string = XMALLOC (sizeof (*string) * ((pline_width * (matrix_edge+ 1)) + 1));
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

            msprintf (string, "%*.0f", rprec,
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
   float tmp;
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
            rprec += floorf (log10f ((float) tmp) + 1);
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
   string = XMALLOC (sizeof (*string) * ((pline_width * (matrix_rows+ 1)) + 1));
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
                   alphabet_no_2_base ((char) i, sigma),
                   alphabet_no_2_base ((char) j, sigma));
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
            
            msprintf (string, "%*.0f", rprec,
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
   int rprec, rprec_idx, pline_width = 1;
   float tmp;
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

      if (tmp < FLOAT_UNDEF)
      {
         if (tmp < 0)
         {
            tmp *= (-1);
            rprec++;
         }
         
         if (tmp > 0)
         {
            rprec += floorf (log10f ((float) tmp) + 1);
         }
      }

      if (rprec > pline_width) 
      {
         pline_width = rprec;
      }       
   }
   rprec = pline_width;
   rprec_idx = floorf (log10f ((float) scheme->G_hairpin_loop_size) + 1);

   /* add up components of a line */
   pline_width += 3; /*: \n*/
   pline_width += rprec_idx;

   /* allocate memory for the string and the undef symbol
      therefore we add "rprec + 1" which is the 2 (for null terminating the
      strings) */
   string = XMALLOC (sizeof (*string) *
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

      if (scheme->G_hairpin_loop[i] < FLOAT_UNDEF)
      {
         msprintf (string, "%*.0f", rprec, scheme->G_hairpin_loop[i]);  
      }
      else
      {
         msprintf (string, "%s", en_undef);
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
   float tmp;
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
               rprec += floorf (log10f (tmp) + 1);
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
   string = (char*) XMALLOC (sizeof (*string) * pline_width);
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
      msprintf (string, "  | %*c", rprec, alphabet_no_2_base ((char) i, sigma));
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
   string = (char*) XMALLOC (sizeof (*string) * ((pline_width * z) + 1));
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
         msprintf (string, "%c", alphabet_no_2_base ((char) j, sigma));
         string++;
         for (k = 0; k < x; k++)
         {
            msprintf (string, " | %*.0f", rprec,
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
   int rprec, rprec_idx, pline_width = 1;
   float tmp;
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

      if (tmp < FLOAT_UNDEF)
      {
         if (tmp < 0)
         {
            tmp *= (-1);
            rprec++;
         }

         if (tmp > 0)
         {
            rprec += floorf (log10f (tmp) + 1);
         }
      }
      
      if (rprec > pline_width) 
      {
         pline_width = rprec;
      }       
   }
   rprec = pline_width;
   rprec_idx = floorf (log10f ((float) scheme->G_bulge_loop_size) + 1);

   /* add up components of a line */
   pline_width += 3; /*: \n*/
   pline_width += rprec_idx;

   /* allocate memory for the string and the undef symbol
      therefore we add "rprec + 1" which is the 2 (for null terminating the
      strings) */
   string = XMALLOC (sizeof (*string) *
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

      if (scheme->G_bulge_loop[i] < FLOAT_UNDEF)
      {
         msprintf (string, "%*.0f", rprec, scheme->G_bulge_loop[i]);
      }
      else
      {
         msprintf (string, "%s", en_undef);         
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
   float tmp;
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
         rprec += floorf (log10f ((float) tmp) + 1);
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

   string = (char*) XMALLOC (sizeof (*string)
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

      msprintf (string, "%*.0f", rprec, scheme->non_gc_penalty_for_bp[i]);
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
   float tmp;
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
         rprec += floorf (log10f ((float) tmp) + 1);
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
   string = (char*) XMALLOC (sizeof (*string)
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

      msprintf (string, ": %*.0f", rprec, scheme->G_tetra_loop[i]);
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
   float tmp;
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
            rprec += floorf (log10f ((float) tmp) + 1);
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
   string = XMALLOC (sizeof (*string) * ((pline_width * (rows + 1)) + 1));
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
      msprintf (string, " | %*c", rprec, alphabet_no_2_base ((char) i, sigma));
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
         
         msprintf (string, "%*.0f", rprec,
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
   float tmp;
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
            rprec += floorf (log10f ((float) tmp) + 1);
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
   string = XMALLOC (sizeof (*string) * ((pline_width * (rows + 1)) + 1));
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
      msprintf (string, " | %*c", rprec, alphabet_no_2_base ((char) i, sigma));
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
         
         msprintf (string, "%*.0f", rprec,
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
                     float** matrix)
{
   int tmp = 0;
   int cval;
   unsigned long i, j;

   /* find largest non-negative no. */
   for (i = 0; i < rows; i++)
   {
      for (j = 0; j < cols; j++)
      {
         if (matrix[i][j] > 0)
         {
            cval = floorf (log10f (matrix[i][j]) + 1);
         }
         else if (matrix[i][j] < 0)
         {
            cval = floorf (log10f (matrix[i][j] * (-1)) + 2);
         }
         else
         {
            cval = 0;
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
   int rprec, rprec_idx, pline_width = 1;
   float tmp;
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

      if (tmp < FLOAT_UNDEF)
      {
         if (tmp < 0)
         {
            tmp *= (-1);
            rprec++;
         }
         
         if (tmp > 0)
         {
            rprec += floorf (log10f (tmp) + 1);
         }
      }
      
      if (rprec > pline_width) 
      {
         pline_width = rprec;
      }       
   }
   rprec = pline_width;
   rprec_idx = floorf (log10f ((float) scheme->G_internal_loop_size) + 1);

   /* add up components of a line */
   pline_width += 3; /*: \n*/
   pline_width += rprec_idx;

   /* allocate memory for the string and the undef symbol
      therefore we add "rprec + 1" which is the 2 (for null terminating the
      strings) */
   string = XMALLOC (sizeof (*string) *
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

      if (scheme->G_internal_loop[i] < FLOAT_UNDEF)
      {
         msprintf (string, "%*.0f", rprec, scheme->G_internal_loop[i]);  
      }
      else
      {
         msprintf (string, "%s", en_undef);
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

   header = XMALLOC (sizeof (*header) * (hsize + 1));
   if (header == NULL)
   {
      return;
   }
   string_start = header;

   msprintf (header, "     ");
   header += 5;
   for (i = 0; i < asize; i++)
   {
      msprintf (header, " | %*c", rprec, alphabet_no_2_base((char) i, sigma));
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
   string = XMALLOC (sizeof (*string) * (store_size + 1));
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
            msprintf (string, "    %c", alphabet_no_2_base((char) k, sigma));
            string += 5;

            for (l = 0; l < asize; l++)
            {
               msprintf (string, " | %*.0f", rprec, scheme->G_int11[i][j][k][l]);
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

   header = XMALLOC (sizeof (*header) * (hsize + 1));
   if (header == NULL)
   {
      return;
   }
   string_start = header;

   msprintf (header, "      ");
   header += 6;
   for (i = 0; i < asize; i++)
   {
      msprintf (header, " | %*c", rprec, alphabet_no_2_base((char) i, sigma));
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
   string = XMALLOC (sizeof (*string) * (store_size + 1));
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
            msprintf (string, "    %c:\n", alphabet_no_2_base((char) k, sigma));
            string += 7;     
            msprintf (string, "%s", header);
            string += hsize;
         
            for (l = 0; l < asize; l++)
            {
               msprintf (string, "     %c", alphabet_no_2_base((char) l,
                                                               sigma));
               string += 6;
               
            for (m = 0; m < asize; m++)
            {
               msprintf (string, " | %*.0f",
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

   header = XMALLOC (sizeof (*header) * (hsize + 1));
   if (header == NULL)
   {
      return;
   }
   string_start = header;

   msprintf (header, "       ");
   header += 7;
   for (i = 0; i < asize; i++)
   {
      msprintf (header, " | %*c", rprec, alphabet_no_2_base((char) i, sigma));
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
   string = XMALLOC (sizeof (*string) * (store_size + 1));
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
            msprintf (string, "    %c:\n", alphabet_no_2_base((char) k, sigma));
            string += 7;

            for (l = 0; l < asize; l++)
            {
               msprintf (string, "     %c:\n", alphabet_no_2_base((char) l,
                                                                  sigma));
               string += 8;
               msprintf (string, "%s", header);
               string += hsize;
         
               for (m = 0; m < asize; m++)
               {
                  msprintf (string, "      %c", alphabet_no_2_base((char) m,
                                                                   sigma));
                  string += 7;
                  
                  for (n = 0; n < asize; n++)
                  {
                     msprintf (string, " | %*.0f",
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

/** @brief Print the mismatch interior energies of a scoring scheme to a stream.
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
   float tmp;
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
               rprec += floorf (log10f ((float) tmp) + 1);
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
   string = (char*) XMALLOC (sizeof (*string) * pline_width);
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
      msprintf (string, "  | %*c", rprec, alphabet_no_2_base ((char) i, sigma));
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
   string = (char*) XMALLOC (sizeof (*string) * ((pline_width * z) + 1));
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
         msprintf (string, "%c", alphabet_no_2_base ((char) j, sigma));
         string++;
         for (k = 0; k < x; k++)
         {
            msprintf (string, " | %*.0f", rprec,
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
