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
 *  @file libcrbapps/fold.c
 *
 *  @brief RNA Secondary Structure Prediction Tool
 *
 *  Module: fold
 *
 *  Library: libcrbapps
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-06-17
 *
 *
 *  Revision History:
 *         - 2008Jun17 bienert: created
 *
 */


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <libcrbbasic/crbbasic.h>
#include <libcrbrna/crbrna.h>
#include "fold_cmdline.h"
#include "fold.h"


typedef struct {
      float v;               /* value */
      char from;             /* previous node */
      ArrayUlong k;          /* k if no previous node */
} struct_cell;

enum directions{
   L = 1,                       /* from left */
   B = 2,                       /* from below */
   D = 4,                       /* diagonal (pair) */
   K = 8,
};

static int 
fold_cmdline_parser_postprocess (const struct fold_args_info* args_info,
                                 unsigned long *seqlen, Alphabet* sigma)
{
   unsigned long i;
   *seqlen = 0;

   /* check input sequence */
   if (args_info->inputs_num == 1)
   {
      THROW_ERROR_MSG ("RNA sequence required as argument, try "
                       "`%s --help` for more information.", get_progname());
      return 1;
   }
   if (args_info->inputs_num != 2)
   {
      THROW_ERROR_MSG ("Only one RNA sequence allowed as argument, try "
                       "`%s --help` for more information.", get_progname());
      return 1;
   }

   /* verify sequence */
   *seqlen = strlen (args_info->inputs[1]);
   for (i = 0; i < *seqlen; i++)
   {
       args_info->inputs[1][i] =
          alphabet_base_2_no (args_info->inputs[1][i], sigma);

      if (args_info->inputs[1][i] == CHAR_UNDEF)
      {
         return 1;
      }
   }

   /* check min length parameter */
   if (args_info->loop_length_given)
   {
      if (args_info->loop_length_arg < 1)
      {
         THROW_ERROR_MSG ("Option \"--loop-length\" requires an integer larger "
                          "than 1 as argument, found: %ld",
                          args_info->loop_length_arg);
         return 1; 
      }
   }

   return 0;
}

static int
calc_cell_nussinov(const unsigned long i,
                   const unsigned long j,
                   const char* sequence,
                   struct_cell** matrix,
                   float** scores)
{
   unsigned long k;

   if (matrix[i][j - 1].from != 0)
   {
      matrix[i][j].v = matrix[i][j - 1].v;
      matrix[i][j].from = L;
   }

   if (matrix[i + 1][j].from != 0)
   {
      if (fabs (matrix[i][j].v - matrix[i + 1][j].v) < 0.0001)
      {
         matrix[i][j].from |= B;
      }
      else if (matrix[i][j].v > matrix[i + 1][j].v)
      {
         matrix[i][j].v = matrix[i + 1][j].v;
         matrix[i][j].from = B;
      }
   }

   if ((scores[(int)sequence[i]][(int)sequence[j]] > 0)
       ||(scores[(int)sequence[i]][(int)sequence[j]] < 0))
   {
      /*mprintf ("(%lu,%lu) %.3f - %.3f + %3.f\n", i, j, matrix[i][j].v,
                                                matrix[i + 1][j -1 ].v, 
                                                scores[(int)sequence[i]]
                                                [(int)sequence[j]]);*/
      if (fabsf (matrix[i][j].v - (  matrix[i + 1][j -1].v 
                                     + scores[(int)sequence[i]]
                                     [(int)sequence[j]])) < 0.0001)
      {
         if (matrix[i + 1][j - 1].from != 0)
         {
            matrix[i][j].from |= D;
         }
      }
      else if (matrix[i][j].v > (matrix[i + 1][j -1].v
                                 + scores[(int)sequence[i]][(int)sequence[j]]))
      {
         matrix[i][j].v = (matrix[i + 1][j - 1].v
                           + scores[(int)sequence[i]][(int)sequence[j]]);
         matrix[i][j].from = D;
      }
   }

   for (k = i + 1; k < j; k++)
   {
      if ((matrix[i][k - 1].from != 0) && (matrix[k][j].from != 0))
      {
         if (fabsf (matrix[i][j].v - (matrix[i][k - 1].v + matrix[k][j].v)) 
             < 0.0001)
         {
            ARRAY_ULONG_PUSH (matrix[i][j].k, k,
                              {
                                 return 1;
                              });
            matrix[i][j].from |= K;
            /*if (! matrix[i][j].from & K)
            {
               matrix[i][j].k = k;
               matrix[i][j].from |= K;
            }*/
         }
         else if (matrix[i][j].v > (matrix[i][k - 1].v + matrix[k][j].v))
         {
            /* mprintf ("K"); */
            ARRAY_RESET (matrix[i][j].k);
            ARRAY_ULONG_PUSH (matrix[i][j].k, k,
                              {
                                 return 1;
                              });           
            
            matrix[i][j].v = matrix[i][k - 1].v + matrix[k][j].v;
            matrix[i][j].from = K;
            /* matrix[i][j].k = k; */
         }
      }
   }
   /* mprintf ("|"); */
   /* 2: 8384 */
   /* 3: 5468 */
   /* 4:  873 */
   /* 5:  106 */
   /* 6:   19 */
   /* 7:    5 */
   /* 8:    3 */
   /* mprintf (" (%lu, %lu): %.3f\n", i, j, matrix[i][j].v); */

   return 0;
}

static int
calc_matrix_nussinov (const char* sequence,
                      struct_cell** matrix,
                      const unsigned long seqlen,
                      const unsigned long lmin,
                      float** scores)
{
   unsigned long l, i, j;

   assert (sequence);
   assert (matrix);

   for (l = 0; l < seqlen; l++)
   {
      for (i = 0; i < seqlen - l; i++)
      {
         j = i + l;

         if (j - i < (lmin + 1))
         {
            /* mprintf ("<l: (i,j) : (%lu,%lu)\n", i, j); */
            matrix[i][j].v = 0;
            matrix[i][j].from = 0;
         }
         else
         {
            /* mprintf (">l: (i,j) : (%lu,%lu)\n", i, j); */
            if(calc_cell_nussinov(i, j, sequence, matrix, scores))
            {
               return 1;
            }
         }
      }
   }

   /*for (i = 0; i < seqlen; i++)
   {
      for (j = 0; j < seqlen; j++)
      {
         mprintf ("%.3f ", matrix[i][j].v);
      }
      mprintf ("\n");
   }
   mprintf ("\n");
   for (i = 0; i < seqlen; i++)
   {
      mprintf (" %c |", transform_number_2_base (sequence[i]));
   }
   mprintf ("\n");
   for (i = 0; i < seqlen; i++)
   {
      for (j = 0; j < seqlen; j++)
      {
         if (matrix[i][j].from != 0)
         {
            if (matrix[i][j].from & L)
            {
               mprintf ("L");               
            }
            else
            {
               mprintf (" ");
            }
            if (matrix[i][j].from & B)
            {
               mprintf ("B");               
            }
            else
            {
               mprintf (" ");
            }
            if (matrix[i][j].from & D)
            {
               mprintf ("D");               
            }
            else
            {
               mprintf (" ");
            }
            mprintf ("|");            
         }
         else
         {
            mprintf ("xxx|");
         }
      }
      mprintf (" %c\n", transform_number_2_base (sequence[i]));
      }

   mprintf ("\n");
   for (i = 0; i < seqlen; i++)
   {
      mprintf ("%2lu ", i);
   }
   mprintf ("\n");
   for (i = 0; i < seqlen; i++)
   {
      mprintf (" %c ", transform_number_2_base (sequence[i]));
   }
   mprintf ("\n");
   for (i = 0; i < seqlen; i++)
   {
      for (j = 0; j < seqlen; j++)
      {
         if (matrix[i][j].k != ULONG_UNDEF)
         {
            mprintf ("%2lu ", matrix[i][j].k);
         }
         else
         {
            mprintf ("xx ");
         }
      }
      mprintf (" %c %2lu\n", transform_number_2_base (sequence[i]), i);
      }*/

   return 0;
}

static void
traceback_matrix (const unsigned long i, const unsigned long j,
                  Str* structure, struct_cell** matrix)
{
   if (i < j)
   {
      /* follow diagonals, first */
      if (matrix[i][j].from & B)
      {
         /* mprintf ("B "); */
         traceback_matrix (i + 1, j, structure, matrix);
      }
      else if (matrix[i][j].from & L)
      {
         /* mprintf ("L "); */
         traceback_matrix (i, j - 1, structure, matrix);
      }
      else if (matrix[i][j].from & D)
      {
         /* mprintf ("D "); */
         /* mfprintf (stderr, "PAIR: %lu,%lu\n", i, j); */
         str_at (structure, i, '(');
         str_at (structure, j, ')');
         traceback_matrix (i + 1, j - 1, structure, matrix);
      }
      else if (ARRAY_CURRENT (matrix[i][j].k) != 0)
      {
         /*mprintf ("K");*/
         /* mprintf ("K: %lu ", ARRAY_SIZE (matrix[i][j].k)); */
         traceback_matrix (i, ARRAY_ACCESS(matrix[i][j].k, 0) - 1,
                           structure,
                           matrix); 
         traceback_matrix (ARRAY_ACCESS(matrix[i][j].k, 0), j,
                           structure,
                           matrix);         
      }

      /*  2: 3611 */
      /*  3:  290 */
      /*  4:   22 */
      /*  5:   10 */
      /*  6:    5 */
      /*  7:    5 */
      /*  8:    2 */
      /*  9:    2 */
      /* 10:    3 */
      /* 12:    1 */
      /* 16:    2 */
      /* 17:    1 */
      /* 19:    1 */
      /* 23:    1 */
      /* 25:    1 */
      /* 30:    1 */
   }
}

/* bienert: uncommented 27/04/2009, seems not to be used anymore */
/* static void */
/* traceback_nussinov (const unsigned long i,  */
/*                     const unsigned long j, */
/*                     Str* structure, */
/*                     const char* sequence, */
/*                     struct_cell** matrix, */
/*                     float** scores) */
/* { */
/*    unsigned long k; */

/*    if (i < j) */
/*    { */
/*       /\*mprintf ("(%c%lu,%c%lu) %.4f ", transform_number_2_base(sequence[i]), */
/*         i, transform_number_2_base(sequence[j]), j, matrix[i][j].v); */
/*         mprintf ("b:%.3f ", matrix[i + 1][j].v);*\/ */
/*       if (fabsf(matrix[i][j].v - matrix[i + 1][j].v) < 0.0001) */
/*       { */
/*          /\*mprintf ("B ");*\/ */
/*          traceback_nussinov (i + 1, j, structure, sequence, matrix, scores); */
/*       } */
/*       else */
/*       { */
/*          /\*mprintf ("l:%.3f ", matrix[i][j - 1].v);*\/ */
/*          if (fabsf(matrix[i][j].v - matrix[i][j - 1].v) < 0.0001) */
/*          { */
/*             /\*mprintf ("L ");*\/ */
/*             traceback_nussinov (i, j - 1, structure, sequence, matrix, scores); */
/*          } */
/*          else */
/*          { */
/* /\*         mprintf ("B: %f == %f + %f\n", matrix[i][j].v, matrix[i+1][j-1].v, *\/ */
/* /\*                      scores[(int)sequence[i]][(int)sequence[j]]); *\/ */
/*             /\*mprintf ("d:%.4f ", matrix[i + 1][j - 1].v);*\/ */
/*             if (fabsf(matrix[i][j].v - (matrix[i + 1][j - 1].v  */
/*                                   + scores[(int)sequence[i]][(int)sequence[j]])) */
/*                 < 0.0001) */
/*             { */
/*                /\* mfprintf (stderr, "PAIR: %lu,%lu\n", i, j); *\/ */
/*                /\*mprintf ("D ");*\/ */
/*                str_at (structure, i, '('); */
/*                str_at (structure, j, ')'); */
/*                traceback_nussinov (i + 1, j - 1, */
/*                                    structure, */
/*                                    sequence, */
/*                                    matrix, */
/*                                    scores); */
/*             } */
/*             else */
/*             { */
/*                for (k = i + 1; k < j; k++) */
/*                { */
/*                   if (fabsf(  matrix[i][j].v  */
/*                            - (matrix[i][k - 1].v + matrix[k][j].v)) < 0.0001) */
/*                   { */
/*                      /\*mprintf ("K%lu ", k);*\/ */
/*                      traceback_nussinov (i, k - 1, */
/*                                          structure, */
/*                                          sequence, */
/*                                          matrix, */
/*                                          scores);  */
/*                      traceback_nussinov (k, j, */
/*                                          structure, */
/*                                          sequence, */
/*                                          matrix, */
/*                                          scores); */
/*                      return; */
/*                   } */
/*                } */
/*                /\*mprintf("N ");*\/ */
/*             } */
/*          } */
/*       } */
/*    } */
/* } */

static Str*
pred_2D_structure_nussinov (const unsigned long l,
                            const char* sequence,
                            const unsigned long seqlen, Alphabet* sigma)
{
   struct_cell** matrix = NULL;
   float** scores = NULL;
   unsigned long i, j;
   int error = 0;
   Str* structure = STR_NEW_CHAR ('.', seqlen);

   assert (sequence);

   if (structure == NULL)
   {
      return structure;
   }

   /* allocate matrix */
   matrix = (struct_cell**) XMALLOC_2D (seqlen,
                                 seqlen,
                                 sizeof (*matrix[0]));
   if (matrix == NULL)
   {
      error = 1;
   }

   for (i = 0; (i < seqlen); i++)
   {
      for (j = 0; j < seqlen; j++)
      {
         matrix[i][j].v = 0.0f;
         matrix[i][j].from = 0;
         ARRAY_ULONG_INIT (matrix[i][j].k, 2);
         if (ARRAY_IS_NULL (matrix[i][j].k))
         {
            i = seqlen;
            j = seqlen;
            error = 1;      
         }
         /* matrix[i][j].k = ULONG_UNDEF; */
      }
   }

   /* fill matrix */
   if (! error)
   {
      scores = create_scoring_matrix (sigma);
      if (scores == NULL)
      {
         error = 1;
      }
   }

   if (! error)
   {
      error = calc_matrix_nussinov (sequence, matrix, seqlen, l, scores);
   }

   /* traceback */
   if ((! error) && (seqlen > 0))
   {
      mprintf ("> MFE: %.3f\n", matrix[0][seqlen - 1].v);
    /*traceback_nussinov (0, seqlen - 1, structure, sequence, matrix, scores);*/
      traceback_matrix (0, seqlen - 1, structure, matrix);
      /* mprintf ("\n"); */
   }

   
   for (i = 0; (i < seqlen); i++)
   {
      for (j = 0; j < seqlen; j++)
      {
         ARRAY_DELETE (matrix[i][j].k);
      }
   }
   XFREE_2D ((void**)matrix);
   XFREE_2D ((void**)scores);

   if (error)
   {
      str_delete (structure);
      return NULL;
   }

   return structure;
}

int
fold_main(const char *cmdline)
{
   struct fold_args_info fold_args;
   Str* structure = NULL;
   unsigned long seqlen, i;
   Alphabet* sigma = NULL;
   int retval = 0;

   /* command line parsing */
   fold_cmdline_parser_init (&fold_args);

   retval = fold_cmdline_parser_string (cmdline, &fold_args, get_progname());

   if (retval == 0)
   {
      retval = fold_cmdline_parser_required ();
   }

   /* postprocess arguments */
   if (retval == 0)
   {
      sigma = ALPHABET_NEW_PAIR ("AUGC", "augc", 4);
      retval = fold_cmdline_parser_postprocess (&fold_args, &seqlen, sigma);
   }

   /* do work */
   if (retval == 0)
   {
      structure = pred_2D_structure_nussinov (fold_args.loop_length_arg,
                                              fold_args.inputs[1],
                                              seqlen,
                                              sigma);
      if (structure == NULL)
      {
         retval = 1;
      }
      else
      {
         for (i = 0; i < seqlen; i++)
         {
            mprintf ("%c", alphabet_no_2_base (fold_args.inputs[1][i], sigma));
         }
         mprintf ("\n");
         mprintf ("%s\n", str_get (structure));
      }
   }

   /* finalise */
   str_delete (structure);
   alphabet_delete (sigma);
   fold_cmdline_parser_free (&fold_args);

   if (retval == 0)
   {
      return EXIT_SUCCESS;
   }
   else
   {
      return EXIT_FAILURE;
   }
}
