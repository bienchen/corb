/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
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
#include "alphabet.h"
#include "fold_cmdline.h"
#include "fold.h"

typedef struct {
      float v;                  /* value */
      char from;                /* previous node */
      unsigned long k;          /* k if no previous node */
} struct_cell;

enum directions{
   L = 1,                       /* from left */
   B = 2,                       /* from below */
   D = 4,                       /* diagonal (pair) */
};

static int 
fold_cmdline_parser_postprocess (const struct fold_args_info* args_info)
{
   unsigned long i;
   unsigned long seqlen = 0;

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
   seqlen = strlen (args_info->inputs[1]);
   for (i = 0; i < seqlen; i++)
   {
       args_info->inputs[1][i] =
          transform_base_2_number (args_info->inputs[1][i]);

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

void
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
      if (fabs (matrix[i][j].v - matrix[i + 1][j].v) < 0.000001)
      {
         matrix[i][j].from |= B;      
      }
      else if (matrix[i][j].v > matrix[i + 1][j].v)
      {
         matrix[i][j].v = matrix[i + 1][j].v;
         matrix[i][j].from = B;
      }
   }

   if (scores[(int)sequence[i]][(int)sequence[j]] != 0)
   {
      if (fabsf (matrix[i][j].v - (  matrix[i + 1][j -1].v 
                                     + scores[(int)sequence[i]]
                                     [(int)sequence[j]])) < 0.000001)
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
      if ((matrix[i][k - 1].from != 0) && (matrix[i][k - 1].from != 0))
      {
         if (fabsf (matrix[i][j].v - (matrix[i][k - 1].v + matrix[k][j].v)) 
             < 0.000001)
         {
            matrix[i][j].k = k;
         }
         else if (matrix[i][j].v > (matrix[i][k - 1].v + matrix[k][j].v))
         {
            matrix[i][j].v = matrix[i][k - 1].v + matrix[k][j].v;
            matrix[i][j].k = k;
         }

      }
   }
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
            calc_cell_nussinov(i, j, sequence, matrix, scores);
         }
      }
   }

   for (i = 0; i < seqlen; i++)
   {
      for (j = 0; j < seqlen; j++)
      {
         mprintf ("%.2f ", matrix[i][j].v);
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
      mprintf ("%c ", transform_number_2_base (sequence[i]));
   }
   mprintf ("\n");
   for (i = 0; i < seqlen; i++)
   {
      for (j = 0; j < seqlen; j++)
      {
         if (matrix[i][j].k != ULONG_UNDEF)
         {
            mprintf ("%lu ", matrix[i][j].k);
         }
         else
         {
            mprintf ("x ");
         }
      }
      mprintf (" %c\n", transform_number_2_base (sequence[i]));
      }

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
         traceback_matrix (i + 1, j, structure, matrix);
      }
      else if (matrix[i][j].from & L)
      {
         traceback_matrix (i, j - 1, structure, matrix);
      }
      else if (matrix[i][j].from & D)
      {
         /* mfprintf (stderr, "PAIR: %lu,%lu\n", i, j); */
         str_at (structure, i, '(');
         str_at (structure, j, ')');
         traceback_matrix (i + 1, j - 1, structure, matrix);
      }
      else if (matrix[i][j].k != ULONG_UNDEF)
      {
         traceback_matrix (i, matrix[i][j].k - 1,
                           structure,
                           matrix); 
         traceback_matrix (matrix[i][j].k, j,
                           structure,
                           matrix);         
      }
   }
}

static void
traceback_nussinov (const unsigned long i, 
                    const unsigned long j,
                    Str* structure,
                    const char* sequence,
                    struct_cell** matrix,
                    float** scores)
{
   unsigned long k;

   if (i < j)
   {
      if (fabsf(matrix[i][j].v - matrix[i + 1][j].v) < 0.000001)
      {
         traceback_nussinov (i + 1, j, structure, sequence, matrix, scores);
      }
      else
      {
         if (fabsf(matrix[i][j].v - matrix[i][j - 1].v) < 0.000001)
         {
            traceback_nussinov (i, j - 1, structure, sequence, matrix, scores);
         }
         else
         {
/*         mprintf ("B: %f == %f + %f\n", matrix[i][j].v, matrix[i+1][j-1].v, */
/*                      scores[(int)sequence[i]][(int)sequence[j]]); */
            
            if (fabsf(matrix[i][j].v - (matrix[i + 1][j - 1].v 
                                  + scores[(int)sequence[i]][(int)sequence[j]]))
                < 0.000001)
            {
               /* mfprintf (stderr, "PAIR: %lu,%lu\n", i, j); */
               str_at (structure, i, '(');
               str_at (structure, j, ')');
               traceback_nussinov (i + 1, j - 1,
                                   structure,
                                   sequence,
                                   matrix,
                                   scores);
            }
            else
            {
               for (k = i + 1; k < j; k++)
               {
                  if (fabsf(  matrix[i][j].v 
                           - (matrix[i][k - 1].v + matrix[k][j].v)) < 0.000001)
                  {
                     traceback_nussinov (i, k - 1,
                                         structure,
                                         sequence,
                                         matrix,
                                         scores); 
                     traceback_nussinov (k, j,
                                         structure,
                                         sequence,
                                         matrix,
                                         scores);
                     return;
                  }
               }
            }
         }
      }
   }
}

static Str*
pred_2D_structure_nussinov (const unsigned long l,
                            const char* sequence,
                            const unsigned long seqlen)
{
   struct_cell** matrix;
   float** scores;
   unsigned long i, j;
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
      str_delete (structure);
      return NULL;
   }

   for (i = 0; i < seqlen; i++)
   {
      for (j = 0; j < seqlen; j++)
      {
         matrix[i][j].v = 0.0f;
         matrix[i][j].from = 0;
         matrix[i][j].k = ULONG_UNDEF;
      }
   }

   /* fill matrix */
   scores = create_scoring_matrix ();
   if (scores == NULL)
   {
      XFREE_2D ((void**)matrix);
      str_delete (structure);
      return NULL;
   }

   if (calc_matrix_nussinov (sequence, matrix, seqlen, l, scores))
   {
      XFREE_2D ((void**)matrix);
      XFREE_2D ((void**)scores);
      str_delete (structure);
      return NULL;
   }

   /* traceback */
   if (seqlen > 0)
   {
      mprintf ("> MFE: %.3f\n", matrix[0][seqlen - 1].v);
    traceback_nussinov (0, seqlen - 1, structure, sequence, matrix, scores);
    /*traceback_matrix (0, seqlen - 1, structure, matrix);*/
   }

   XFREE_2D ((void**)matrix);
   XFREE_2D ((void**)scores);

   return structure;
}

int
fold_main(const char *cmdline)
{
   struct fold_args_info fold_args;
   Str* structure;
   unsigned long seqlen, i;
   int retval = 0;

   /* command line parsing */
   fold_cmdline_parser_init (&fold_args);

   retval = fold_cmdline_parser_string (cmdline, &fold_args, get_progname());

   if (retval == 0)
   {
      retval = fold_cmdline_parser_required ();
   }

   /* postprocess arguments */
   seqlen = strlen (fold_args.inputs[1]);
   if (retval == 0)
   { 
      retval = fold_cmdline_parser_postprocess (&fold_args);
   }

   /* do work */
   if (retval == 0)
   { 
      structure = pred_2D_structure_nussinov (fold_args.loop_length_arg,
                                              fold_args.inputs[1],
                                              seqlen);
      if (structure == NULL)
      {
         retval = 1;
      }
      else
      {
         for (i = 0; i < seqlen; i++)
         {
            mprintf ("%c", transform_number_2_base (fold_args.inputs[1][i]));
         }
         mprintf ("\n");
         mprintf ("%s\n", str_get (structure));
      }
   }

   /* finalise */
   if (retval == 0)
   { 
      str_delete (structure);
   }
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
