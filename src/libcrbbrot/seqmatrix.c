/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbrot/seqmatrix.c
 *
 *  @brief Sequence matrix for SCMF
 *
 *  Module: seqmatrix
 *
 *  Library: libcrbbrot
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-03-05
 *
 *
 *  Revision History:
 *         - 2008Mar05 bienert: created
 *
 */


#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <libcrbbasic/crbbasic.h>
#include <libcrbrna/crbrna.h>
#include "preset.h"
#include "seqmatrix.h"

/* use until something better is found */
#define SM_ROWS 4               /* HP: 2 letters, RNA: 4 */
#define R_GAS 1 /*8.314472*/

struct SeqMatrix {
      char* fixed_sites;
      float** matrix[2];
      short int curr_matrix;
      unsigned long* pairlist;
      char* sequence;
      size_t rows;
      size_t cols;   
};


/**********************   Constructors and destructors   **********************/

/** @brief Create a new sequence matrix.
 *
 * The constructor for @c SeqMatrix objects. If compiled with memory checking
 * enabled, @c file and @c line should point to the position where the function
 * was called. Both parameters are automatically set by using the macro
 * @c SEQMATRIX_NEW.\n
 * Returns @c NULL on error.
 *
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
SeqMatrix*
seqmatrix_new (const char* file, const int line)
{
   /* allocate 1 object */
   SeqMatrix* sm = XOBJ_MALLOC(sizeof (SeqMatrix), file, line);

   if (sm != NULL)
   {
      sm->fixed_sites = NULL;
      sm->pairlist    = NULL;
      sm->sequence    = NULL;
      sm->matrix[0]   = NULL;
      sm->matrix[1]   = NULL;
      sm->curr_matrix = 0;
   }

   return sm;
}

/** @brief Delete a sequence matrix.
 *
 * The destructor for @c SeqMatrix objects.
 *
 * @param[in] sm object to be freed.
 */
void
seqmatrix_delete (SeqMatrix* sm)
{
   if (sm != NULL)
   {
      XFREE    (sm->fixed_sites);
      XFREE_2D ((void**)sm->matrix[0]);
      XFREE_2D ((void**)sm->matrix[1]);
      XFREE    (sm->pairlist);
      XFREE    (sm->sequence);
      XFREE    (sm);
   }
}


/*********************************   Access   *********************************/
/** @brief Check if a row is set as 'fixed'.
 *
 * Check if a certain row of a sequence matrix is preempted with a fixed site.\n
 * Returns 'true', if the site is set, 'false' otherwise.
 *
 * @param[in] row The row to check.
 * @param[in] sm The sequence matrix
 */
__inline__ bool
seqmatrix_is_col_fixed (unsigned long col, SeqMatrix* sm)
{
   assert (sm);
   assert (sm->matrix[sm->curr_matrix]);
   assert (sm->fixed_sites);
   /* assert (); */

   if (sm->fixed_sites[(col / CHAR_BIT)] & (1 << (col % CHAR_BIT)))
   {
      return true;
   }

   return false;
}


/********************************   Altering   ********************************/

/** @brief Initialise a sequence matrix.
 *
 * On initialisation, all probabilities of the matrix are set. Additionally, a
 * set of sites can be fixed to exclude them from the simulation.\n
 * Returns 0 on success, @c ERR_SM_ALLOC indicating memory problems.
 *
 * @param[in] size Size of the structure.
 * @param[in] presettings List of fixed sites.
 * @param[in] sm Sequence matrix to initialise.
 */
int
seqmatrix_init (const unsigned long* pairs,
                const unsigned long size,
                const PresetArray* presettings,
                SeqMatrix* sm)
{
   unsigned long i, j;
   unsigned long pslen;
   unsigned long row, col;
   /* float whobble = 0.001f; */

   assert (sm);
   assert (sm->fixed_sites == NULL);
   assert (sm->pairlist    == NULL);
   assert (sm->matrix[0]   == NULL);
   assert (sm->matrix[1]   == NULL);

   /* allocate matrix */
   sm->rows = SM_ROWS;
   sm->cols = size;

   sm->matrix[0] = (float**) XMALLOC_2D (sm->rows,sm->cols,
                                         sizeof (float));
   if (sm->matrix[0] == NULL)
   {
      return ERR_SM_ALLOC;
   }
   sm->curr_matrix = 0;

   sm->matrix[1] = (float**) XMALLOC_2D (sm->rows,sm->cols,
                                         sizeof (float));
   if (sm->matrix[1] == NULL)
   {
      return ERR_SM_ALLOC;
   }

   /* init sm */
   if (size > 0)
   {
      sm->matrix[0][0][0] = 1.0f / SM_ROWS; /* even distributed init */
      sm->matrix[1][0][0] = 0.0f;
      for (i = 1; i < size; i++)
      {
         sm->matrix[0][0][i] = sm->matrix[0][0][0]; /* even distri */
         sm->matrix[1][0][i] = sm->matrix[1][0][0];
      }

      for (i = 1; i < sm->rows; i++)
      {
         memcpy (sm->matrix[0][i],
                 sm->matrix[0][i - 1],
                 sizeof (float) * size);/* even distri */
         memcpy (sm->matrix[1][i],
                 sm->matrix[1][i - 1],
                 sizeof (float) * size);
      }

      /*srand(997654329);
      for (j = 0; j < size; j++)
      {
         for (i = 0; i < SM_ROWS; i++)
         {
            sm->matrix[0][i][j] = rand();
            sm->matrix[1][0][j] += sm->matrix[0][i][j];
         }
         for (i = 0; i < SM_ROWS; i++)
         {
            sm->matrix[0][i][j] = sm->matrix[0][i][j] / sm->matrix[1][0][j];
         }
         sm->matrix[1][0][j] = 0.0f;
         }*/

   }

   /* copy pairlist */
   sm->pairlist = XCALLOC (size, sizeof (*(sm->pairlist)));

   if (pairs != NULL)
   {
      for (i = 0; i < size; i++)
      {
         sm->pairlist[i] = pairs[i]; /* xxx use memcpy instead? */
      }
   }

   /* run over list of presettings and set sites */
   /* allocate memory & init fixed sites */
   sm->fixed_sites = XCALLOC ((size / (sizeof(*(sm->fixed_sites))*CHAR_BIT))+1,
                              sizeof(*(sm->fixed_sites)));

   if (sm->fixed_sites == NULL)
   {
      return ERR_SM_ALLOC;
   }

   if (presettings != NULL)
   {
      /* use, if rows are adjacent in memory
        init_row[0] = 0.0f;
        for (i = 1; i < SM_ROWS; i++)
        {
        init_row[i] = init_row[0];
        }
      */

      pslen = presetarray_get_length (presettings);
      for (i = 0; i < pslen; i++)
      {
         row = presetarray_get_ith_base (i, presettings);
         col = presetarray_get_ith_pos (i, presettings);

         /* set site as fixed */
         sm->fixed_sites[(col / CHAR_BIT)] |= 
            1 << (col % CHAR_BIT);

         /* set probability in matrix to 1 */
         /* use, if rows are the connected memory areas
           memcpy (sm->matrix[row], init_row, sizeof (*(sm->matrix)) * SM_ROWS);
         */
         if (size > 0)
         {
            for (j = 0; j < SM_ROWS; j++)
            {
               sm->matrix[0][j][col] = 0.0f;
            }
            sm->matrix[0][row][col] = 1.0f;
            sm->matrix[1][row][col] = 1.0f;
         }
      }
   }

   return 0;
}


/** @brief Update a row of a sequence matrix
 *
 * Update a row of a sequence matrix in a SCMF simulation.
 * Returns...
 *
 * @params[in] sm Sequence matrix.
 */
static __inline__ int
sequence_matrix_calc_eeff_row_scmf (unsigned long col,
                                    SeqMatrix* sm,
                                    float t,
                                    float** scores)
{
   unsigned long j;
   unsigned long l;
   unsigned long k;
   short int new_matrix = 0;
   float tmp_neg;               /* negative energy contribution */
   float tmp_het;               /* heterogenity term */
   float het_count;
   float het_rate;

   het_rate = ((-1) * ((logf (1 / 0.000001f)) / (sm->cols)));

   if (sm->curr_matrix == 0)
   {
      new_matrix = 1;
   }

   /*fprintf (stderr, "C: %lu Pairs with %lu\n", col, sm->pairlist[col]); */

   for (j = 0; j < sm->rows; j++)
   {
      /* update cell */
      sm->matrix[new_matrix][j][col] = 0.0f;
     /*  mfprintf (stderr, "   R:%lu ", j); */

      /* calculate contribution of wanted interaction (if any) */
      if (sm->pairlist[col])
      {
         for (l = 0; l < sm->rows; l++)
         {
            if (col < sm->pairlist[col])
            {
               sm->matrix[new_matrix][j][col] +=
                  (sm->matrix[sm->curr_matrix][l][sm->pairlist[col] - 1] 
                   * scores[j][l]);
               /*mfprintf (stderr, "(%2.3f*%2.3f:%2.6f) ", scores[j][l],
                         sm->matrix[sm->curr_matrix][l][sm->pairlist[col] - 1],
                         sm->matrix[new_matrix][j][col]);*/
            }
            else
            {
               sm->matrix[new_matrix][j][col] +=
                  (sm->matrix[sm->curr_matrix][l][sm->pairlist[col] - 1] 
                   * scores[l][j]);
               /*mfprintf(stderr, "(%2.3f*%2.3f:%2.6f) ", scores[l][j],
                         sm->matrix[sm->curr_matrix][l][sm->pairlist[col] - 1],
                         sm->matrix[new_matrix][j][col]);*/          
            }
         }
      }

      /* calculate contribution of unwanted pairs */
      tmp_neg = 0.0f;
      tmp_het = 0.0f;
      het_count = 0.0f;
      for (k = 0; k < sm->cols; k++)
      {
         if ((k != col) && ((k + 1) != sm->pairlist[col])) /*without k!=col ?*/
         {
            for (l = 0; l < sm->rows; l++)
            {
               if (col < k)
               {
                  tmp_neg +=   sm->matrix[sm->curr_matrix][l][k] 
                     * scores[j][l];
               }
               else
               {
                  tmp_neg +=   sm->matrix[sm->curr_matrix][l][k] 
                     * scores[l][j];
               }
            }
            
            /* heterogenity term */
            /* static version: use window*/
            /*   if (  ((k < col)&&(col - k < 4)) *//*3 5*//*
                ||((k > col)&&(k - col < 4)))
            {
              tmp_het += (1.079 * sm->matrix[sm->curr_matrix][j][k]);
              }*/
            /* variable version 1: linear decrease
            tmp_het += (((float)(sm->cols - k) / sm->cols)
                        * sm->matrix[sm->curr_matrix][j][k]);
                        het_count += (float)(sm->cols - k) / sm->cols;*/

            /* variable version 2: exp decrease */
            /* not k u idiot: distance to current position */
            if (col > k)
            {
               tmp_het += (sm->matrix[sm->curr_matrix][j][k]
                             * expf(het_rate * (col - (k+1))));
               het_count += expf ( (het_rate* (col - (k+1))));
               /* mfprintf (stderr, ":%.3f", tmp_het); */
            }
            else
            {
               tmp_het += (sm->matrix[sm->curr_matrix][j][k]
                            * expf(het_rate * (k - (col+1))));
               het_count += (expf (het_rate * (k - (col+1))));

               /* mfprintf (stderr, ":%.3f", tmp_het); */
            }
         }
      }
      tmp_neg = (tmp_neg / sm->cols) * (-1.25f); /* 1.25 */
      /* mfprintf (stderr, "tmp_neg: %3f ", tmp_neg); */
      tmp_het = (tmp_het / het_count) * (3.0f); /* 0.1747 */
      /* mfprintf (stderr, "tmp_het: %3f ", tmp_het); */
      sm->matrix[new_matrix][j][col] += tmp_neg;
      sm->matrix[new_matrix][j][col] += tmp_het;
      sm->matrix[new_matrix][j][col] =
         expf ((-1.0f) * (sm->matrix[new_matrix][j][col]/(R_GAS * t)));
      /* mfprintf (stderr, "\n"); */
   }
   /* mfprintf (stderr, "\n"); */

/*T: 0.100002 S: -0.318406*/
/*3220130003312003000003022030232000000032312003002133200030003220312013320003*/

/*3.0 | -0.056336*/
/*3302013000322000000003321021302000000031203000002233200000003223320131220000*/

/*6.0 | -0.116815*/
/*3020213003122000000003302032132000000032032000002133200000003220320313120000*/


   return 0;
}


/** @brief Update the columns of a sequence matrix
 *
 * Update cols of a sequence matrix in a SCMF simulation.
 * Returns...
 *
 * @params[in] sm Sequence matrix.
 */
static __inline__ int
sequence_matrix_calc_eeff_col_scmf (SeqMatrix* sm, float t, float** scores)
{
   unsigned long i;
   int error = 0;

   /* for each col */
   i = 0;
   while ((!error) && (i < sm->cols))
   {
      if (!seqmatrix_is_col_fixed (i, sm))
      {
         error = sequence_matrix_calc_eeff_row_scmf (i, sm, t, scores);        
      }
      /*else
      {
         mfprintf (stderr, "C: %lu Pairs with %lu\n", i, sm->pairlist[i]);
         }*/
      i++;
   }

   if (sm->curr_matrix == 0)
   {
      sm->curr_matrix = 1;
   }
   else
   {
      sm->curr_matrix = 0;
   }

   return error;
}

/** @brief Transform a sequence matrix into an unambigouos sequence.
 *
 * Creates a sequence out of a sequence matrix in an iterative simulation
 * process. As a start, all unambigouos columns are fixed. Then the matrix is
 * used in a new SCMF process with the "most" unambigouos of the undecided
 * columns fixed by majority voting. This procedure is repeated until... In
 * the end the matrix is compressed into a sequence.
 *
 * @params[in] sm The sequence matrix.
 */
int
seqmatrix_collate_is (float fthresh,
                      unsigned long steps,
                      float temp,
                      float** scores,
                      SeqMatrix* sm,
                      Alphabet* sigma)
{
   unsigned long i, j, k, largest_amb_col = 0;
   float largest_amb;
   int retval = 0;

   assert (sm);
   assert (sm->sequence == NULL);
   assert (fthresh <= 1.0f);

   /* Approach: find unambigouos sites and fixate 'em */
   /*           find the largest of the ambigouos sites */
   /*           until all sites are fixed */
   /* mfprintf (stderr, "IN\n"); */

   while ((largest_amb_col != sm->cols + 1) && (! retval))
   {
      largest_amb_col = sm->cols + 1;
      largest_amb = 0.0f;

      /* for all columns */
      for (j = 0; j < sm->cols; j++)
      {
         if (!seqmatrix_is_col_fixed (j, sm))
         {
            /* for all rows */
            for (i = 0; i < sm->rows; i++)
            {
               if (sm->matrix[sm->curr_matrix][i][j] >= fthresh)
               {
                  /* mfprintf (stderr, "GT: %lu %lu %.3f\n", i, j, */
/*                   sm->matrix[sm->curr_matrix][i][j]); */
                  /* unambigouos site found, fixate it */
                  sm->fixed_sites[(j / CHAR_BIT)] |= 
                     1 << (j % CHAR_BIT);
                  /* set remaining pos to 0! */
                  for (k = 0; k < sm->rows; k++)
                  {
                     sm->matrix[0][k][j] = 0.0f;
                     sm->matrix[1][k][j] = 0.0f;
                  }
                  sm->matrix[0][i][j] = 1.0f;
                  sm->matrix[1][i][j] = 1.0f;
                  i = sm->rows + 1;
               }
            }
         }
      }
      
      /* seqmatrix_print_2_stderr (3, sm); */
      
      /* find largest ambigouos site */
      for (j = 0; j < sm->cols; j++)
      {
         /* check if site is fixed */
         if (!seqmatrix_is_col_fixed (j, sm))
         {
            for (i = 0; i < sm->rows; i++)
            {
               if (sm->matrix[sm->curr_matrix][i][j] > largest_amb)
               {
                  largest_amb = sm->matrix[sm->curr_matrix][i][j];
                  largest_amb_col = j;
               }
            }
         }
         /*else
         {
            mfprintf (stderr, "Fixed site: %lu\n", j);
         }*/
      }
      
      /* set site to 1/0 and fixate it */
      if (largest_amb_col < sm->cols + 1)
      {
         sm->fixed_sites[(largest_amb_col / CHAR_BIT)] |=
            1 << (largest_amb_col % CHAR_BIT);
         for (k = 0; k < sm->rows; k++)
         {
            if (sm->matrix[sm->curr_matrix][k][largest_amb_col] == largest_amb)
            {
               sm->matrix[0][k][largest_amb_col] = 1.0f;
               sm->matrix[1][k][largest_amb_col] = 1.0f;
            }
            else
            {
               sm->matrix[0][k][largest_amb_col] = 0.0f;
               sm->matrix[1][k][largest_amb_col] = 0.0f;
            }
         }
         /* simulate */
         retval = sequence_matrix_simulate_scmf (steps,
                                                 temp,
                                                 sm,
                                                 scores);
      }
   }

   if (! retval)
   {
      retval = seqmatrix_collate_mv (sm, sigma);
   }

   return retval;
}

/** @brief Transform a sequence matrix into an unambigouos sequence.
 *
 * Collates all rows of a column to a single representative. Thereby the row
 * with the highest share is chosen. The sequence is stored in the SeqMatrix
 * object an can be accessed via...\n
 * Returns @c ERR_SM_ALLOC on problems allocating memory for the sequence. 0
 * otherwise.
 *
 * @params[in] sm The sequence matrix.
 */
int
seqmatrix_collate_mv (SeqMatrix* sm, Alphabet* sigma)
{
   unsigned long i,j;
   float curr_max_prob;

   assert (sm);
   assert (sm->sequence == NULL);

   /* try to init sequence */
   sm->sequence = XMALLOC (sm->cols * sizeof (*(sm->sequence)));
   if (sm->sequence == NULL)
   {
      return ERR_SM_ALLOC;
   }


   /* for all columns */
   for (j = 0; j < sm->cols; j++)
   {
      curr_max_prob = 0.0f;
      for (i = 0; i < sm->rows; i++)
      {
         /* find highest number */
         if (curr_max_prob < sm->matrix[sm->curr_matrix][i][j])
         {
            /* write position to seq */
            curr_max_prob = sm->matrix[sm->curr_matrix][i][j];
            sm->sequence[j] = alphabet_no_2_base(i, sigma);
/*3203211002132030030103203032032001030032132003002303201303003212300321320003*/
         }
      }
   }

   return 0;
}

/** @brief Perform a SCMF simulation on a sequence matrix.
 *
 * Calculate the mean force field for a sequence matrix and update cells. This
 * is done either until the system converges or for a specified number of
 * steps. For results of the simulation, the matrix itself has to be
 * interpreted.\n
 * Returns ...
 *
 * @params[in] steps No. of max. simulation steps.
 * @params[in] sm The sequence matrix.
 */
int
sequence_matrix_simulate_scmf (const unsigned long steps,
                               const float t_init,
                               SeqMatrix* sm,
                               float** scores)
{
   unsigned long t;
   int error= 0;
   unsigned long j, i;
   float row_sum;
   float T = t_init;
   float c_rate = 1.0f;
   float c_port = 0.0f;
   short int m = 0;
   float s;

   assert (sm);
   assert (scores && scores[0]);

   if (sm->curr_matrix == 0)
   {
     m = 1;
   }

   if (steps > 0)
   {
      c_rate = expf ((-1) * ((logf (t_init / 0.1f)) / (steps - 1)));
      /* c_port = (-1) * (t_init - 0.1f) / (steps - 1); */
   }

   /* perform for a certain number of steps */
   t = 0;
   while ((!error) && (t < steps))
   {
      s = 0.0f;

      /* mfprintf (stderr, "Step %lu Temp %3f cool %3f:\n", t, T, c_port); */
      /* calculate Eeff */
      error = sequence_matrix_calc_eeff_col_scmf (sm, T, scores);
      /* seqmatrix_print_2_stderr (6, sm); */
      /* update matrix */
      /* for all columns */
      for (j = 0; j < sm->cols; j++)
      {
         /* which are not fixed */
         if (!seqmatrix_is_col_fixed (j, sm))
         {
            /* calc sum of col */
            row_sum = 0.0f;
            for (i = 0; i < sm->rows; i++)
            {
               row_sum += sm->matrix[sm->curr_matrix][i][j];
            }

            /* for each row */
            for (i = 0; i < sm->rows; i++)
            {
               /* avoid oscilation by Pnew = uPcomp + (1 - u)Pold) */
               sm->matrix[sm->curr_matrix][i][j] = 
                  sm->matrix[sm->curr_matrix][i][j] / row_sum;
               /* avoid oscilation by Pnew = uPcomp + (1 - u)Pold) */
               sm->matrix[sm->curr_matrix][i][j] = 
                  (0.2 * sm->matrix[sm->curr_matrix][i][j])
                  + (0.8 * sm->matrix[m][i][j]);

               /* calculate "entropy", ignore fixed sites since ln(1) = 0 */
               s += (sm->matrix[sm->curr_matrix][i][j]
                     * logf (sm->matrix[sm->curr_matrix][i][j]));
            }
         }
         /*else
         {
            mfprintf (stderr, "SFixed site: %lu - ", j);
            for (i = 0; i < sm->rows; i++)
            {
               mfprintf (stderr, "%.3f ", sm->matrix[sm->curr_matrix][i][j]);
            }
            mfprintf (stderr, "\n");
            }*/
      }

      s = s/sm->cols;
      /* mfprintf (stdout, "T: %3f S: %3f\n", T, s); */

      t++;
      /* mfprintf (stderr, "\n"); */

     /* seqmatrix_print_2_stderr (3, sm); */
      /* if (s > -0.01f) return error; */
     T = (T * c_rate) + c_port;
   }

   /* transform matrix to flat sequence */
   /* error = seqmatrix_collate_is (0.99, scores, sm); */
   /*error = seqmatrix_collate_mv (sm);*/

   return error;
}


/*********************************   Output   *********************************/

/** @brief Print the sequence of a sequence matrix to a stream.
 *
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] sm The sequence matrix.
 */
void
seqmatrix_fprintf_sequence (FILE* stream, SeqMatrix* sm)
{
   unsigned long i;

   assert (sm);
   assert (sm->sequence);

   /* mfprintf (stream, "%lu\n", (unsigned long) sm->cols); */
   for (i = 0; i < sm->cols; i++)
   {
      mfprintf (stream, "%c", sm->sequence[i]);
   }
}

/** @brief Print the sequence of a sequence matrix to stdout.
 *
 * @params[in] sm The sequence matrix.
 */
void
seqmatrix_printf_sequence (SeqMatrix* sm)
{
   seqmatrix_fprintf_sequence (stdout, sm);
}

/** @brief Print a sequence matrix to a stream.
 *
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] p  Precision of the numbers.
 * @params[in] sm The sequence matrix.
 */
void
seqmatrix_fprintf (FILE* stream, const int p, const SeqMatrix* sm)
{
   unsigned long line_width = 0;
   unsigned long i, j;
   int rprec = 0;               /* precision for row number */
   int cprec = 0;               /* precision for cells */
   float tmp;
   char* string;
   char* string_start;

   assert (sm);
   assert (p <= FLT_DIG); /* check that p stays in the precision of float */
   
   /* to circumvent problems with multiple threads, interfering prints...
      we first create the whole matrix in the memory */
   
   /* find the widest line */
   for (i = 0; i < sm->rows; i++)
   {
      for (j = 0; j < sm->cols; j++)
      {
         if (sm->matrix[sm->curr_matrix][i][j] < 0.0f)
         {
            tmp = (-1) * sm->matrix[sm->curr_matrix][i][j];
            rprec = 1;
         }
         else
         {
            tmp = sm->matrix[sm->curr_matrix][i][j];
            rprec = 0;
         }

         if (tmp > 1.0f)
         {
            rprec += floor (log10 (tmp) + 1.0f);
         }
         else
         {
            rprec += 1;
         }

         if (rprec > cprec) 
         {
            cprec = rprec;
         }      
      }
   }
   cprec += p + 1;

   if (sm->rows > 1.0f)
   {
      rprec = floor (log10 (sm->rows) + 1.0f);
   }
   else
   {
      rprec = 1;
   }

   /* add up components of a line */
   line_width += rprec;
   line_width += 4;              /* + ':  |' row start */
   line_width += ((3 + cprec) * sm->cols); /* + '.' + '  |' */
   line_width += 1;              /* + '\n' */

   /* alloc memory for the string */
   string = XMALLOC (sizeof (char) * (((line_width) * sm->rows) + 1));
   string_start = string;

   /* write matrix */
   for (i = 0; i < sm->rows; i++)
   {
      msprintf (string, "%*lu:  | ", rprec, i);
      string += rprec + 4;

      for (j = 0; j < sm->cols; j++)
      {
         msprintf (string, " %*.*f |",
                   cprec,
                   p,
                   sm->matrix[sm->curr_matrix][i][j]);
         string += 3 + cprec;
      }

      string[0] = '\n';
      string++;
   }   
   string[0] = '\0';

   string = string_start;
   mfprintf (stream, "%s\n", string);

   XFREE (string);
}


/** @brief Print a sequence matrix to stdout.
 *
 * @params[in] p  Precision of the numbers.
 * @params[in] sm The sequence matrix.
 */
void
seqmatrix_print_2_stdout (const int p, const SeqMatrix* sm)
{
   seqmatrix_fprintf (stdout, p, sm);
}


/** @brief Print a sequence matrix to stderr.
 *
 * @params[in] p  Precision of the numbers.
 * @params[in] sm The sequence matrix.
 */
void
seqmatrix_print_2_stderr (const int p, const SeqMatrix* sm)
{
   seqmatrix_fprintf (stderr, p, sm);
}

/*3110002002222000000003333022222000000033333000002222200000003333331110020000*/
/*1000222002222000000003333022222000000033333000002222200000003333333311100000*/
/*1002222002222000000003333022222000000033333000002222200000003333333331100000*/
/*1133333330000333333331111300000333333311111333330000033333331111122222003333*/
/*1113333330000333333331111300000333333311111333330000033333331111133330003333*/
/*1113333330000333333331111300000333333311111333330000033333331111122220003333*/
/*1113333330000333333331111300000333333311111333330000033333331111122220003333*/
/*1113333330000333333331111300000333333311111333330000033333331111133330003333*/
/*1113333330000333333331111300000333333311111333330000033333331111122220003333*/
/*1113333330000333333331111300000333333311111333330000033333331111122220003333*/
/*1113333330000333333331111300000333333311111333330000033333331111122220003333*/
/*1113333330000333333331111300000333333311111333330000033333331111122220003333*/
/*0022222002222000000003333022222000000033333000002222200000003333333333110000*/
/*2222222002222000000003333022222000000033333000002222200000003333333333330000*/
/*2222222002222000000003333022222000000033333000002222200000003333333333330000*/
/*2222222002222000000003333022222000000033333000002222200000003333333333330000*/
/*2222222002222000000003333022222000000033333000002222200000003333333333330000*/
/*3311111002222000000003333022222000000033333000002222200000003333300000220000*/
/*CCUUUUUAAGGGGAAAAAAAACCCCAGGGGGAAAAAAACCCCCAAAAAGGGGGAAAAAAACCCCCAAAAAGGAAAA*/
/*(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....*/
/*3203211002132030030103203032032001030032132003002303201303003212300321320003*/
/*1133220001332300001133220001332100013332201033012213300031132203313322000133*/
/*3102233002233000133002233002233002130022331000331022300331002331022331020330*/
/*0320312002312003003003023023023003003023123003001320330033002132030213210300*/



/*GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA*/
/*3102233002233001233002233002231002331002331002331022310302102331022331020310*/

/*100k steps: 3102233002233000133002233002233002130022331000331022300331002331022331020330*/
