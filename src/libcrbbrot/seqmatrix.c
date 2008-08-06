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
 * ToDo: seqmatrix_init via macro because it allocates memory!!!
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
#include "seqmatrix.h"

/* use until something better is found */
#define SM_ROWS 4               /* HP: 2 letters, RNA: 4 */
#define R_GAS 1 
#define R_GAS_NN 8.314472

/* matrix orders */
enum seqmatrix_matrix_no{
   F_Mtrx = 0,      /* first matrix */
   S_Mtrx,          /* second matrix */
   No_Of_Mtrx       /* overall number of matrices */
};

struct SeqMatrix {
      char* fixed_sites;
      float** matrix[No_Of_Mtrx];
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
   unsigned long i;

   /* allocate 1 object */
   SeqMatrix* sm = XOBJ_MALLOC(sizeof (SeqMatrix), file, line);

   if (sm != NULL)
   {
      sm->fixed_sites  = NULL;
      sm->pairlist     = NULL;
      sm->sequence     = NULL;

      for (i = 0; i < No_Of_Mtrx; i++)
      {
         sm->matrix[i] = NULL;
      }
      sm->curr_matrix = F_Mtrx;
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
   unsigned long i;

   if (sm != NULL)
   {
      XFREE    (sm->fixed_sites);
      XFREE    (sm->pairlist);
      XFREE    (sm->sequence);

      for (i = 0; i < No_Of_Mtrx; i++)
      {
         XFREE_2D ((void**)sm->matrix[i]);
      }

      XFREE    (sm);
   }
}


/*********************************   Access   *********************************/

/** @brief Check if a column is set as 'fixed'.
 *
 * Check if a certain column of a sequence matrix is preempted with a fixed
 * site.\n
 * Returns 'true', if the site is set, 'false' otherwise.
 *
 * @param[in] col The column to check.
 * @param[in] sm The sequence matrix
 */
__inline__ bool
seqmatrix_is_col_fixed (const unsigned long col, const SeqMatrix* sm)
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

/** @brief get the width of a sequence matrix.
 *
 * Returns the number of columns of a sequence matrix.
 *
 * @param[in] sm sequence matrix.
 */
unsigned long
seqmatrix_get_width (const SeqMatrix* sm)
{
   assert (sm);

   return sm->cols;
}

/********************************   Altering   ********************************/

/** @brief Fix a certain column in a sequence matrix.
 *
 * Sets a certain cell of a sequenc ematrix to 1 and marks it to be preserved.
 * 
 * @param[in] row row to be fixed.
 * @param[in] col column to be fixed.
 * @param[in] sm sequence matrix.
 */
void
seqmatrix_fix_col (const unsigned long row, const unsigned long col,
                   SeqMatrix* sm)
{
   unsigned long i;

   assert (sm);
   assert (sm->fixed_sites);
   assert (col < sm->cols);
   assert (row < sm->rows);

   /* fix site */
   sm->fixed_sites[(col / CHAR_BIT)] |= 1 << (col % CHAR_BIT);

   /* set everything to 0 */
   for (i = 0; i < sm->rows; i++)
   {
      sm->matrix[F_Mtrx][i][col] = 0.0f;
      sm->matrix[S_Mtrx][i][col] = 0.0f;
   }

   /* set demand to 1 */
   sm->matrix[F_Mtrx][row][col] = 1.0f;
   sm->matrix[S_Mtrx][row][col] = 1.0f;
}


/** @brief Initialise a sequence matrix.
 *
 * On initialisation, all probabilities of the matrix are set. Additionally, a
 * set of sites can be fixed to exclude them from the simulation.\n
 * Returns 0 on success, @c ERR_SM_ALLOC indicating memory problems.
 *
 * @param[in] pairs base pairs.
 * @param[in] size Size of the structure.
 * @param[in] sm Sequence matrix to initialise.
 */
int
seqmatrix_init (const unsigned long* pairs,
                const unsigned long size,
                SeqMatrix* sm)
{
   unsigned long i/* , j */;
   /* float whobble = 0.001f; */

   assert (sm);
   assert (sm->fixed_sites == NULL);
   assert (sm->pairlist    == NULL);
   assert (sm->matrix[F_Mtrx]   == NULL);
   assert (sm->matrix[S_Mtrx]   == NULL);

   /* allocate matrix */
   sm->rows = SM_ROWS;
   sm->cols = size;

   sm->matrix[F_Mtrx] = (float**) XMALLOC_2D (sm->rows,sm->cols,
                                              sizeof (float));
   if (sm->matrix[F_Mtrx] == NULL)
   {
      return ERR_SM_ALLOC;
   }
   sm->curr_matrix = F_Mtrx;

   sm->matrix[S_Mtrx] = (float**) XMALLOC_2D (sm->rows,sm->cols,
                                         sizeof (float));
   if (sm->matrix[S_Mtrx] == NULL)
   {
      return ERR_SM_ALLOC;
   }

   /* init sm */
   if (size > 0)
   {
      sm->matrix[F_Mtrx][0][0] = 1.0f / SM_ROWS; /* even distributed init */
      sm->matrix[S_Mtrx][0][0] = 0.0f;
      for (i = 1; i < size; i++)
      {
         sm->matrix[F_Mtrx][0][i] = sm->matrix[F_Mtrx][0][0]; /* even distri */
         sm->matrix[S_Mtrx][0][i] = sm->matrix[S_Mtrx][0][0];
      }

      for (i = 1; i < sm->rows; i++)
      {
         memcpy (sm->matrix[F_Mtrx][i],
                 sm->matrix[F_Mtrx][i - 1],
                 sizeof (float) * size);/* even distri */
         memcpy (sm->matrix[S_Mtrx][i],
                 sm->matrix[S_Mtrx][i - 1],
                 sizeof (float) * size);
      }

      /*srand(997654329);
      for (j = 0; j < size; j++)
      {
         for (i = 0; i < SM_ROWS; i++)
         {
            sm->matrix[F_Mtrx][i][j] = rand();
            sm->matrix[S_Mtrx][0][j] += sm->matrix[F_Mtrx][i][j];
         }
         for (i = 0; i < SM_ROWS; i++)
         {
            sm->matrix[F_Mtrx][i][j] = sm->matrix[F_Mtrx][i][j] / sm->matrix[S_Mtrx][0][j];
         }
         sm->matrix[S_Mtrx][0][j] = 0.0f;
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
sequence_matrix_calc_eeff_row_scmf (const unsigned long col,
                                    SeqMatrix* sm,
                                    const float t,
                                    float** scores)
{
   unsigned long j;
   unsigned long l;
   unsigned long k;
   short int new_matrix = F_Mtrx;
   float tmp_neg;               /* negative energy contribution */
   float tmp_het;               /* heterogenity term */
   float het_count;
   float het_rate;

   het_rate = ((-1) * ((logf (1 / 0.000001f)) / (sm->cols)));

   if (sm->curr_matrix == F_Mtrx)
   {
      new_matrix = S_Mtrx;
   }

   /*fprintf (stderr, "C: %lu Pairs with %lu\n", col, sm->pairlist[col]); */

   for (j = 0; j < sm->rows; j++)
   {
      /* update cell */
      sm->matrix[new_matrix][j][col] = 0.0f;

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
            }
            else
            {
               sm->matrix[new_matrix][j][col] +=
                  (sm->matrix[sm->curr_matrix][l][sm->pairlist[col] - 1] 
                   * scores[l][j]);
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
      tmp_het = (tmp_het / het_count) * (3.0f); /* 0.1747 */
      sm->matrix[new_matrix][j][col] += tmp_neg;
      sm->matrix[new_matrix][j][col] += tmp_het;
      sm->matrix[new_matrix][j][col] =
         expf ((-1.0f) * (sm->matrix[new_matrix][j][col]/(R_GAS * t)));
   }

   return 0;
}

/** @brief Update a row of a sequence matrix using the NN
 *
 * Update a row of a sequence matrix in a SCMF simulation.
 * Returns...
 *
 * @params[in] sm Sequence matrix.
 */
static __inline__ int
sequence_matrix_calc_eeff_row_scmf_nn (const unsigned long col,
                                       SeqMatrix* sm,
                                       const float t __attribute__((unused)),
                                       const NN_scores* scores,
                                       const Alphabet* sigma)
{
   short int new_matrix = F_Mtrx;    /* matrix to be filled */
   unsigned long j;                  /* current row */
   unsigned long l, k, m;            /* current allowed base pair */
   unsigned long bp_allowed;         /* no. of allowed base pairs */
   unsigned long alpha_size;         /* size of the alphabet */
   char bi, bj, bip1, bjm1;          /* container for bases i and j */
   float tmp_neg;                    /* negative interaction term */
   float tmp_het;                    /* heterogenity term */
   float het_rate;                   /* decrease of influence */
   float het_count;
   long G_stack_score;
   float update_prob;

   het_rate = ((-1) * ((logf (1 / 0.000001f)) / (sm->cols)));

   bp_allowed = nn_scores_no_allowed_basepairs (scores);
   alpha_size = alphabet_size (sigma);

   if (sm->curr_matrix == F_Mtrx)
   {
      new_matrix = S_Mtrx;
   }

   /* for each row */
   for (j = 0; j < sm->rows; j++)
   {
      /* update cell */
      sm->matrix[new_matrix][j][col] = 0.0f;

      /* calculate contribution of wanted interaction (if any) */
      if (sm->pairlist[col])
      {
         /* mfprintf (stderr, "Current: %lu %c\n", col, */
/*                    alphabet_no_2_base (j,sigma)); */
         if (col < sm->pairlist[col])
         {  /* we are at the "i part" of a base pair */
            /* 5' - ii+1
                    jj-1 - 3' */
            /* decide if we got to pairs or one and a mismatch */
            /* Ask Andrew: All this does not work properly for immedeate base
               pairs: ()!!! */
            /* We do not need to check here if i+1 has a base pair at all
               because sm->pairlist[col] >= 2 (pairs with pos 1),
               thus sm->pairlist[col+1] == 0 never matches anything */
            if (sm->pairlist[col+1] == (sm->pairlist[col] - 1))
            {
               /* for all allowed base pairs */
               for (l = 0; l < bp_allowed; l++)
               {
                  nn_scores_get_allowed_basepair (l, &bi, &bj, scores);
                  /* does pair match with current base? */
                  if (j == (unsigned) bi)
                  {
                     /* pair each pair with all allowed pairs */
                     for (k = 0; k < bp_allowed; k++)
                     {
                        nn_scores_get_allowed_basepair (k, &bip1, &bjm1,
                                                        scores);
                        G_stack_score = nn_scores_get_G_stack (bi, bj,
                                                               bjm1, bip1,
                                                               scores);
                        update_prob =
                 sm->matrix[sm->curr_matrix][(int)bj][sm->pairlist[col] - 1]
               * sm->matrix[sm->curr_matrix][(int)bip1][col + 1]
               * sm->matrix[sm->curr_matrix][(int)bjm1][sm->pairlist[col] - 2];
                        /*mfprintf (stderr, "  count: %c%c, %c%c; "
                                          "score: %ld; "
                                          "uprob: %f; "
                                          "cscore: %.2f\n",
                                  alphabet_no_2_base (bi, sigma),
                                  alphabet_no_2_base (bj, sigma),
                                  alphabet_no_2_base (bjm1, sigma),
                                  alphabet_no_2_base (bip1, sigma),
                                  G_stack_score,
                                  update_prob,
                                  sm->matrix[new_matrix][j][col]);*/
                        sm->matrix[new_matrix][j][col] +=
                        (update_prob * G_stack_score);
                     }
                  }
               }
            }
            else /* use mismatch stacking params */
            {
               /* for all allowed base pairs */
               for (l = 0; l < bp_allowed; l++)
               {
                  nn_scores_get_allowed_basepair (l, &bi, &bj, scores);
                  /* does pair match with current base? */
                  if (j == (unsigned) bi)
                  {
                     /* pair with each possible pair */
                     for (k = 0; k < alpha_size; k++)
                     {
                        for (m = 0; m < alpha_size; m++)
                        {
                           G_stack_score = nn_scores_get_G_mm_stack (bi, bj,
                                                                     m, k,
                                                                     scores);
                           update_prob =
                sm->matrix[sm->curr_matrix][(int)bj][sm->pairlist[col] - 1]
              * sm->matrix[sm->curr_matrix][k][col+1]
              * sm->matrix[sm->curr_matrix][m][sm->pairlist[col] - 2];
                           /*mfprintf (stderr, "  count: %c%c, %c%c; "
                                     "score: %ld; "
                                     "uprob: %f; "
                                     "cscore: %.2f\n",
                                     alphabet_no_2_base (bi, sigma),
                                     alphabet_no_2_base (bj, sigma),
                                     alphabet_no_2_base (m, sigma),
                                     alphabet_no_2_base (k, sigma),
                                     G_stack_score,
                                     update_prob,
                                     sm->matrix[new_matrix][j][col]);*/
                           sm->matrix[new_matrix][j][col] +=
                              (update_prob * G_stack_score);
                        }
                     }
                  }
               }
            }
         }
         else /* we are at the "j part" of a base pair */
         {
            /* decide if we got to pairs or one and a mismatch */
            if (sm->pairlist[col-1] == (sm->pairlist[col] + 1))
            {
               /* for all allowed base pairs */
               for (l = 0; l < bp_allowed; l++)
               {
                   nn_scores_get_allowed_basepair (l, &bi, &bj, scores);
                   /* does pair match with current base? */
                   if (j == (unsigned) bj)
                   {
                      /* pair each pair with all allowed pairs */
                      for (k = 0; k < bp_allowed; k++)
                      {
                         nn_scores_get_allowed_basepair (k, &bip1, &bjm1,
                                                         scores);
                         G_stack_score = nn_scores_get_G_stack (bi,bj,
                                                                bjm1,bip1,
                                                                scores);
                         update_prob =
                     sm->matrix[sm->curr_matrix][(int)bi][sm->pairlist[col] - 1]
                   * sm->matrix[sm->curr_matrix][(int)bip1][sm->pairlist[col]]
                   * sm->matrix[sm->curr_matrix][(int)bjm1][col - 1];
                         /*mfprintf (stderr, "  count: %c%c, %c%c "
                           "score: %ld; "
                           "uprob: %f; "
                           "cscore: %.2f\n",
                           alphabet_no_2_base (bi, sigma),
                           alphabet_no_2_base (bj, sigma),
                           alphabet_no_2_base (bjm1, sigma),
                           alphabet_no_2_base (bip1, sigma),
                           G_stack_score,
                           update_prob,
                           sm->matrix[new_matrix][j][col]);*/
                         sm->matrix[new_matrix][j][col] +=
                            (update_prob * G_stack_score);
                     }
                  }
               }
            }
            else /* use mismatch stacking params */
            {
               /* for all allowed base pairs */
               for (l = 0; l < bp_allowed; l++)
               {
                  nn_scores_get_allowed_basepair (l, &bi, &bj, scores);
                  /* does pair match with current base? */
                  if (j == (unsigned) bj)
                  {
                     /* pair with each possible pair */
                     for (k = 0; k < alpha_size; k++)
                     {
                        for (m = 0; m < alpha_size; m++)
                        {
                           G_stack_score = nn_scores_get_G_mm_stack (bi, bj,
                                                                     m, k,
                                                                     scores);
                           update_prob =
                   sm->matrix[sm->curr_matrix][(int)bi][sm->pairlist[col] - 1]
                 * sm->matrix[sm->curr_matrix][k][sm->pairlist[col]]
                 * sm->matrix[sm->curr_matrix][m][col - 1];
                           /*mfprintf (stderr, "  count: %c%c, %c%c "
                                     "score: %ld; "
                                     "uprob: %f "
                                     "cscore: %.2f\n",
                                     alphabet_no_2_base (bi, sigma),
                                     alphabet_no_2_base (bj, sigma),
                                     alphabet_no_2_base (m, sigma),
                                     alphabet_no_2_base (k, sigma),
                                     G_stack_score,
                                     update_prob,
                                     sm->matrix[new_matrix][j][col]);*/
                           sm->matrix[new_matrix][j][col] += 
                              (update_prob * G_stack_score);
                        }
                     }
                  }
               }
            }
         }
      }

      /* calculate contribution of unwanted pairs */
      tmp_neg   = 0.0f;
      tmp_het   = 0.0f;
      het_count = 0.0f;
      for (k = 0; k < sm->cols; k++)
      {
         if ((k != col) && ((k + 1) != sm->pairlist[col]))
         {
            if (col < k)        /* col is i, k is j (base pair is i < j) */
            {
               /* we only consider allowe base pairs */
               for (l = 0; l < bp_allowed; l++)
               {
                  nn_scores_get_allowed_basepair (l, &bi, &bj, scores);
                  
                  /* does pair match with current base? */
                  if (j == (unsigned) bi)
                  {
                     for (m = 0; m < bp_allowed; m++)
                     {
                       nn_scores_get_allowed_basepair (m, &bip1, &bjm1, scores);
                       update_prob =
                 sm->matrix[sm->curr_matrix][(int)bj][k]
               * sm->matrix[sm->curr_matrix][(int)bip1][col+1]
               * sm->matrix[sm->curr_matrix][(int)bjm1][k-1];
                       tmp_neg += (nn_scores_get_G_stack (bi, bj,
                                                          bjm1, bip1,
                                                          scores)
                                   * update_prob);
                       /*mfprintf (stderr, "col: %lu j: %c k: %lu bp: %c%c, "
                                 "%c%c up: %f ne: %f\n",
                                 col,
                                 alphabet_no_2_base (j, sigma),
                                 k,
                                 alphabet_no_2_base (bi, sigma),
                                 alphabet_no_2_base (bj, sigma), 
                                 alphabet_no_2_base (bip1, sigma),
                                 alphabet_no_2_base (bjm1, sigma),
                                 update_prob,
                                 tmp_neg);*/
                     }
                  }
               }
            }
            else /* k is i, col is j */
            {
               /* we only consider allowe base pairs */
               for (l = 0; l < bp_allowed; l++)
               {
                  nn_scores_get_allowed_basepair (l, &bi, &bj, scores);

                  if (j == (unsigned) bj)
                  {
                     for (m = 0; m < bp_allowed; m++)
                     {
                       nn_scores_get_allowed_basepair (m, &bip1, &bjm1, scores);
                       update_prob =
                            sm->matrix[sm->curr_matrix][(int)bi][k]
                          * sm->matrix[sm->curr_matrix][(int)bip1][k+1]
                          * sm->matrix[sm->curr_matrix][(int)bjm1][col-1];
                       tmp_neg += (nn_scores_get_G_stack (bi, bj,
                                                          bjm1, bip1,
                                                          scores)
                                   * update_prob);
                       /*mfprintf (stderr, "col: %lu j: %c k: %lu bp: %c%c, "
                                "%c%c up: %f ne: %f\n",
                                col,
                                alphabet_no_2_base (j, sigma),
                                k,
                                alphabet_no_2_base (bi, sigma),
                                alphabet_no_2_base (bj, sigma),
                                alphabet_no_2_base (bip1, sigma),
                                alphabet_no_2_base (bjm1, sigma),
                                update_prob,
                                tmp_neg);*/
                     }
                  }
               }            
            }

            /* heterogenity term */
            if (col > k)
            {
               tmp_het += (sm->matrix[sm->curr_matrix][j][k]
                           * expf(het_rate * (col - (k+1))));
               het_count += expf ( (het_rate* (col - (k+1))));
            }
            else
            {
               tmp_het += (sm->matrix[sm->curr_matrix][j][k]
                           * expf(het_rate * (k - (col+1))));
               het_count += (expf (het_rate * (k - (col+1))));
            }
         }        
      }

      /* calculate ??? */
      tmp_neg = (tmp_neg / sm->cols) * (-1.25f);
      /* mfprintf (stderr, "het: %f ", tmp_het); */
      tmp_het = (tmp_het ) * (3.0f);
      /* tmp_het = tmp_het * 10.0; */ /* 10: to large */
      /*mfprintf (stdout, "%lu:%c sm: %f neg: %f het: %f\n",
        col, alphabet_no_2_base (j, sigma), sm->matrix[new_matrix][j][col],
        tmp_neg, tmp_het);*/
      sm->matrix[new_matrix][j][col] += tmp_neg;
      sm->matrix[new_matrix][j][col] += tmp_het;
      sm->matrix[new_matrix][j][col] =
        expf ((-1.0f) * (sm->matrix[new_matrix][j][col]/(R_GAS_NN * t)));
   }

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
sequence_matrix_calc_eeff_col_scmf (SeqMatrix* sm,
                                    const float t,
                                    float** scores)
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

   if (sm->curr_matrix == F_Mtrx)
   {
      sm->curr_matrix = S_Mtrx;
   }
   else
   {
      sm->curr_matrix = F_Mtrx;
   }

   return error;
}

/** @brief Update the columns of a sequence matrix using the NN
 *
 * Update cols of a sequence matrix in a SCMF simulation.
 * Returns...
 *
 * @params[in] sm Sequence matrix.
 */
static __inline__ int
sequence_matrix_calc_eeff_col_scmf_nn (SeqMatrix* sm,
                                       const float t,
                                       const NN_scores* scores,
                                       const Alphabet* sigma)
{
   int error = 0;
   unsigned long i; 

   /* for each col */
   i = 0;

   while ((!error) && (i < sm->cols))
   {
      /* calculate all rows of a non-fixed column */
      if (!seqmatrix_is_col_fixed (i, sm))
      {
        error = sequence_matrix_calc_eeff_row_scmf_nn (i, sm, t, scores, sigma);
      }

      i++;
   }

   /* switch matrices */
   if (sm->curr_matrix == F_Mtrx)
   {
      sm->curr_matrix = S_Mtrx;
   }
   else
   {
      sm->curr_matrix = F_Mtrx;
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
seqmatrix_collate_is (const float fthresh,
                      const unsigned long steps,
                      const float temp,
                      float** scores,
                      SeqMatrix* sm,
                      const Alphabet* sigma)
{
   unsigned long i, j, /* k, */ largest_amb_col = 0, largest_amb_row = 0;
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
                  /* unambigouos site found, fixate it */
                  seqmatrix_fix_col (i, j, sm);

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
                  largest_amb_row = i;
               }
            }
         }

      }
      
      /* set site to 1/0 and fixate it */
      if (largest_amb_col < sm->cols + 1)
      {
         seqmatrix_fix_col (largest_amb_row, largest_amb_col, sm);

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

int
seqmatrix_collate_is_nn (const float fthresh,
                      const unsigned long steps,
                      const float temp,
                      const NN_scores* scores,
                      SeqMatrix* sm,
                      const Alphabet* sigma)
{
   unsigned long i, j, /* k, */ largest_amb_col = 0, largest_amb_row = 0;
   float largest_amb;
   int retval = 0;

   assert (sm);
   assert (sm->sequence == NULL);
   assert (fthresh <= 1.0f);

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
                  /* unambigouos site found, fixate it */
                  seqmatrix_fix_col (i, j, sm);

                  i = sm->rows + 1;
               }
            }
         }
      }
      
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
                  largest_amb_row = i;
               }
            }
         }
      }
      
      /* set site to 1/0 and fixate it */
      if (largest_amb_col < sm->cols + 1)
      {
         seqmatrix_fix_col (largest_amb_row, largest_amb_col, sm);
/*UGCGCACUAGGACAAUACACAGUCCAGGGGGAAAUAGACCCCCAGAUAGGGGCAUAGAAAGCCCCGUGCGCAAAGA
 */

         /* simulate */
         retval = sequence_matrix_simulate_scmf_nn (steps,
                                                    temp,
                                                    sm,
                                                    scores,
                                                    sigma);
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
seqmatrix_collate_mv (SeqMatrix* sm, const Alphabet* sigma)
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
      sm->sequence[j] = 0;
      curr_max_prob = -1.0f;
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
   float col_sum;
   float T = t_init;
   float c_rate = 1.0f;
   float c_port = 0.0f;
   short int m = 0;
   float s;

   assert (sm);
   assert (scores && scores[0]);

   if (sm->curr_matrix == F_Mtrx)
   {
     m = S_Mtrx;
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

      /* calculate Eeff */
      error = sequence_matrix_calc_eeff_col_scmf (sm, T, scores);

      /* update matrix */
      /* for all columns */
      for (j = 0; j < sm->cols; j++)
      {
         /* which are not fixed */
         if (!seqmatrix_is_col_fixed (j, sm))
         {
            /* calc sum of col */
            col_sum = 0.0f;
            for (i = 0; i < sm->rows; i++)
            {
               col_sum += sm->matrix[sm->curr_matrix][i][j];
            }

            /* for each row */
            for (i = 0; i < sm->rows; i++)
            {
               sm->matrix[sm->curr_matrix][i][j] = 
                  sm->matrix[sm->curr_matrix][i][j] / col_sum;
               /* avoid oscilation by Pnew = uPcomp + (1 - u)Pold) */
               sm->matrix[sm->curr_matrix][i][j] = 
                  (0.2 * sm->matrix[sm->curr_matrix][i][j])
                  + (0.8 * sm->matrix[m][i][j]);

               /* calculate "entropy", ignore fixed sites since ln(1) = 0 */
               s += (sm->matrix[sm->curr_matrix][i][j]
                     * logf (sm->matrix[sm->curr_matrix][i][j]));
            }
         }
      }

      s = s/sm->cols;

      t++;

      /* if (s > -0.01f) return error; */
      T = (T * c_rate) + c_port;
   }

   return error;
}

/** @brief Perform a SCMF simulation on a sequence matrix using the NN.
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
sequence_matrix_simulate_scmf_nn (const unsigned long steps,
                                  const float t_init,
                                  SeqMatrix* sm,
                                  const NN_scores* scores,
                                  const Alphabet* sigma)
{
   int error= 0;
   short int m = 0;             /* matrix to use */
   float c_rate = 1.0f;         /* cooling factor */
   float c_port = 0.0f;         /* cooling summand */
   unsigned long t;             /* time */
   unsigned long i, j;          /* iterator */
   float T = t_init;            /* current temperature */
   float s;                     /* matrix entropy */
   float col_sum;               

   assert (sm);
   assert (scores);

   if (sm->curr_matrix == F_Mtrx)
   {
     m = S_Mtrx;
   }

   /* determine cooling rate */
   if (steps > 0)
   {
      c_rate = expf ((-1) * ((logf (t_init / 1.0f)) / (steps - 1)));
   }

   /*mfprintf (stderr, "Using NN model: %lu steps, init. T = %.2f, "
             "cooling rate = %f\n", steps,
             t_init, c_rate);*/

   /* perform for a certain number of steps */
   t = 0;
   while ((!error) && (t < steps))
   {
      s = 0.0f;

      /* calculate Eeff */
      error = sequence_matrix_calc_eeff_col_scmf_nn (sm, T, scores, sigma);

      /* update matrix */
      /* for all columns */      
      for (j = 0; j < sm->cols; j++)
      {
         /* which are not fixed */
         if (!seqmatrix_is_col_fixed (j, sm))
         {
            /* calc sum of col */
            col_sum = 0.0f;
            for (i = 0; i < sm->rows; i++)
            {
               col_sum += sm->matrix[sm->curr_matrix][i][j];
            }

            /* for each row */
            for (i = 0; i < sm->rows; i++)
            {
               sm->matrix[sm->curr_matrix][i][j] = 
                  sm->matrix[sm->curr_matrix][i][j] / col_sum;
               /* avoid oscilation by Pnew = uPcomp + (1 - u)Pold) */
               sm->matrix[sm->curr_matrix][i][j] = 
                  (0.6 * sm->matrix[sm->curr_matrix][i][j])
                  + (0.4 * sm->matrix[m][i][j]);

               /* calculate "entropy", ignore fixed sites since ln(1) = 0 */
               s += (sm->matrix[sm->curr_matrix][i][j]
                     * logf (sm->matrix[sm->curr_matrix][i][j]));
            }
         }
      }

      s = (s / sm->cols) * (-1.0f);

      if (s < 0.05f) 
      {
         mfprintf (stdout, "Entropy dropout: %f\n", s);
         return error;
      }

      /* mfprintf (stdout, "%lu %f\n", t, s); */

      T = (T * c_rate) + c_port;
      t++;
   }
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

   assert (sm != NULL);
   assert (sm->sequence != NULL);

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
/*UGCGCACUAGGACAAUACACAGUCCAGGGGCAAAUAGAGCCCCAGAUAGGGGCAUAGAAAGCCCCGUGCGCAAAGA*/

/*3203211002132030030103203032032001030032132003002303201303003212300321320003*/
/*1133220001332300001133220001332100013332201033012213300031132203313322000133*/
/*3102233002233000133002233002233002130022331000331022300331002331022331020330*/
/*0320312002312003003003023023023003003023123003001320330033002132030213210300*/



/*GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA*/
/*3102233002233001233002233002231002331002331002331022310302102331022331020310*/

/*100k steps: 3102233002233000133002233002233002130022331000331022300331002331022331020330*/

/*

match with blastn, refseq_rna
UGCGCACUAGGACAAUACACAGUCCAGGGGCAAAUAGAGCCCCAGATAGGGGCATAGAAAGCCCCGTGCGCAAAGA

ACGCGTGATCCTGTTATGTGTCAGGTCCCCGTTTATCTCGGGGTCTATCCCCGTATCTTTCGGGGCACGCGTTTCT


Ref seq
GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA

*/
