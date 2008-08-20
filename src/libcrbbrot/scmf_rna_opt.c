/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbrot/scmf_rna_opt.c
 *
 *  @brief Calculate per cell effective energy for scmf.
 *
 *  Module: scmf_rna_opt
 *
 *  Library: libcrbbrot
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-08-08
 *
 *
 *  Revision History:
 *         - 2008Aug08 bienert: created
 *
 */

#include <config.h>
#include <math.h>
#include <libcrbbasic/crbbasic.h>
#include <libcrbrna/crbrna.h>
#include "seqmatrix.h"
#include "scmf_rna_opt.h"

struct Scmf_Rna_Opt_data {
      void* scores;
      Alphabet* sigma;    /* alphabet */
      char* seq;          /* sequence */
};

/** @brief Create new data object for cell energy calculations.
 *
 * Creates Scmf_Rna_Opt_data object...
 */
Scmf_Rna_Opt_data*
scmf_rna_opt_data_new (const char* file, const int line)
{
   /* allocate 1 object */
   Scmf_Rna_Opt_data* cedat = XOBJ_MALLOC(sizeof (Scmf_Rna_Opt_data), file, line);

   if (cedat != NULL)
   {
      cedat->scores = NULL;
      cedat->sigma     = NULL;
      cedat->seq       = NULL;
   }

   return cedat;
}

/** @brief Delete a cell energy data object.
 *
 */
void
scmf_rna_opt_data_delete_nn (Scmf_Rna_Opt_data* cedat)
{
   if (cedat != NULL)
   {
      alphabet_delete (cedat->sigma);
      nn_scores_delete (cedat->scores);
      XFREE (cedat->seq);
      XFREE (cedat);
   }
}

/** @brief Delete a cell energy data object.
 *
 */
void
scmf_rna_opt_data_delete_nussi (Scmf_Rna_Opt_data* cedat)
{
   if (cedat != NULL)
   {
      alphabet_delete (cedat->sigma);
      XFREE_2D ((void**)cedat->scores);
      XFREE (cedat->seq);
      XFREE (cedat);
   }
}

/** @brief Initialise a new cell energy data object and initialise it.
 */
Scmf_Rna_Opt_data*
scmf_rna_opt_data_new_nn (const unsigned long seqlen,
                           const char* file, const int line)
{
   Scmf_Rna_Opt_data* this = scmf_rna_opt_data_new (file, line);

   if (this != NULL)
   {
      this->seq = XOBJ_MALLOC ((seqlen + 1) * sizeof (*(this->seq)),
                               file, line);
      if (this->seq != NULL)
      {
         this->seq[seqlen] = '\0';
      }
      else
      {
         scmf_rna_opt_data_delete_nn (this);
         return NULL;
      }

      this->sigma = ALPHABET_NEW_SINGLE (RNA_ALPHABET, strlen(RNA_ALPHABET)/2);
      if (this->sigma == NULL)
      {
         scmf_rna_opt_data_delete_nn (this);
         return NULL;
      }

      this->scores = NN_SCORES_NEW_INIT(this->sigma);
      if (this->scores == NULL)
      {
         scmf_rna_opt_data_delete_nn (this);
         return NULL;
      }
   }

   return this;
}

/** @brief Initialise a new cell energy data object and initialise it.
 */
Scmf_Rna_Opt_data*
scmf_rna_opt_data_new_nussi (const unsigned long seqlen,
                             const char* file, const int line)
{
   Scmf_Rna_Opt_data* this = scmf_rna_opt_data_new (file, line);

   if (this != NULL)
   {
      this->seq = XOBJ_MALLOC ((seqlen + 1) * sizeof (*(this->seq)),
                               file, line);
      if (this->seq != NULL)
      {
         this->seq[seqlen] = '\0';
      }
      else
      {
         scmf_rna_opt_data_delete_nussi (this);
         return NULL;
      }

      this->sigma = ALPHABET_NEW_SINGLE (RNA_ALPHABET, strlen(RNA_ALPHABET)/2);
      if (this->sigma == NULL)
      {
         scmf_rna_opt_data_delete_nussi (this);
         return NULL;
      }

      this->scores = create_scoring_matrix (this->sigma);
      if (this->scores == NULL)
      {
         scmf_rna_opt_data_delete_nussi (this);
         return NULL;
      }
   }

   return this;
}

/** @brief Set scores for the nearest neighbour model and an alphabet.
 *
 */
void
scmf_rna_opt_data_set (NN_scores* scores, Alphabet* sigma,
                      Scmf_Rna_Opt_data* cedat)
{
   assert (cedat);
   assert (scores);
   assert (sigma);

   cedat->scores = scores;
   cedat->sigma = sigma;
}



int
scmf_rna_opt_data_transform_row_2_base (const unsigned long row,
                                        const unsigned long col,
                                        void* data)
{
   Scmf_Rna_Opt_data* cont;
   assert (data);

   cont = (Scmf_Rna_Opt_data*) data;

   cont->seq[col] = alphabet_no_2_base(row, cont->sigma);

   return 0;
}

/** @brief Retriev alphabet component from data container.
 */
/* Alphabet* */
/* scmf_rna_opt_data_get_alphabet (Scmf_Rna_Opt_data* cont) */
/* { */
/*    assert (cont->sigma); */

/*    return cont->sigma; */
/* } */

/** @brief calculate energy using the nussinov model.
 *
 * Calculate the energy for a cell of a sequence matrix using the nussinov
 * energy model. Supposed to be placed in the inner most loop.
 * The energy has mainly three terms: Interaction energy (actually no. of
 * Hbonds), negative design energy (negative interaction energy to all possible
 * partners) and a heterogenity term.\n
 * Returns the energy value for the cell.
 *
 * @params[in] row Current row.
 * @params[in] col Current column.
 * @params[in] sco Scoring matrix.
 * @params[in] sm Sequence matrix.
 */
float
scmf_rna_opt_calc_nussinov (const unsigned long row,
                           const unsigned long col,
                           void* sco,
                           SeqMatrix* sm)
{
   float cell = 0.0f;             /* cell to be calculated */
   unsigned long interaction = 0; /* col of cell interacts with col of sm */
   unsigned long i, j;            /* indeces */
   unsigned long rows;            /* no. of rows of the matrix */
   unsigned long cols;            /* no. of cols of the matrix */
   float tmp_neg = 0.0f;          /* tmp for negative design term */
   float tmp_het = 0.0f;          /* tmp for heterogenity energy */
   float het_count = 0.0f;        /* heterogenity normalisation factor */
   float het_rate;
   float** scores;
   Scmf_Rna_Opt_data* test;

   assert (sm);
   assert (row < seqmatrix_get_rows (sm));
   assert (col < seqmatrix_get_width (sm));
   assert (sco);

   test = (Scmf_Rna_Opt_data*) sco;
   scores = (float**) test->scores;

   rows = seqmatrix_get_rows (sm);
   cols = seqmatrix_get_width (sm);

   het_rate = ((-1) * ((logf (1 / 0.000001f)) / (cols)));

   /* calculate contribution of wanted interaction (if any) */
   interaction = seqmatrix_col_interacts_with (col, sm);
   if (interaction)
   {
      for (i = 0; i < rows; i++)
      {
         if (col < interaction)
         {
            cell += (seqmatrix_get_probability (i, (interaction - 1), sm)
                     * scores[row][i]);
         }
         else
         {
            cell += (seqmatrix_get_probability (i, (interaction - 1), sm)
                     * scores[i][row]);
         }
      }
   }

   /* calculate contribution of unwanted pairs */
   for (j = 0; j < cols; j++)
   {
      if ((j != col) && ((j + 1) != interaction))
      {
         for (i = 0; i < rows; i++)
         {
            if (col < j)
            {
               tmp_neg += (seqmatrix_get_probability (i, j, sm) 
                           * scores[row][i]);
            }
            else
            {
               tmp_neg += (seqmatrix_get_probability (i, j, sm) 
                           * scores[i][row]);
            }
         }
            
         /* heterogenity term */
         if (col > j)
         {
            tmp_het += (seqmatrix_get_probability (row, j, sm)
                        * expf(het_rate * (col - (j+1))));
            het_count += expf ( (het_rate* (col - (j+1))));
         }
         else
         {
            tmp_het += (seqmatrix_get_probability (row, j, sm)
                        * expf(het_rate * (j - (col+1))));
            het_count += (expf (het_rate * (j - (col+1))));
         }
      }
   }

   tmp_neg = (tmp_neg / cols) * (-1.25f); /* 1.25 */
   tmp_het = (tmp_het / het_count) * (3.0f); /* 0.1747 */
   cell += tmp_neg;
   cell += tmp_het;
   /*cell =  expf ((-1.0f) * (cell/(sm->gas_constant*t)));*/
   
   /* seqmatrix_set_cell (cell, row, col, sm); */
   
   return cell;
}

char*
scmf_rna_opt_data_get_seq (Scmf_Rna_Opt_data* this)
{
   return this->seq;
}

/** @brief calculate energy using the Nearest Neighbour energy model.
 *
 * Calculate the energy for a cell of a sequence matrix using the Nearest
 * Neighbour scoring scheme. Supposed to be placed in the inner most loop.
 * The energy has mainly three terms: Interaction energy, negative design
 * energy (negative interaction energy to all possible partners) and a
 * heterogenity term.\n
 * Returns the energy vaslue for the cell.
 *
 * @params[in] row Current row.
 * @params[in] col Current column.
 * @params[in] sco Nearest Neighbour scores.
 * @params[in] sm Sequence matrix.
 */
float
scmf_rna_opt_calc_nn (const unsigned long row,
                     const unsigned long col,
                     void* sco,
                     SeqMatrix* sm)
{
   float cell = 0;                /* cell to be calculated */
   unsigned long cols;            /* no. of cols of the matrix */
   unsigned long bp_allowed;      /* no. of allowed base pairs */
   unsigned long alpha_size;      /* size of the alphabet */
   unsigned long interaction = 0; /* col of cell interacts with col of sm */
   unsigned long l, k, m;         /* indeces */
   char bi, bj, bip1, bjm1;       /* container for bases i and j */
   long G_stack_score;            /* fetch Gibb's free energy for a stack */
   float update_prob;             /* probability component for cell energy */
   float tmp_neg = 0.0f;          /* negative interaction term */
   float tmp_het = 0.0f;          /* heterogenity term */
   float het_count = 0.0f;
   float het_rate;
   Scmf_Rna_Opt_data* cedat;

   assert (sm);
   assert (row < seqmatrix_get_rows (sm));
   assert (col < seqmatrix_get_width (sm));
   assert (sco);

   cedat = (Scmf_Rna_Opt_data*) sco;

   cols = seqmatrix_get_width (sm);
   bp_allowed = nn_scores_no_allowed_basepairs (cedat->scores);
   alpha_size = alphabet_size (cedat->sigma);

   het_rate = ((-1) * ((logf (1 / 0.000001f)) / (cols)));

   /* calculate contribution of wanted interaction (if any) */
   interaction = seqmatrix_col_interacts_with (col, sm);
   if (interaction)
   {
      if (col < interaction)
      {  /* we are at the "i part" of a base pair */
         /* 5' - ii+1
                 jj-1 - 3' */
         /* decide if we got to pairs or one and a mismatch */
         /* Ask Andrew: All this does not work properly for immedeate base
            pairs: ()!!! */
         /* We do not need to check here if i+1 has a base pair at all because
            seqmatrix_col_interacts_with (col, sm) >= 2 (pairs with pos 1),
            thus seqmatrix_col_interacts_with (col+1, sm) == 0 never matches
            anything */
         if (seqmatrix_col_interacts_with ((col + 1), sm) == (interaction - 1))
         {
            /* for all allowed base pairs */
            for (l = 0; l < bp_allowed; l++)
            {
               nn_scores_get_allowed_basepair (l, &bi, &bj, cedat->scores);

               /* does pair match with current base? */
               if (row == (unsigned) bi)
               {
                  /* pair each pair with all allowed pairs */
                  for (k = 0; k < bp_allowed; k++)
                  {
                     nn_scores_get_allowed_basepair (k, &bip1, &bjm1,
                                                     cedat->scores);
                     G_stack_score = nn_scores_get_G_stack (bi, bj, bjm1, bip1,
                                                            cedat->scores);
                     update_prob = 
                          seqmatrix_get_probability (bj, (interaction - 1), sm)
                        * seqmatrix_get_probability (bip1, (col + 1), sm)
                        * seqmatrix_get_probability (bjm1,(interaction - 2),sm);

                     cell += (update_prob * G_stack_score);
                  }
               }
            }
         }
         else /* use mismatch stacking params */
         {
            /* for all allowed base pairs */
            for (l = 0; l < bp_allowed; l++)
            {
               nn_scores_get_allowed_basepair (l, &bi, &bj, cedat->scores);
               
               /* does pair match with current base? */
               if (row == (unsigned) bi)
               {
                  for (k = 0; k < alpha_size; k++)
                  {
                     for (m = 0; m < alpha_size; m++)
                     {
                        G_stack_score = nn_scores_get_G_mm_stack (bi, bj, m, k,
                                                              cedat->scores);
                        update_prob =
                           seqmatrix_get_probability (bj, (interaction - 1), sm)
                         * seqmatrix_get_probability (k, (col + 1), sm)
                         * seqmatrix_get_probability (m, (interaction - 2), sm);

                        cell += (update_prob * G_stack_score);
                     }
                  }
               }
            }
         }
      }
      else /* we are at the "j part" of a base pair */
      {
         /* decide if we got to pairs or one and a mismatch */
         if (seqmatrix_col_interacts_with ((col - 1), sm) == (interaction + 1))
         {
            /* for all allowed base pairs */
            for (l = 0; l < bp_allowed; l++)
            {
               nn_scores_get_allowed_basepair (l, &bi, &bj, cedat->scores);

               /* does pair match with current base? */
               if (row == (unsigned) bj)
               {
                  /* pair each pair with all allowed pairs */
                  for (k = 0; k < bp_allowed; k++)
                  {
                     nn_scores_get_allowed_basepair (k, &bip1, &bjm1,
                                                     cedat->scores);

                     G_stack_score = nn_scores_get_G_stack (bi,bj, bjm1,bip1,
                                                            cedat->scores);

                     update_prob = 
                          seqmatrix_get_probability (bi, (interaction - 1), sm)
                        * seqmatrix_get_probability (bip1, interaction, sm)
                        * seqmatrix_get_probability (bjm1, (col - 1), sm);

                     cell += (update_prob * G_stack_score);                 
                  }
               }
            }
         }
         else /* use mismatch stacking params */
         {
            /* for all allowed base pairs */
            for (l = 0; l < bp_allowed; l++)
            {
               nn_scores_get_allowed_basepair (l, &bi, &bj, cedat->scores);

               /* does pair match with current base? */
               if (row == (unsigned) bj)
               {
                  /* pair with each possible pair */
                  for (k = 0; k < alpha_size; k++)
                  {
                     for (m = 0; m < alpha_size; m++)
                     {
                        G_stack_score = nn_scores_get_G_mm_stack (bi, bj,  m, k,
                                                              cedat->scores);

                        update_prob =
                          seqmatrix_get_probability (bi, (interaction - 1), sm)
                        * seqmatrix_get_probability (k, interaction, sm)
                        * seqmatrix_get_probability (m, (col - 1), sm);

                        cell += (update_prob * G_stack_score);
                     }
                  }
               }
            }
         }
      }
   }

   /* calculate contribution of unwanted pairs */
   for (k = 0; k < cols; k++)
   {
      if ((k != col) && ((k + 1) != interaction))
      {
         if (col < k) /* col is i, k is j (base pair is i < j) */
         {
            /* we only consider allowed base pairs */
            for (l = 0; l < bp_allowed; l++)
            {
               nn_scores_get_allowed_basepair (l, &bi, &bj, cedat->scores);

               /* does pair match with current base? */
               if (row == (unsigned) bi)
               {
                  for (m = 0; m < bp_allowed; m++)
                  {
                     nn_scores_get_allowed_basepair (m, &bip1, &bjm1,
                                                     cedat->scores);

                     update_prob =
                          seqmatrix_get_probability (bj, k, sm)
                        * seqmatrix_get_probability (bip1, (col + 1), sm)
                        * seqmatrix_get_probability (bjm1, (k - 1), sm);

                     tmp_neg += (nn_scores_get_G_stack (bi, bj, bjm1, bip1,
                                                        cedat->scores)
                                 * update_prob);
                  }
               }               
            }
         }
         else         /* k is i, col is j */
         {
            /* we only consider allowe base pairs */
            for (l = 0; l < bp_allowed; l++)
            {
               nn_scores_get_allowed_basepair (l, &bi, &bj, cedat->scores);
               
               if (row == (unsigned) bj)
               {
                  for (m = 0; m < bp_allowed; m++)
                  {
                     nn_scores_get_allowed_basepair (m, &bip1, &bjm1,
                                                     cedat->scores);

                     update_prob = 
                          seqmatrix_get_probability (bi, k, sm)
                        * seqmatrix_get_probability (bip1, (k + 1), sm)
                        * seqmatrix_get_probability (bjm1, (col - 1), sm);

                     tmp_neg += (nn_scores_get_G_stack (bi, bj, bjm1, bip1,
                                                        cedat->scores)
                                 * update_prob);
                  }
               }
            }
         }

         /* heterogenity term */
         if (col > k)
         {
            tmp_het += (seqmatrix_get_probability (row, k, sm)
                        * expf(het_rate * (col - (k+1))));
            het_count += expf ( (het_rate* (col - (k+1))));
         }
         else
         {
            tmp_het += (seqmatrix_get_probability (row, k, sm)
                        * expf(het_rate * (k - (col+1))));
            het_count += (expf (het_rate * (k - (col+1))));
         }         
      }
   }

   tmp_neg = (tmp_neg / cols) * (-1.25f);
   tmp_het = (tmp_het ) * (3.0f);
   cell += tmp_neg;
   cell += tmp_het;

   return cell;
}
