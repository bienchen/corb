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
      Rna* rna;           /* sequence and pairlist */
      char** bp_allowed;  /* mapping of base to allowed pairing partners */
      float het_rate;     /* rate of decrease for the heterogenity term */
      float* en_neg;      /* negative energies for incremental updates */
      float** en_neg2;    /* negative energies for incremental updates */
      float** en_neg_35;  /* neg.en. for 5' 3' direction */
};

/** @brief Create new data object for cell energy calculations.
 *
 * Creates Scmf_Rna_Opt_data object...
 */
Scmf_Rna_Opt_data*
scmf_rna_opt_data_new (const char* file, const int line)
{
   /* allocate 1 object */
   Scmf_Rna_Opt_data* cedat = XOBJ_MALLOC(sizeof (Scmf_Rna_Opt_data),
                                          file, line);

   if (cedat != NULL)
   {
      cedat->scores     = NULL;
      cedat->sigma      = NULL;
      cedat->rna        = NULL;
      cedat->bp_allowed = NULL;
      cedat->het_rate   = 0.0f;
      cedat->en_neg     = NULL;
      cedat->en_neg2    = NULL;
      cedat->en_neg_35   = NULL;
   }

   return cedat;
}

void
scmf_rna_opt_data_delete (Scmf_Rna_Opt_data* cedat)
{
   if (cedat != NULL)
   {
      assert (cedat->scores == NULL);
      assert (cedat->bp_allowed == NULL);

      alphabet_delete (cedat->sigma);
      rna_delete (cedat->rna);
      XFREE (cedat->en_neg);
      XFREE_2D ((void**)cedat->en_neg2);
      XFREE_2D ((void**)cedat->en_neg_35);
      XFREE (cedat);
   }
}

/** @brief Initialise a new cell energy data object.
 */
Scmf_Rna_Opt_data*
scmf_rna_opt_data_new_init (const char* structure,
                            const unsigned long seqlen,
                            const char* alpha_string,
                            const unsigned long alpha_size,
                            float het_rate,
                            const char* file, const int line)
{
   int error;

   Scmf_Rna_Opt_data* this = scmf_rna_opt_data_new (file, line);

   if (this != NULL)
   {
      this->het_rate = het_rate;

      /* init alphabet */
      this->sigma = ALPHABET_NEW_SINGLE (alpha_string, alpha_size);
      if (this->sigma == NULL)
      {
         scmf_rna_opt_data_delete (this);
         return NULL;
      }

      /* init Rna object */
      this->rna = RNA_NEW;
      if (this->rna == NULL)
      {
         scmf_rna_opt_data_delete (this);
         return NULL;
      }
      error = RNA_ALLOC_SEQUENCE(seqlen, this->rna);
      if (error)
      {
         scmf_rna_opt_data_delete (this);
         return NULL;         
      }
      
      /* init negative energies */
      this->en_neg = XCALLOC (alphabet_size (this->sigma),
                              sizeof (*(this->en_neg)));

      this->en_neg2 = (float**) XMALLOC_2D (alphabet_size (this->sigma),
                                           alphabet_size (this->sigma),
                                             sizeof (float));
      this->en_neg_35 = (float**) XMALLOC_2D (alphabet_size (this->sigma),
                                           alphabet_size (this->sigma),
                                             sizeof (float));

      error = RNA_INIT_PAIRLIST_VIENNA(structure, seqlen, this->rna);
   }

   return this;
}

/* int */
/* scmf_rna_opt_data_init_negative_design_energies (void* data, */
/*                                                  SeqMatrix* sm) */
/* { */
/*    unsigned long i, j, k, l, m, cols, alpha, allowed_bp; */
/*    char bj, bip1, bjm1; */
/*    float prob; */
/*    Scmf_Rna_Opt_data* this; */

/*    assert (sm); */
/*    assert (data); */

/*    this = (Scmf_Rna_Opt_data*) data; */

/*    cols = seqmatrix_get_width (sm); */
/*    alpha = alphabet_size (this->sigma); */
/*    allowed_bp = nn_scores_no_allowed_basepairs (this->scores); */
/*    i = 0; */

/*    for (k = 0; k < alpha; k++) */
/*    { */
/*       this->en_neg[k] = 0; */
/*    } */

/*    /\*memset(this->en_neg, 0, alpha * sizeof (*(this->en_neg)));*\/ */

/*    /\* idea: init neg. energies as neg. energies for the first col *\/ */
/*    /\* update during matrix calculation *\/ */
/*    /\* what about interactions? *\/ */
/*    for (j = 1; j < cols; j++) */
/*    { */
/*       for (k = 0; k < alpha; k++) */
/*       { */
/*          /\* we only consider allowed base pairs *\/ */
/*          for (l = 0; this->bp_allowed[k][l] != 0; l++) */
/*          { */
/*             bj = this->bp_allowed[k][l] - 1;  */
            
/*             for (m = 0; m < allowed_bp; m++) */
/*             { */
/*                nn_scores_get_allowed_basepair (m, &bip1, &bjm1, */
/*                                                this->scores); */
               
/*                prob = seqmatrix_get_probability (bj, j, sm) */
/*                   * seqmatrix_get_probability (bip1, 1, sm) */
/*                   * seqmatrix_get_probability (bjm1, (j - 1), sm); */
               
/*                this->en_neg[k] += (nn_scores_get_G_stack (k, bj, bjm1, bip1, */
/*                                                           this->scores) */
/*                                    * prob); */
/*             } */
/*          } */
/*       } */
/*    } */

/*    mfprintf (stderr, "E: %f\n", this->en_neg[0]); */
/*    mfprintf (stderr, "E: %f\n", this->en_neg[1]); */
/*    mfprintf (stderr, "E: %f\n", this->en_neg[2]); */
/*    mfprintf (stderr, "E: %f\n", this->en_neg[3]); */

/*    return 0; */
/* } */

int
scmf_rna_opt_data_init_negative_design_energies_alt (void* data,
                                                 SeqMatrix* sm)
{
   unsigned long i, j, k, l, m, cols, alpha, allowed_bp;
   char bj, bip1, bjm1;
   float prob;
   Scmf_Rna_Opt_data* this;

   assert (sm);
   assert (data);

   this = (Scmf_Rna_Opt_data*) data;

   cols = seqmatrix_get_width (sm);
   alpha = alphabet_size (this->sigma);
   allowed_bp = nn_scores_no_allowed_basepairs (this->scores);
   i = 0;

   memset(this->en_neg2[0], 0, (alpha * alpha) * sizeof (*(this->en_neg2)));
   memset(this->en_neg_35[0], 0, (alpha * alpha) * sizeof (*(this->en_neg_35)));
   memset(this->en_neg, 0, alpha * sizeof (*(this->en_neg)));

   /*for (k = 0; k < alpha; k++)
   {
      for (i = 0; i < alpha; i++)
         mfprintf (stderr, "%lu: %f ", k, this->en_neg2[k][i]);
      mfprintf (stderr, "\n");
      }*/

   /* idea: calculate all neg. contributions for the first column without
            incorporating the probabilities from this column into the
            calculation. The probabilities are then used during cell
            evaluation. */
   for (j = 1; j < cols; j++)
   {
      for (k = 0; k < alpha; k++)
      {
         /* we only consider allowed base pairs */
         for (l = 0; this->bp_allowed[k][l] != 0; l++)
         {
            bj = this->bp_allowed[k][l] - 1;
            
            for (m = 0; m < allowed_bp; m++)
            {
               nn_scores_get_allowed_basepair (m, &bip1, &bjm1,
                                               this->scores);

               prob = seqmatrix_get_probability (bj,    j,      sm)
                    * seqmatrix_get_probability (bjm1, (j - 1), sm);
        
               this->en_neg2[k][(int)bip1] += (nn_scores_get_G_stack (k, bj,
                                                                 bjm1, bip1,
                                                                 this->scores)
                                          * prob);
               this->en_neg[k] += (nn_scores_get_G_stack (k, bj,
                                                                 bjm1, bip1,
                                                                 this->scores)
                                          * prob);
            }
         }
      }
   }

   /*mfprintf (stderr, "E: %f\n", this->en_neg2[0][0]
             +this->en_neg2[0][1]
             +this->en_neg2[0][2]
             +this->en_neg2[0][3]);
   mfprintf (stderr, "E: %f\n", this->en_neg[0]);
   mfprintf (stderr, "E: %f\n", this->en_neg2[1][0]
             +this->en_neg2[1][1]
             +this->en_neg2[1][2]
             +this->en_neg2[1][3]);
   mfprintf (stderr, "E: %f\n", this->en_neg[1]);
   mfprintf (stderr, "E: %f\n", this->en_neg2[2][0]
             +this->en_neg2[2][1]
             +this->en_neg2[2][2]
             +this->en_neg2[2][3]);
   mfprintf (stderr, "E: %f\n", this->en_neg[2]);
   mfprintf (stderr, "E: %f\n", this->en_neg2[3][0]
             +this->en_neg2[3][1]
             +this->en_neg2[3][2]
             +this->en_neg2[3][3]);
             mfprintf (stderr, "E: %f\n", this->en_neg[3]);*/
   return 0;
}

static __inline__ void
scmf_rna_opt_iterate_neg_design_term (unsigned long row,
                                      unsigned long col,
                                      Scmf_Rna_Opt_data* this,
                                      SeqMatrix* sm)
{
   unsigned long j, k;
   char bi, bj, bip1, bjm1;
   unsigned long allowed_bp;
   float prob;

   assert (this);
   assert (sm);

   allowed_bp = nn_scores_no_allowed_basepairs (this->scores);

   if ((col + 1) < seqmatrix_get_width(sm))
   {
      /* update negative energy for nex col */
      /* subtract col, col+1 contribution (5' -> 3') */
      for (j = 0; this->bp_allowed[row][j] != 0; j++)
      {
         bj = this->bp_allowed[row][j] - 1;
         
         for (k = 0; k < allowed_bp; k++)
         {
            nn_scores_get_allowed_basepair (k, &bip1, &bjm1,
                                            this->scores);
            
            prob = seqmatrix_get_probability (bj, (col + 1), sm)
               * seqmatrix_get_probability (bjm1, col, sm);
            
            this->en_neg2[row][(int)bip1] -= (nn_scores_get_G_stack (row, bj,
                                                                     bjm1, bip1,
                                                                     this->scores)
                                              * prob);
         }
      }

      /* update negative energy for next col */
      /* add col+1, col contribution (5' <- 3') */
      for (j = 0; this->bp_allowed[row][j] != 0; j++)
      {
         bi = this->bp_allowed[row][j] - 1;
         
         for (k = 0; k < allowed_bp; k++)
         {
            nn_scores_get_allowed_basepair (k, &bip1, &bjm1,
                                            this->scores);
            
            prob = seqmatrix_get_probability (bi, col, sm)
               * seqmatrix_get_probability (bip1, (col + 1), sm);
            
            this->en_neg_35[row][(int)bjm1] += (nn_scores_get_G_stack (bi, row,
                                                                       bjm1, bip1,
                                                                       this->scores)
                                                * prob);
         }
      }
   }
}

int
scmf_rna_opt_data_update_neg_design_energy (void* data,
                                            unsigned long col,
                                            SeqMatrix* sm)
{
   Scmf_Rna_Opt_data* this;
   unsigned long rows, i;

   assert (sm);
   assert (data);
   
   rows = seqmatrix_get_rows (sm);

   this = (Scmf_Rna_Opt_data*) data;

   for (i = 0; i < rows; i++)
   {
      scmf_rna_opt_iterate_neg_design_term (i, col, this, sm);
   }

   return 0;
}

void
scmf_rna_opt_data_set_scores (void* scores,
                              Scmf_Rna_Opt_data* cedat)
{
   assert (cedat);

   cedat->scores = scores;
}

void
scmf_rna_opt_data_set_bp_allowed (char** bp_allowed,
                                  Scmf_Rna_Opt_data* cedat)
{
   assert (cedat);

   cedat->bp_allowed = bp_allowed;
}

int
scmf_rna_opt_data_transform_row_2_base (const unsigned long row,
                                        const unsigned long col,
                                        void* data)
{
   Scmf_Rna_Opt_data* cont;
   assert (data);

   cont = (Scmf_Rna_Opt_data*) data;

   rna_set_sequence_base (alphabet_no_2_base(row, cont->sigma), col, cont->rna);

   return 0;
}

/** @brief Retriev alphabet component from data container.
 */
Alphabet*
scmf_rna_opt_data_get_alphabet (Scmf_Rna_Opt_data* cont)
{
   assert (cont);
   assert (cont->sigma);

   return cont->sigma;
}

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
   unsigned long interaction;     /* col of cell interacts with col of sm */
   unsigned long i, j;            /* indeces */
   unsigned long rows;            /* no. of rows of the matrix */
   unsigned long cols;            /* no. of cols of the matrix */
   float tmp_neg = 0.0f;          /* tmp for negative design term */
   float tmp_het = 0.0f;          /* tmp for heterogenity energy */
   float het_count = 0.0f;        /* heterogenity normalisation factor */
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

   /* calculate contribution of wanted interaction (if any) */
   /* interaction = seqmatrix_col_interacts_with (col, sm);*/
   interaction = rna_base_pairs_with (col, test->rna);
   if (interaction != NOT_PAIRED)
   {
      for (i = 0; i < rows; i++)
      {
         if (col < interaction)
         {
            cell += (seqmatrix_get_probability (i, interaction, sm)
                     * scores[row][i]);
         }
         else
         {
            cell += (seqmatrix_get_probability (i, interaction, sm)
                     * scores[i][row]);
         }
      }
   }

   /* calculate contribution of unwanted pairs */
   for (j = 0; j < cols; j++)
   {
      if ((j != col) && (j != interaction))
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
                        * expf(test->het_rate * (col - (j+1))));
            het_count += expf ( (test->het_rate* (col - (j+1))));
         }
         else
         {
            tmp_het += (seqmatrix_get_probability (row, j, sm)
                        * expf(test->het_rate * (j - (col+1))));
            het_count += (expf (test->het_rate * (j - (col+1))));
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
   /* return this->seq */
   return rna_get_sequence (this->rna);
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
   unsigned long allowed_bp;      /* no. of allowed base pairs */
   unsigned long alpha_size;      /* size of the alphabet */
   unsigned long interaction;     /* col of cell interacts with col of sm */
   unsigned long l, k, m;         /* indeces */
   char bi, bj, bip1, bjm1;       /* container for bases i and j */
   long G_stack_score;            /* fetch Gibb's free energy for a stack */
   float update_prob;             /* probability component for cell energy */
   float tmp_neg = 0.0f;          /* negative interaction term */
   float tmp_het = 0.0f;          /* heterogenity term */
   float het_count = 0.0f;
   Scmf_Rna_Opt_data* cedat;
/*    long tmp; */
/*    long tmp4; */
/*    float tmp2 = 0; */
/*    float tmp3 = 0; */
/*    float tmp4; */

   assert (sm);
   assert (row < seqmatrix_get_rows (sm));
   assert (col < seqmatrix_get_width (sm));
   assert (sco);

   cedat = (Scmf_Rna_Opt_data*) sco;

   cols = seqmatrix_get_width (sm);
   allowed_bp = nn_scores_no_allowed_basepairs (cedat->scores);
   alpha_size = alphabet_size (cedat->sigma);

   /* calculate contribution of wanted interaction (if any) */
   interaction = rna_base_pairs_with (col, cedat->rna);

   if ((interaction != NOT_PAIRED) && (col < interaction))
   {  /* we are at the "i part" of a base pair */
      /* 5' - ii+1
         jj-1 - 3' */
      /* decide if we got two pairs or one and a mismatch */
      /* All this does not work properly for immedeate base pairs: ()!!! */
      if (rna_base_pairs_with ((col + 1), cedat->rna) == interaction - 1)
      {
         /* for all allowed base pairs */
         for (l = 0; cedat->bp_allowed[row][l] != 0; l++)
         {
            bj = cedat->bp_allowed[row][l] - 1;
            
            /* pair each pair with all allowed pairs */
            for (k = 0; k < allowed_bp; k++)
            {
               nn_scores_get_allowed_basepair (k, &bip1, &bjm1, cedat->scores);
   
               G_stack_score = nn_scores_get_G_stack (row, bj, bjm1, bip1,
                                                      cedat->scores);
  
               update_prob = 
                  seqmatrix_get_probability (bj, interaction, sm)
                  * seqmatrix_get_probability (bip1, (col + 1), sm)
                  * seqmatrix_get_probability (bjm1,(interaction - 1),sm);
   
               cell += (update_prob * G_stack_score);
            }
         }
      }
      else /* use mismatch stacking params */
      {
         /* for all allowed base pairs */
         for (l = 0; cedat->bp_allowed[row][l] != 0; l++)
         {
            bj = cedat->bp_allowed[row][l] - 1;            
            
            for (k = 0; k < alpha_size; k++)
            {
               for (m = 0; m < alpha_size; m++)
               {
                  G_stack_score = nn_scores_get_G_mm_stack (row, bj, m, k,
                                                            cedat->scores);
                  update_prob =
                     seqmatrix_get_probability (bj, interaction, sm)
                     * seqmatrix_get_probability (k, (col + 1), sm)
                     * seqmatrix_get_probability (m, (interaction - 1), sm);
                  
                  cell += (update_prob * G_stack_score);

               }
            }
         }
      }
   }
   else if (interaction != NOT_PAIRED) /*we are at the "j part" of a base pair*/
   {
      /* decide if we got two pairs or one and a mismatch */
      if (rna_base_pairs_with ((col - 1), cedat->rna) == (interaction + 1))
      {
         /* for all allowed base pairs */
         for (l = 0; cedat->bp_allowed[row][l] != 0; l++)
         {
            bi = cedat->bp_allowed[row][l] - 1; 
            
            /* pair each pair with all allowed pairs */
            for (k = 0; k < allowed_bp; k++)
            {
               nn_scores_get_allowed_basepair (k, &bip1, &bjm1,
                                               cedat->scores);
               
               G_stack_score = nn_scores_get_G_stack (bi, row, bjm1,bip1,
                                                      cedat->scores);
               
               update_prob = 
                  seqmatrix_get_probability (bi, interaction, sm)
                  * seqmatrix_get_probability (bip1, (interaction + 1), sm)
                  * seqmatrix_get_probability (bjm1, (col - 1), sm);
               
               cell += (update_prob * G_stack_score);                 
            }
         }
      }
      else /* use mismatch stacking params */
      {
         /* for all allowed base pairs */
         for (l = 0; cedat->bp_allowed[row][l] != 0; l++)
         {
            bi = cedat->bp_allowed[row][l] - 1; 
            
            /* pair with each possible pair */
            for (k = 0; k < alpha_size; k++)
            {
               for (m = 0; m < alpha_size; m++)
               {
                  G_stack_score = nn_scores_get_G_mm_stack (bi, row,  m, k,
                                                            cedat->scores);
                  
                  update_prob =
                     seqmatrix_get_probability (bi, interaction, sm)
                     * seqmatrix_get_probability (k, (interaction + 1), sm)
                     * seqmatrix_get_probability (m, (col - 1), sm);
                  
                  cell += (update_prob * G_stack_score);
               }
            }
         }
      }
   }

   /* SB NEW! 08-09-10 */
   /* get neg.energy for current row n' col */
   tmp_neg = 0;
   /* SB: Cleanup 08-09-12 */
   if (col > 0)
   {
      /* handle neg.interactions upstream */
      for (k = 0; k < alpha_size; k++)
      {
         tmp_neg += (  cedat->en_neg_35[row][k]
                       * seqmatrix_get_probability (k, (col - 1), sm));
      }
   }
   if ((col + 1) < cols)
   {
      /* handle downstream interactions */
      for (k = 0; k < alpha_size; k++)
      {
         tmp_neg += (  cedat->en_neg2[row][k]
                       * seqmatrix_get_probability (k, (col + 1), sm));
      }      
   }

   /* care about interactions */
   /* interacts with downstream position */
   if ((interaction != NOT_PAIRED) && (col < interaction))
   {
      for (l = 0; cedat->bp_allowed[row][l] != 0; l++)
      {
         bj = cedat->bp_allowed[row][l] - 1;
         
         for (k = 0; k < allowed_bp; k++)
         {
            nn_scores_get_allowed_basepair (k, &bip1, &bjm1,
                                            cedat->scores);
            
            G_stack_score = nn_scores_get_G_stack (row, bj, bjm1, bip1,
                                                   cedat->scores);
            update_prob =
               seqmatrix_get_probability (bj, interaction, sm)
               * seqmatrix_get_probability (bip1, (col + 1), sm)
               * seqmatrix_get_probability (bjm1,(interaction - 1),sm);
            
            tmp_neg -= (update_prob * G_stack_score);
         }
      }
   }
   /* interacts with upstream position */
   else if (interaction != NOT_PAIRED)
   {
      for (l = 0; cedat->bp_allowed[row][l] != 0; l++)
      {
         bi = cedat->bp_allowed[row][l] - 1;
         
         for (k = 0; k < allowed_bp; k++)
         {
            nn_scores_get_allowed_basepair (k, &bip1, &bjm1,
                                            cedat->scores);
            
            G_stack_score = nn_scores_get_G_stack (bi, row, bjm1,bip1,
                                                   cedat->scores);
            
            update_prob =
               seqmatrix_get_probability (bi, interaction, sm)
               * seqmatrix_get_probability (bip1, (interaction + 1), sm)
               * seqmatrix_get_probability (bjm1, (col - 1), sm);
            
            tmp_neg -= (update_prob * G_stack_score);
         }
      }
   }   
   /* SB: 08-09-12 */
   
   scmf_rna_opt_iterate_neg_design_term (row, col, cedat, sm);

   /* SB END - 08-09-10 */


   /* calculate contribution of unwanted pairs */
   /*tmp = 0; tmp2 = 0.0f; tmp3 = 0.0f; */
   for (k = 0; k < cols; k++)
   {
      if ((k != col) && (k != interaction))
      {
/*SNIP*/
/*          if (col < k) /\* col is i, k is j (base pair is i < j) *\/ */
/*          { */
/*             /\* we only consider allowed base pairs *\/ */
/*             /\* tmp = 0; *\/ */
/*             for (l = 0; cedat->bp_allowed[row][l] != 0; l++) */
/*             { */
/*                bj = cedat->bp_allowed[row][l] - 1; */

/*                for (m = 0; m < allowed_bp; m++) */
/*                { */
/*                   nn_scores_get_allowed_basepair (m, &bip1, &bjm1, */
/*                                                   cedat->scores); */
                  
/*                   update_prob = */
/*                      seqmatrix_get_probability (bj, k, sm) */
/*                      * seqmatrix_get_probability (bip1, (col + 1), sm) */
/*                      * seqmatrix_get_probability (bjm1, (k - 1), sm); */
                  
/*                   tmp_neg += (nn_scores_get_G_stack (row, bj, bjm1, bip1, */
/*                                                      cedat->scores) */
/*                               * update_prob); */
/*                } */
/*             } */
/*             /\* mfprintf (stdout, "y: %.2f ", tmp2); *\/ */
/*          } */
/*          else         /\* k is i, col is j *\/ */
/*          { */
/*             /\* we only consider allowed base pairs *\/ */
/*             for (l = 0; cedat->bp_allowed[row][l] != 0; l++) */
/*             { */
/*                bi = cedat->bp_allowed[row][l] - 1; */
/*                for (m = 0; m < allowed_bp; m++) */
/*                { */
/*                   nn_scores_get_allowed_basepair (m, &bip1, &bjm1, */
/*                                                   cedat->scores); */
                  
/*                   update_prob = */
/*                      seqmatrix_get_probability (bi, k, sm) */
/*                      * seqmatrix_get_probability (bip1, (k + 1), sm) */
/*                      * seqmatrix_get_probability (bjm1, (col - 1), sm); */
                  
/*                   tmp_neg += (nn_scores_get_G_stack (bi, row, bjm1, bip1, */
/*                                                      cedat->scores) */
/*                               * update_prob); */

/*                } */
/*             } */
/*             /\* mfprintf (stdout, "x: %.2f ", tmp2); *\/ */
/*          } */
/*SNIP*/
         /* heterogenity term */
         if (col > k)
         {
            tmp_het += (seqmatrix_get_probability (row, k, sm)
                        * expf(cedat->het_rate * (col - (k+1))));
            het_count += expf ( (cedat->het_rate* (col - (k+1))));
         }
         else
         {
            tmp_het += (seqmatrix_get_probability (row, k, sm)
                        * expf(cedat->het_rate * (k - (col+1))));
            het_count += (expf (cedat->het_rate * (k - (col+1))));
         }         
      }
   }

   /*mfprintf (stderr, "c: %lu r: %lu Eneg: %.2f cell: %.2f\n",
     col, row, tmp_neg, cell);*/

   /*mprintf ("c: %lu r: %lu Eneg: %f\n", col, row, tmp_neg);*/

   tmp_neg = (tmp_neg / cols) * (-1.25f);
   tmp_het = (tmp_het ) * (3.0f);
   cell += tmp_neg;
   cell += tmp_het;

   return cell;
}


/* Tests nn:
--steps 100 -t10
(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....
ACGUGGCAAGGGCAACAGAUAGCCCAGGGGCACAGAGAGCCCCAGAUAGGGGCAUAGAGAGCCCCGCCACGUAGAG

Entropy dropout: 0.036480
Entropy dropout: 0.018240
Entropy dropout: -0.000000

real    0m49.946s
user    0m49.859s
sys     0m0.056s
================================================================================

--steps 100 -t10
(((((....)))))..(...).((()))
CAGGCAGAAGCCUGAAGACACAGGGCCC

Entropy dropout: 0.049463
Entropy dropout: -0.000000

real    0m2.254s
user    0m2.236s
sys     0m0.000s
================================================================================

--steps 100 -t10
...(((...)))(((())))...(((...(...))))
ACAGACACAGUCGGGGCCCCAGAGGCAGAGAGACGCC

Entropy dropout: 0.037461
Entropy dropout: -0.000000

real    0m6.744s
user    0m6.728s
sys     0m0.008s
================================================================================

--steps 100 -t10
(((((...(((...)))(((())))...(((...(...))))(((((....)))))..(...).((()))...)))))
GUGGCACAGACAUAGUCGGGGCCCCAUAGACAUAGAAACGUCGGGGCAAGAGCCCCAAGAGACAGGGCCCAGAGCCAC

Entropy dropout: 0.035544
Entropy dropout: 0.017773
Entropy dropout: -0.000000

real    0m39.088s
user    0m38.926s
sys     0m0.068s
================================================================================

*/

/* Test nussi
CGACGUUAAGUCGACAACAUACGACACGACGAUCACAACGUCGAACAAGCACGAUCACAACGUGCAACGUCGAAAC
*/
