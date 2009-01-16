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
 *  ToDo:
 *      - change term cols into n_sites -> log. unlinked from matrix geometry
 *      - change rows to n_states
 *      - change bp_allowed to char?
 *      - heterogenity: use scalar products of whole sites -> sort of
 *        heterogenity on multi sites, not single (gundolf)
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
      float het_scale;    /* scaling factor for het term */
      float neg_scale;    /* scaling factor for negative design */
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
   Scmf_Rna_Opt_data* cedat = XOBJ_MALLOC(sizeof (*cedat), file, line);

   if (cedat != NULL)
   {
      cedat->scores     = NULL;
      cedat->sigma      = NULL;
      cedat->rna        = NULL;
      cedat->bp_allowed = NULL;
      cedat->het_rate   = 0.0f;
      cedat->en_neg     = NULL;
      cedat->en_neg2    = NULL;
      cedat->en_neg_35  = NULL;
      cedat->het_scale  = 1.0f;
      cedat->neg_scale  = 1.0f;
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
   int error = 0;

   Scmf_Rna_Opt_data* this = scmf_rna_opt_data_new (file, line);

   if (this != NULL)
   {
      this->het_rate = het_rate;

      /* mfprintf (stderr, "Het-Rate: %f %f\n", het_rate, (logf (0.01) / logf (expf(1)))/het_rate); */

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
                                            sizeof (**(this->en_neg2)));
      this->en_neg_35 = (float**) XMALLOC_2D (alphabet_size (this->sigma),
                                           alphabet_size (this->sigma),
                                              sizeof (**(this->en_neg_35)));

      error = RNA_INIT_PAIRLIST_VIENNA(structure, seqlen, this->rna);
   }

   if (error)
   {
      scmf_rna_opt_data_delete (this);
      this = NULL;
   }

   return this;
}

int
scmf_rna_opt_data_secstruct_init (Scmf_Rna_Opt_data* this)
{
   return RNA_SECSTRUCT_INIT (this->rna);
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
   unsigned long /*i,*/ j, k, l, m, cols, alpha, allowed_bp;
   char bj, bip1, bjm1;
   float prob;
   Scmf_Rna_Opt_data* this;

   assert (sm);
   assert (data);

   this = (Scmf_Rna_Opt_data*) data;

   cols = seqmatrix_get_width (sm);
   alpha = alphabet_size (this->sigma);
   allowed_bp = nn_scores_no_allowed_basepairs (this->scores);
   /*i = 0;*/

   memset(this->en_neg2[0], 0, (alpha * alpha) * sizeof (**(this->en_neg2)));
   memset(this->en_neg_35[0], 0, (alpha * alpha) * sizeof(**(this->en_neg_35)));
   memset(this->en_neg, 0, alpha * sizeof (*(this->en_neg)));

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
      /* update negative energy for next col */
      /* subtract col, col+1 contribution (5' -> 3') */
      for (j = 0; this->bp_allowed[row][j] != 0; j++)
      {
         bj = (char)(this->bp_allowed[row][j] - 1);
         
         for (k = 0; k < allowed_bp; k++)
         {
            nn_scores_get_allowed_basepair (k, &bip1, &bjm1,
                                            this->scores);
            
            prob = seqmatrix_get_probability (bj, (col + 1), sm)
               * seqmatrix_get_probability (bjm1, col, sm);
            
            this->en_neg2[row][(int)bip1] -= (nn_scores_get_G_stack ((char)row,
                                                                     bj,
                                                                     bjm1, bip1,
                                                                     this->scores)
                                              * prob);
         }
      }

      /* update negative energy for next col */
      /* add col+1, col contribution (5' <- 3') */
      for (j = 0; this->bp_allowed[row][j] != 0; j++)
      {
         bi = (char)(this->bp_allowed[row][j] - 1);
         
         for (k = 0; k < allowed_bp; k++)
         {
            nn_scores_get_allowed_basepair (k, &bip1, &bjm1,
                                            this->scores);
            
            prob = seqmatrix_get_probability (bi, col, sm)
               * seqmatrix_get_probability (bip1, (col + 1), sm);
            
            this->en_neg_35[row][(int)bjm1] += (nn_scores_get_G_stack (bi,
                                                                      (char)row,
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
scmf_rna_opt_data_set_scales (float neg, float het, Scmf_Rna_Opt_data* cedat)
{
   assert (cedat);

   cedat->neg_scale = neg;
   cedat->het_scale = het;
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

   rna_set_sequence_base (alphabet_no_2_base((char)row, cont->sigma),
                          col, cont->rna);

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
scmf_rna_opt_calc_simplenn (const unsigned long row,
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
            bj = (char)(cedat->bp_allowed[row][l] - 1);
            
            /* pair each pair with all allowed pairs */
            for (k = 0; k < allowed_bp; k++)
            {
               nn_scores_get_allowed_basepair (k, &bip1, &bjm1, cedat->scores);
   
               G_stack_score = nn_scores_get_G_stack ((char)row, bj, bjm1, bip1,
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
            bj = (char)(cedat->bp_allowed[row][l] - 1);            
            
            for (k = 0; k < alpha_size; k++)
            {
               for (m = 0; m < alpha_size; m++)
               {
                  G_stack_score = nn_scores_get_G_mm_stack ((char)row, bj,
                                                            (char)m, (char)k,
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
            bi = (char)(cedat->bp_allowed[row][l] - 1); 
            
            /* pair each pair with all allowed pairs */
            for (k = 0; k < allowed_bp; k++)
            {
               nn_scores_get_allowed_basepair (k, &bip1, &bjm1,
                                               cedat->scores);
               
               G_stack_score = nn_scores_get_G_stack (bi, (char)row, bjm1,bip1,
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
            bi = (char)(cedat->bp_allowed[row][l] - 1); 
            
            /* pair with each possible pair */
            for (k = 0; k < alpha_size; k++)
            {
               for (m = 0; m < alpha_size; m++)
               {
                  G_stack_score = nn_scores_get_G_mm_stack (bi, (char)row,
                                                            (char)m, (char)k,
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
         bj = (char)(cedat->bp_allowed[row][l] - 1);
         
         for (k = 0; k < allowed_bp; k++)
         {
            nn_scores_get_allowed_basepair (k, &bip1, &bjm1,
                                            cedat->scores);
            
            G_stack_score = nn_scores_get_G_stack ((char)row, bj, bjm1, bip1,
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
         bi = (char)(cedat->bp_allowed[row][l] - 1);
         
         for (k = 0; k < allowed_bp; k++)
         {
            nn_scores_get_allowed_basepair (k, &bip1, &bjm1,
                                            cedat->scores);
            
            G_stack_score = nn_scores_get_G_stack (bi, (char) row, bjm1,bip1,
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


/* RNA design using the full NN model.
*/
static void
scmf_rna_opt_calc_hairpin (const unsigned long row,
                           const unsigned long hairpin,
                           SeqMatrix* sm,
                           Scmf_Rna_Opt_data* this)
{
   unsigned long start, end, size;
   char bpp, bi, bj; /* base pair partner */
   unsigned long k, l, m;
   unsigned long alpha_size;
   unsigned long n_tetra_loops;
   float cell5p, cell3p;
   float update_prob5p, update_prob3p;
   unsigned long allowed_bp = nn_scores_no_allowed_basepairs (this->scores);
   const char* t_loop;

   alpha_size = alphabet_size (this->sigma);

   rna_secstruct_get_geometry_hairpin(&start, &end, &size, hairpin, this->rna);
/*  mfprintf (stdout, "process hairpin loop %lu, %c %.2f: %lu - %lu | %lu\n", */
/*              hairpin, alphabet_no_2_base(row, this->sigma), */
/*              t, start, end, size); */

   /* closing base pair base */
   cell5p = 0.0f;
   cell3p = 0.0f;
   /* iterate all possible base pairs with the current base! */
   for (k = 0; this->bp_allowed[row][k] != 0; k++)
   {
      bpp = (char) (this->bp_allowed[row][k] - 1);
      
      /* for all possible unpaired bases */
      for (l = 0; l < alpha_size; l++)
      {
         update_prob5p =   seqmatrix_get_probability(bpp, end,     sm)
                         * seqmatrix_get_probability(l,   end - 1, sm);

         update_prob3p =   seqmatrix_get_probability(bpp, start,     sm)
                         * seqmatrix_get_probability(l,   start + 1, sm);

         for (m = 0; m < alpha_size; m++)
         {
            cell5p += (update_prob5p *
                       seqmatrix_get_probability(m, start + 1, sm)
               * nn_scores_get_G_hairpin_mismatch (row, bpp, m, l, size,
                                                   this->scores));

            cell3p += (update_prob3p *
                       seqmatrix_get_probability(m, end - 1, sm)
               * nn_scores_get_G_hairpin_mismatch (bpp, row, l, m, size,
                                                   this->scores));
         }
      }
   }
   
   /* store eeff */
   /* 4 bases make one mismatch/ closing bp */
   cell5p = cell5p / 4; /* SB 08-12-12 */
   cell3p = cell3p / 4; /* SB 08-12-12 */
   seqmatrix_add_2_eeff (cell5p, row, start, sm);
   seqmatrix_add_2_eeff (cell3p, row,   end, sm);

   /* process opening "base pair" (i+1, j-1) */
   cell5p = 0.0f;
   cell3p = 0.0f;
   /* for all possible closing base pairs */
   for (k = 0; k < allowed_bp; k++)
   {
      /* get bases and base probability*/
      nn_scores_get_allowed_basepair (k, &bi, &bj, this->scores);
      update_prob5p =  seqmatrix_get_probability(bi, start,     sm)
                     * seqmatrix_get_probability(bj,   end, sm);

      /* for all possible unpaired bases */
      for (l = 0; l < alpha_size; l++)
      {
         cell5p += (update_prob5p * seqmatrix_get_probability(l, end - 1, sm))
            * nn_scores_get_G_hairpin_mismatch (bi, bj, row, l, size,
                                                      this->scores);
         
         cell3p += (update_prob5p * seqmatrix_get_probability(l, start + 1, sm))
            * nn_scores_get_G_hairpin_mismatch (bi, bj, l, row, size,
                                                      this->scores);
      }
   }
   /* 4 bases make one mismatch/ closing bp */
   cell5p = cell5p / 4; /* SB 08-12-12 */
   cell3p = cell3p / 4; /* SB 08-12-12 */
   seqmatrix_add_2_eeff (cell5p, row, start + 1, sm);
   seqmatrix_add_2_eeff (cell3p, row, end - 1  , sm);

   /* modeling tetraloops */
   if (size == nn_scores_get_size_tetra_loop(this->scores))
   {
      size = nn_scores_get_size_tetra_loop_full(this->scores);
      /* loop over all tetraloops */
      n_tetra_loops = nn_scores_get_no_of_tetra_loops (this->scores);
      for (k = 0; k < n_tetra_loops; k++)
      {
         /* fetch current loop */
         t_loop = nn_scores_get_tetra_loop (k, this->scores);
         update_prob5p = 1.0;

         /* calc eeff on base of ALL involved bases 
            this is the loop size and its closing bp */
         for (l = 0; l < size; l++)
         {
            update_prob5p *= seqmatrix_get_probability(t_loop[l],
                                                       start + l,
                                                       sm);
         }
         /* tetra loop score is for whole loop, so we spread the value */
         cell5p = update_prob5p *
            nn_scores_get_G_tetra_loop (t_loop, 0,this->scores) / size;

         /* add eeff to matching cells */
         for (l = 0; l < size; l++)
         {
            if ((unsigned) t_loop[l] == row)
            {
               cell3p = cell5p / seqmatrix_get_probability(t_loop[l],
                                                           start + l,
                                                           sm);
               seqmatrix_add_2_eeff (cell3p, row, start + l, sm);
            }
         }         
      }
   }
}

static void
scmf_rna_opt_calc_ext_loop (const unsigned long row,
                            SeqMatrix* sm,
                            Scmf_Rna_Opt_data* this)
{
   unsigned long n;                   /* general count */
   unsigned long k, l, m;
   unsigned long p5pos, p3pos, fbpos; /* pos of 5', 3' and free base */
   unsigned long alpha_size = alphabet_size (this->sigma);
   unsigned long allowed_bp = nn_scores_no_allowed_basepairs (this->scores);
   char bpp, bi, bj;
   float cell5p, cell3p;
   float p5p, p3p;                    /* probability */

   /* penalty for non gc basepair initiating a stem */
   n = rna_secstruct_get_noof_stems_extloop (this->rna);

   for (k = 0; k < n; k++)
   {
      cell5p = 0.0f;
      cell3p = 0.0f;
      rna_secstruct_get_i_stem_extloop (&p5pos, &p3pos, k, this->rna);

      /* for all possible closing bp */
      for (l = 0; this->bp_allowed[row][l] != 0; l++)
      {
         bpp = (char) (this->bp_allowed[row][l] - 1);

         cell5p += (seqmatrix_get_probability(bpp, p3pos, sm)
                     * nn_scores_get_G_non_gc_penalty_for_bp (row, bpp,
                                                             this->scores));

         cell3p += (seqmatrix_get_probability(bpp, p5pos, sm)
                    * nn_scores_get_G_non_gc_penalty_for_bp (bpp, row,
                                                             this->scores));
         /*mfprintf (stderr, "NGCP: %f\n", 
                   nn_scores_get_G_non_gc_penalty_for_bp(bpp, row,
                   this->scores));*/
      }

      /* usually the non_gc_penalty counts per pair, we
         spread it over the involved bases */
      cell5p = cell5p / 2; /* SB 08-12-11 */
      cell3p = cell3p / 2; /* SB 08-12-11 */
      seqmatrix_add_2_eeff (cell5p, row, p5pos, sm);
      seqmatrix_add_2_eeff (cell3p, row, p3pos, sm);
   }

   /* 5' dangle */
   n = rna_secstruct_get_noof_5pdangles_extloop (this->rna);
   for (k = 0; k < n; k++)
   {
      cell5p = 0.0f;
      cell3p = 0.0f;

      rna_secstruct_get_i_5pdangle_extloop (&p5pos, &p3pos, &fbpos, k,
                                            this->rna);

      /* design bp: for all bp allowed with 'row' and all free bases */
      for (l = 0; this->bp_allowed[row][l] != 0; l++)
      {
         bpp = (char) (this->bp_allowed[row][l] - 1);

         p5p = seqmatrix_get_probability(bpp, p3pos, sm);
         p3p = seqmatrix_get_probability(bpp, p5pos, sm);

         for (m = 0; m < alpha_size; m++)
         {
            cell5p += (p5p * seqmatrix_get_probability(m, fbpos, sm)
                       * nn_scores_get_G_dangle5 (row, bpp, m, this->scores));
            cell3p += (p3p * seqmatrix_get_probability(m, fbpos, sm)
                       * nn_scores_get_G_dangle5 (bpp, row, m, this->scores));
         }
      }
      /* 3 bases are involved in 5p dangle */
      cell5p = cell5p / 3; /* SB 08-12-11 */
      cell3p = cell3p / 3; /* SB 08-12-11 */
      seqmatrix_add_2_eeff (cell5p, row, p5pos, sm);
      seqmatrix_add_2_eeff (cell3p, row, p3pos, sm);

      /* design free base: for all allowed bp with free base 'row' */
      cell5p = 0.0f;
      for (l = 0; l < allowed_bp; l++)
      {
         nn_scores_get_allowed_basepair (l, &bi, &bj, this->scores);

         cell5p += (seqmatrix_get_probability(bi, p5pos, sm)
                    * seqmatrix_get_probability(bj, p3pos, sm)
                    * nn_scores_get_G_dangle5 (bi, bj, row, this->scores));
      }
      /* 3 bases are involved in 5p dangle */
      cell5p = cell5p / 3; /* SB 08-12-11 */
      seqmatrix_add_2_eeff (cell5p, row, fbpos, sm);
   }

   /* 3' dangle */
   n = rna_secstruct_get_noof_3pdangles_extloop (this->rna);
   for (k = 0; k < n; k++)
   {
      cell5p = 0.0f;
      cell3p = 0.0f;

      rna_secstruct_get_i_3pdangle_extloop (&p5pos, &p3pos, &fbpos, k,
                                            this->rna);

      /* design bp: for all bp allowed with 'row' and all free bases */
      for (l = 0; this->bp_allowed[row][l] != 0; l++)
      {
         bpp = (char) (this->bp_allowed[row][l] - 1);

         p5p = seqmatrix_get_probability(bpp, p3pos, sm);
         p3p = seqmatrix_get_probability(bpp, p5pos, sm);

         for (m = 0; m < alpha_size; m++)
         {
            cell5p += (p5p * seqmatrix_get_probability(m, fbpos, sm)
                       * nn_scores_get_G_dangle3 (row, bpp, m, this->scores));
            cell3p += (p3p * seqmatrix_get_probability(m, fbpos, sm)
                       * nn_scores_get_G_dangle3 (bpp, row, m, this->scores));
         }
      }
      /* 3 bases are involved in 3p dangle */
      cell5p = cell5p / 3; /* SB 08-12-11 */
      cell3p = cell3p / 3; /* SB 08-12-11 */
      seqmatrix_add_2_eeff (cell5p, row, p5pos, sm);
      seqmatrix_add_2_eeff (cell3p, row, p3pos, sm);

      /* design free base: for all allowed bp with free base 'row' */
      cell3p = 0.0f;
      for (l = 0; l < allowed_bp; l++)
      {
         nn_scores_get_allowed_basepair (l, &bi, &bj, this->scores);

         cell3p += (seqmatrix_get_probability(bi, p5pos, sm)
                    * seqmatrix_get_probability(bj, p3pos, sm)
                    * nn_scores_get_G_dangle3 (bi, bj, row, this->scores));
      }
      /* 3 bases are involved in 3p dangle */
      cell3p = cell3p / 3; /* SB 08-12-11 */
      seqmatrix_add_2_eeff (cell3p, row, fbpos, sm);
   }
}

static void
scmf_rna_opt_calc_multi_loop (const unsigned long row,
                              const unsigned long loop,
                              SeqMatrix* sm,
                              Scmf_Rna_Opt_data* this)
{
   unsigned long n;                   /* count */
   unsigned long k, l, m;             /* iterator */
   float cell5p, cell3p;              /* tmp. stor. */
   unsigned long p5pos, p3pos, fbpos; /* pos of 5', 3' and free base */
   char bpp, bi, bj;
   float p5p, p3p;                    /* probabilities */
   unsigned long alpha_size = alphabet_size (this->sigma);
   unsigned long allowed_bp = nn_scores_no_allowed_basepairs (this->scores);

   /* penalty for non gc basepair initiating a stem */
   n = rna_secstruct_get_i_noof_stems_multiloop (loop, this->rna);

   for (k = 0; k < n; k++)
   {
      cell5p = 0.0f;
      cell3p = 0.0f;

      rna_secstruct_get_i_stem_multiloop (&p5pos, &p3pos, k, loop, this->rna);

      /* for all possible closing bp */
      for (l = 0; this->bp_allowed[row][l] != 0; l++)
      {
         bpp = (char) (this->bp_allowed[row][l] - 1);

         cell5p += (seqmatrix_get_probability(bpp, p3pos, sm)
                     * nn_scores_get_G_non_gc_penalty_for_bp (row, bpp,
                                                             this->scores));

         cell3p += (seqmatrix_get_probability(bpp, p5pos, sm)
                    * nn_scores_get_G_non_gc_penalty_for_bp (bpp, row,
                                                             this->scores));
      }
      /* usually the non_gc_penalty counts per pair, we
         spread it over the involved bases */
      cell5p = cell5p / 2; /* SB 08-12-11 */
      cell3p = cell3p / 2; /* SB 08-12-11 */
      seqmatrix_add_2_eeff (cell5p, row, p5pos, sm);
      seqmatrix_add_2_eeff (cell3p, row, p3pos, sm);
   }

   /* 5' dangle */
   n = rna_secstruct_get_i_noof_5pdangles_multiloop (loop, this->rna);
   for (k = 0; k < n; k++)
   {
      cell5p = 0.0f;
      cell3p = 0.0f;

      rna_secstruct_get_i_5pdangle_multiloop (&p5pos, &p3pos, &fbpos, k, loop,
                                              this->rna);

      /* design bp: for all bp allowed with 'row' and all free bases */
      for (l = 0; this->bp_allowed[row][l] != 0; l++)
      {
         bpp = (char) (this->bp_allowed[row][l] - 1);

         p5p = seqmatrix_get_probability(bpp, p3pos, sm);
         p3p = seqmatrix_get_probability(bpp, p5pos, sm);

         for (m = 0; m < alpha_size; m++)
         {
            cell5p += (p5p * seqmatrix_get_probability(m, fbpos, sm)
                       * nn_scores_get_G_dangle5 (row, bpp, m, this->scores));
            cell3p += (p3p * seqmatrix_get_probability(m, fbpos, sm)
                       * nn_scores_get_G_dangle5 (bpp, row, m, this->scores));
         }
      }
      /* all 3 bases get equal contribution */
      cell5p = cell5p / 3; /* SB 08-12-11 */
      cell3p = cell3p / 3; /* SB 08-12-11 */
      seqmatrix_add_2_eeff (cell5p, row, p5pos, sm);
      seqmatrix_add_2_eeff (cell3p, row, p3pos, sm);

      /* design free base: for all allowed bp with free base 'row' */
      cell5p = 0.0f;
      for (l = 0; l < allowed_bp; l++)
      {
         nn_scores_get_allowed_basepair (l, &bi, &bj, this->scores);

         cell5p += (seqmatrix_get_probability(bi, p5pos, sm)
                    * seqmatrix_get_probability(bj, p3pos, sm)
                    * nn_scores_get_G_dangle5 (bi, bj, row, this->scores));
      }
      cell5p = cell5p / 3; /* SB 08-12-11 */
      seqmatrix_add_2_eeff (cell5p, row, fbpos, sm);
   }

   /* 3' dangle */
   n = rna_secstruct_get_i_noof_3pdangles_multiloop (loop, this->rna);
   for (k = 0; k < n; k++)
   {
      cell5p = 0.0f;
      cell3p = 0.0f;

      rna_secstruct_get_i_3pdangle_multiloop (&p5pos, &p3pos, &fbpos, k, loop,
                                              this->rna);

      /* design bp: for all bp allowed with 'row' and all free bases */
      for (l = 0; this->bp_allowed[row][l] != 0; l++)
      {
         bpp = (char) (this->bp_allowed[row][l] - 1);

         p5p = seqmatrix_get_probability(bpp, p3pos, sm);
         p3p = seqmatrix_get_probability(bpp, p5pos, sm);

         for (m = 0; m < alpha_size; m++)
         {
            cell5p += (p5p * seqmatrix_get_probability(m, fbpos, sm)
                       * nn_scores_get_G_dangle3 (row, bpp, m, this->scores));
            cell3p += (p3p * seqmatrix_get_probability(m, fbpos, sm)
                       * nn_scores_get_G_dangle3 (bpp, row, m, this->scores));
         }
      }
      /* all 3 bases get equal contribution */
      cell5p = cell5p / 3; /* SB 08-12-11 */
      cell3p = cell3p / 3; /* SB 08-12-11 */
      seqmatrix_add_2_eeff (cell5p, row, p5pos, sm);
      seqmatrix_add_2_eeff (cell3p, row, p3pos, sm);

      /* design free base: for all allowed bp with free base 'row' */
      cell3p = 0.0f;
      for (l = 0; l < allowed_bp; l++)
      {
         nn_scores_get_allowed_basepair (l, &bi, &bj, this->scores);

         cell3p += (seqmatrix_get_probability(bi, p5pos, sm)
                    * seqmatrix_get_probability(bj, p3pos, sm)
                    * nn_scores_get_G_dangle3 (bi, bj, row, this->scores));

         /*mprintf ("%lu %lu %c%c - %c: %11f %4i %f\n", row, k,
                  alphabet_no_2_base (bi, this->sigma),
                  alphabet_no_2_base (bj, this->sigma),
                  alphabet_no_2_base (row, this->sigma),
                  cell3p,
                  nn_scores_get_G_dangle3 (bi, bj, row, this->scores),
                  (seqmatrix_get_probability(bi, p5pos, sm)
                  * seqmatrix_get_probability(bj, p3pos, sm)));*/
      }
      cell3p = cell3p / 3; /* SB 08-12-11 */
      seqmatrix_add_2_eeff (cell3p, row, fbpos, sm);
   }
}

static void
scmf_rna_opt_calc_bulge (const unsigned long row,
                         const unsigned long loop,
                         SeqMatrix* sm,
                         Scmf_Rna_Opt_data* this)
{
   unsigned long i1pos;
   unsigned long j1pos;
   unsigned long i2pos;
   unsigned long j2pos;
   unsigned long size;
   unsigned long k, l;           /* iterator */
   unsigned long allowed_bp = nn_scores_no_allowed_basepairs (this->scores);
   char bpp, bi, bj;
   float cell_i1, cell_j1, cell_i2, cell_j2;
   float p;                     /* probability */

   /* fetch loop geometry */
   rna_secstruct_get_geometry_bulge (&i1pos, &j1pos, &i2pos, &j2pos, &size,
                                     loop,
                                     this->rna);

   /* loop over allowed bp */
   cell_i1 = 0.0f;
   cell_j1 = 0.0f;
   cell_i2 = 0.0f;
   cell_j2 = 0.0f;
   for (k = 0; this->bp_allowed[row][k] != 0; k++)
   {
      bpp = (char) (this->bp_allowed[row][k] - 1);

      /* combine with all possible bp */
      for (l = 0; l < allowed_bp; l++)
      {
         nn_scores_get_allowed_basepair (l, &bi, &bj, this->scores);

         /* design first pair */
         p = seqmatrix_get_probability(bi, i2pos, sm) 
            * seqmatrix_get_probability(bj, j2pos, sm);

         cell_i1 += (seqmatrix_get_probability(bpp, j1pos, sm) * p
                     * nn_scores_get_G_bulge_stack (row, bpp, bj, bi,
                                                    size, this->scores));

         cell_j1 += (seqmatrix_get_probability(bpp, i1pos, sm) * p
                     * nn_scores_get_G_bulge_stack (bpp, row, bj, bi,
                                                    size, this->scores));

         /* design 2nd pair */
         p = seqmatrix_get_probability(bi, i1pos, sm)
            * seqmatrix_get_probability(bj, j1pos, sm);

         cell_i2 += (seqmatrix_get_probability(bpp, j2pos, sm) * p
                     * nn_scores_get_G_bulge_stack (bi, bj, bpp, row,
                                                    size, this->scores));

         cell_j2 += (seqmatrix_get_probability(bpp, i2pos, sm) * p
                     * nn_scores_get_G_bulge_stack (bi, bj, row, bpp,
                                                    size, this->scores));
      }
   }

   /* we have 2 base pairs, so we spread on 4 */
   cell_i1 = cell_i1 / 4; /* SB 08-12-12 */
   cell_j1 = cell_j1 / 4; /* SB 08-12-12 */
   cell_i2 = cell_i2 / 4; /* SB 08-12-12 */
   cell_j2 = cell_j2 / 4; /* SB 08-12-12 */
   seqmatrix_add_2_eeff (cell_i1, row, i1pos, sm);
   seqmatrix_add_2_eeff (cell_j1, row, j1pos, sm);
   seqmatrix_add_2_eeff (cell_i2, row, i2pos, sm);
   seqmatrix_add_2_eeff (cell_j2, row, j2pos, sm);
}

static void
scmf_rna_opt_calc_stack (const unsigned long row,
                         const unsigned long stack,
                         SeqMatrix* sm,
                         Scmf_Rna_Opt_data* this)
{
   unsigned long i, j;          /* base i and j of a pair */
   unsigned long k, l;          /* iterator */
   char bpp, bi, bj;
   unsigned long allowed_bp = nn_scores_no_allowed_basepairs (this->scores);
   float cell_i, cell_j, cell_ip1, cell_jm1;
   float p;                     /* probability */

   /* get base pair */
   rna_secstruct_get_i_geometry_stack (&i, &j, stack, this->rna);

   /* for all allowed pairs */
   cell_i   = 0.0f;
   cell_j   = 0.0f;
   cell_ip1 = 0.0f;
   cell_jm1 = 0.0f;
   for (k = 0; this->bp_allowed[row][k] != 0; k++)
   {
      bpp = (char) (this->bp_allowed[row][k] - 1);

      /* for all possible pairs */
      for (l = 0; l < allowed_bp; l++)
      {
         nn_scores_get_allowed_basepair (l, &bi, &bj, this->scores);
      
         p = seqmatrix_get_probability( bi, i + 1, sm)
            * seqmatrix_get_probability(bj, j - 1, sm);

         cell_i += (seqmatrix_get_probability(bpp, j, sm) * p
                  * nn_scores_get_G_stack (row, bpp, bj, bi, this->scores));

         cell_j += (seqmatrix_get_probability(bpp, i, sm) * p
                    * nn_scores_get_G_stack (bpp, row, bj, bi, this->scores));

         p = seqmatrix_get_probability( bi, i, sm)
            * seqmatrix_get_probability(bj, j, sm);

         cell_ip1 += (seqmatrix_get_probability(bpp, j - 1, sm) * p
                      * nn_scores_get_G_stack (bi, bj, bpp, row, this->scores));

         cell_jm1 += (seqmatrix_get_probability(bpp, i + 1, sm) * p
                      * nn_scores_get_G_stack (bi, bj, row, bpp, this->scores));
      }
   }
   /* we distribute the energy values on 4 bases */
   cell_i   = cell_i   / 4; /* SB 08-12-12 */
   cell_j   = cell_j   / 4; /* SB 08-12-12 */
   cell_ip1 = cell_ip1 / 4; /* SB 08-12-12 */
   cell_jm1 = cell_jm1 / 4; /* SB 08-12-12 */
   seqmatrix_add_2_eeff (cell_i,   row, i,     sm);
   seqmatrix_add_2_eeff (cell_j,   row, j,     sm);
   seqmatrix_add_2_eeff (cell_ip1, row, i + 1, sm);
   seqmatrix_add_2_eeff (cell_jm1, row, j - 1, sm);
}

static __inline__ void
scmf_rna_opt_calc_internal (const unsigned long row,
                            unsigned long allowed_bp,
                            unsigned long alpha_size,
                            unsigned long pi1,
                            unsigned long pj1,
                            unsigned long pi2,
                            unsigned long pj2,
                            SeqMatrix* sm,
                            Scmf_Rna_Opt_data* this)
{
   unsigned long k, l, m;
   char bpp, bi, bj;
   float p;
   float cell_i1, cell_j1, cell_j2, cell_i2;

   /*mfprintf (stderr,"i1: %lu, j1: %lu, i2: %lu, j2: %lu\n", i1, j1, i2, j2);*/
   
   /* design paired bases */
   cell_i1 = 0.0f;
   cell_j1 = 0.0f;
   cell_j2 = 0.0f;
   cell_i2 = 0.0f;
   /* for all allowed pairs */
   for (k = 0; this->bp_allowed[row][k] != 0; k++)
   {
      bpp = (char) (this->bp_allowed[row][k] - 1);

      /* for all bases */
      for (l = 0; l < alpha_size; l++)
      {
         /* for all bases */
         for (m = 0; m < alpha_size; m++)
         {
            /* i1 */
            p = seqmatrix_get_probability(bpp, pj1, sm)
               * seqmatrix_get_probability(l, pi1 + 1, sm)
               * seqmatrix_get_probability(m, pj1 - 1, sm);

            cell_i1 += p * nn_scores_get_G_mismatch_interior (row,
                                                              bpp,
                                                              l,
                                                              m,
                                                              this->scores);
            
            /* j1 */
            p = seqmatrix_get_probability(bpp, pi1, sm)
               * seqmatrix_get_probability(l, pi1 + 1, sm)
               * seqmatrix_get_probability(m, pj1 - 1, sm);

            cell_j1 += p * nn_scores_get_G_mismatch_interior (bpp,
                                                              row,
                                                              l,
                                                              m,
                                                              this->scores);

            /* j2 */
            p = seqmatrix_get_probability(bpp, pi2, sm)
               * seqmatrix_get_probability(l, pj2 + 1, sm)
               * seqmatrix_get_probability(m, pi2 - 1, sm);

            cell_j2 += p * nn_scores_get_G_mismatch_interior (row,
                                                              bpp,
                                                              l,
                                                              m,
                                                              this->scores);

            /* i2 */
            p = seqmatrix_get_probability(bpp, pj2, sm)
               * seqmatrix_get_probability(l, pj2 + 1, sm)
               * seqmatrix_get_probability(m, pi2 - 1, sm);

            cell_i2 += p * nn_scores_get_G_mismatch_interior (bpp,
                                                              row,
                                                              l,
                                                              m,
                                                              this->scores);
         }
      }
   }
   /* for each site of the loop 4 bases are involved */
   cell_i1 = cell_i1 / 4;       /* SB 08-12-12 */
   cell_j1 = cell_j1 / 4;       /* SB 08-12-12 */
   cell_i2 = cell_i2 / 4;       /* SB 08-12-12 */
   cell_j2 = cell_j2 / 4;       /* SB 08-12-12 */
   seqmatrix_add_2_eeff (cell_i1, row, pi1, sm);
   seqmatrix_add_2_eeff (cell_j1, row, pj1, sm);
   seqmatrix_add_2_eeff (cell_j2, row, pj2, sm);
   seqmatrix_add_2_eeff (cell_i2, row, pi2, sm);

   /* design unpaired bases */
   cell_i1 = 0.0f;
   cell_j1 = 0.0f;
   cell_j2 = 0.0f;
   cell_i2 = 0.0f;
   /* for all possible pairs */
   for (k = 0; k < allowed_bp; k++)
   {
      nn_scores_get_allowed_basepair (k, &bi, &bj, this->scores);

      /* for all bases */
      for (l = 0; l < alpha_size; l++)
      {
         /* i1 + 1 */
         p = seqmatrix_get_probability(bi, pi1, sm)
            * seqmatrix_get_probability(bj, pj1, sm)
            * seqmatrix_get_probability(l, pj1 - 1, sm);

         cell_i1 += p * nn_scores_get_G_mismatch_interior (bi,
                                                           bj,
                                                           row,
                                                           l,
                                                           this->scores);

         /* j1 - 1 */
         p = seqmatrix_get_probability(bi, pi1, sm)
            * seqmatrix_get_probability(bj, pj1, sm)
            * seqmatrix_get_probability(l, pi1 + 1, sm);

         cell_j1 += p * nn_scores_get_G_mismatch_interior (bi,
                                                           bj,
                                                           l,
                                                           row,
                                                           this->scores);

         /* j2 + 1 */
         p = seqmatrix_get_probability(bj, pj2, sm)
            * seqmatrix_get_probability(bi, pi2, sm)
            * seqmatrix_get_probability(l, pi2 - 1, sm);

         cell_j2 += p * nn_scores_get_G_mismatch_interior (bj,
                                                           bi,
                                                           row,
                                                           l,
                                                           this->scores);

         /* i2 - 1 */
         p = seqmatrix_get_probability(bj, pj2, sm)
            * seqmatrix_get_probability(bi, pi2, sm)
            * seqmatrix_get_probability(l, pj2 + 1, sm);
        
         cell_i2 += p * nn_scores_get_G_mismatch_interior (bj,
                                                           bi,
                                                           l,
                                                           row,
                                                           this->scores);
      }
   }
   /* for each site of the loop 4 bases are involved */
   cell_i1 = cell_i1 / 4;       /* SB 08-12-12 */
   cell_j1 = cell_j1 / 4;       /* SB 08-12-12 */
   cell_i2 = cell_i2 / 4;       /* SB 08-12-12 */
   cell_j2 = cell_j2 / 4;       /* SB 08-12-12 */
   seqmatrix_add_2_eeff (cell_i1, row, pi1 + 1, sm);
   seqmatrix_add_2_eeff (cell_j1, row, pj1 - 1, sm);
   seqmatrix_add_2_eeff (cell_j2, row, pj2 + 1, sm);
   seqmatrix_add_2_eeff (cell_i2, row, pi2 - 1, sm);
}

static __inline__ void
scmf_rna_opt_calc_int22 (const unsigned long row,
                         unsigned long allowed_bp,
                         unsigned long alpha_size,
                         unsigned long pi1,
                         unsigned long pj1,
                         unsigned long pi2,
                         unsigned long pj2,
                         SeqMatrix* sm,
                         Scmf_Rna_Opt_data* this)
{
   unsigned long k, l, m, n, o, p;
   char bpp, bi, bj, bi2, bj2;
   float pr, p_bi1pi2m, p_bj2p, p_bp1, p_bp2;
   float cell_i1, cell_j1, cell_i2, cell_j2;

   /*mfprintf (stderr,"i1: %lu, j1: %lu, i2: %lu, j2: %lu\n", i1, j1, i2, j2);*/

   /* design paired bases */
   cell_i1 = 0.0f;
   cell_j1 = 0.0f;
   cell_i2 = 0.0f;
   cell_j2 = 0.0f;
   /* for all allowed pairs */
   for (k = 0; this->bp_allowed[row][k] != 0; k++)
   {
      bpp = (char) (this->bp_allowed[row][k] - 1);

      /* for all possible pairs */
      for (l = 0; l < allowed_bp; l++)
      {
         nn_scores_get_allowed_basepair (l, &bi2, &bj2, this->scores);

         p_bp2 = seqmatrix_get_probability(bj2, pj2, sm)
            * seqmatrix_get_probability(bi2, pi2, sm);

         p_bp1 = seqmatrix_get_probability(bi2, pi1, sm)
            * seqmatrix_get_probability(bj2, pj1, sm);

         /* for all bases */
         for (m = 0; m < alpha_size; m++)
         {
            /* for all bases */
            for (n = 0; n < alpha_size; n++)
            {
               p_bi1pi2m = seqmatrix_get_probability(m, pi1 + 1, sm)
                  * seqmatrix_get_probability(n, pi2 - 1, sm);

               /* for all bases */
               for (o = 0; o < alpha_size; o++)
               {
                  p_bj2p = p_bi1pi2m * seqmatrix_get_probability(o, pj2 + 1, sm);

                  /* for all bases */
                  for (p = 0; p < alpha_size; p++)
                  {
                     /* i1 */
                     pr =p_bj2p * p_bp2
                        * seqmatrix_get_probability(p, pj1 - 1, sm)
                        * seqmatrix_get_probability(bpp, pj1, sm);

                     cell_i1 += pr *
                        nn_scores_get_G_internal_2x2_loop (row, bpp,
                                                           m, n,
                                                           bj2, bi2,
                                                           o, p,
                                                           this->scores);

                     /* j1 */
                     pr = p_bj2p * p_bp2
                        * seqmatrix_get_probability(p, pj1 - 1, sm)
                        * seqmatrix_get_probability(bpp, pi1, sm);

                     cell_j1 += pr *
                        nn_scores_get_G_internal_2x2_loop (bpp, row,
                                                           m, n,
                                                           bj2, bi2,
                                                           o, p,
                                                           this->scores);

                     /* i2 */
                     pr = p_bj2p * p_bp1
                        * seqmatrix_get_probability(p, pj1 - 1, sm)
                        * seqmatrix_get_probability(bpp, pj2, sm);

                     cell_i2 += pr *
                        nn_scores_get_G_internal_2x2_loop (bi2, bj2,
                                                           m, n,
                                                           bpp, row,
                                                           o, p,
                                                           this->scores);

                     /* j2 */
                     pr = p_bj2p *  p_bp1
                        * seqmatrix_get_probability(p, pj1 - 1, sm)
                        * seqmatrix_get_probability(bpp, pi2, sm);

                     cell_j2 += pr *
                        nn_scores_get_G_internal_2x2_loop (bi2, bj2,
                                                           m, n,
                                                           row, bpp,
                                                           o, p,
                                                           this->scores);

                  }
               }
            }
         }
      }
   }
   /* 8 bases are involved in a 2x2 loop */
   cell_i1 = cell_i1 / 8; /* SB 08-12-12 */
   cell_j1 = cell_j1 / 8; /* SB 08-12-12 */
   cell_i2 = cell_i2 / 8; /* SB 08-12-12 */
   cell_j2 = cell_j2 / 8; /* SB 08-12-12 */
   seqmatrix_add_2_eeff (cell_i1, row, pi1, sm);
   seqmatrix_add_2_eeff (cell_j1, row, pj1, sm);
   seqmatrix_add_2_eeff (cell_i2, row, pi2, sm);
   seqmatrix_add_2_eeff (cell_j2, row, pj2, sm);

   /* design unpaired bases */
   cell_i1 = 0.0f;
   cell_i2 = 0.0f;
   cell_j2 = 0.0f;
   cell_j1 = 0.0f;
   /* for all possible pairs */
   for (k = 0; k < allowed_bp; k++)
   {
      nn_scores_get_allowed_basepair (k, &bi, &bj, this->scores);

      p_bp1 = seqmatrix_get_probability(bi, pi1, sm)
         * seqmatrix_get_probability(bj, pj1, sm);

      /* for all possible pairs */
      for (l = 0; l < allowed_bp; l++)
      {
         nn_scores_get_allowed_basepair (l, &bi2, &bj2, this->scores);

         p_bp2 =  seqmatrix_get_probability(bj2, pj2, sm)
            * seqmatrix_get_probability(bi2, pi2, sm);

         /* for all bases */
         for (m = 0; m < alpha_size; m++)
         {
            /* for all bases */
            for (n = 0; n < alpha_size; n++)
            {
               /* for all bases */
               for (o = 0; o < alpha_size; o++)
               {
                  /* i1 + 1 */
                  pr = p_bp1 * p_bp2
                     * seqmatrix_get_probability(m, pi2 - 1, sm)
                     * seqmatrix_get_probability(n, pj2 + 1, sm)
                     * seqmatrix_get_probability(o, pj1 - 1, sm);

                  cell_i1 += pr *
                        nn_scores_get_G_internal_2x2_loop (bi, bj,
                                                           row, m,
                                                           bj2, bi2,
                                                           n, o,
                                                           this->scores);

                  /* i2 - 1 */
                  pr = p_bp1 * p_bp2
                     * seqmatrix_get_probability(m, pi1 + 1, sm)
                     * seqmatrix_get_probability(n, pj2 + 1, sm)
                     * seqmatrix_get_probability(o, pj1 - 1, sm);

                  cell_i2 += pr *
                        nn_scores_get_G_internal_2x2_loop (bi, bj,
                                                           m, row,
                                                           bj2, bi2,
                                                           n, o,
                                                           this->scores);

                  /* j2 + 1 */
                  pr = p_bp1 * p_bp2
                     * seqmatrix_get_probability(m, pi1 + 1, sm)
                     * seqmatrix_get_probability(n, pi2 - 1, sm)
                     * seqmatrix_get_probability(o, pj1 - 1, sm);

                  cell_j2 += pr *
                        nn_scores_get_G_internal_2x2_loop (bi, bj,
                                                           m, n,
                                                           bj2, bi2,
                                                           row, o,
                                                           this->scores);

                  /* j1 - 1 */
                  pr = p_bp1 * p_bp2
                     * seqmatrix_get_probability(m, pi1 + 1, sm)
                     * seqmatrix_get_probability(n, pi2 - 1, sm)
                     * seqmatrix_get_probability(o, pj2 + 1, sm);
                  
                  cell_j1 += pr *
                        nn_scores_get_G_internal_2x2_loop (bi, bj,
                                                           m, n,
                                                           bj2, bi2,
                                                           o, row,
                                                           this->scores);
               }
            }
         }
      }
   }
   /* 8 bases are involved in a 2x2 loop */
   cell_i1 = cell_i1 / 8; /* SB 08-12-12 */
   cell_j1 = cell_j1 / 8; /* SB 08-12-12 */
   cell_i2 = cell_i2 / 8; /* SB 08-12-12 */
   cell_j2 = cell_j2 / 8; /* SB 08-12-12 */
   seqmatrix_add_2_eeff (cell_i1, row, pi1 + 1, sm);
   seqmatrix_add_2_eeff (cell_i2, row, pi2 - 1, sm);
   seqmatrix_add_2_eeff (cell_j2, row, pj2 + 1, sm);
   seqmatrix_add_2_eeff (cell_j1, row, pj1 - 1, sm);
}

static void
scmf_rna_opt_calc_int12 (const unsigned long row,
                         unsigned long allowed_bp,
                         unsigned long alpha_size,
                         unsigned long pi1,
                         unsigned long pj1,
                         unsigned long pi2,
                         unsigned long pj2,
                         SeqMatrix* sm,
                         Scmf_Rna_Opt_data* this)
{
   unsigned long k, l, m, n, o;
   char bpp, bi, bj, bi2, bj2;
   float p, p_bp2, p_bp1, p_bb;
   float cell_i1, cell_j1, cell_i2, cell_j2;

   /* for all allowed pairs */
   cell_i1 = 0.0f;
   cell_j1 = 0.0f;
   cell_i2 = 0.0f;
   cell_j2 = 0.0f;
   for (k = 0; this->bp_allowed[row][k] != 0; k++)
   {
      bpp = (char) (this->bp_allowed[row][k] - 1);

      /* for all possible pairs */
      for (l = 0; l < allowed_bp; l++)
      {
         nn_scores_get_allowed_basepair (l, &bi, &bj, this->scores);

         p_bp2 = seqmatrix_get_probability(bi, pi2, sm)
            * seqmatrix_get_probability(bj, pj2, sm);

         p_bp1 = seqmatrix_get_probability(bi, pi1, sm)
            * seqmatrix_get_probability(bj, pj1, sm);

         /* for all bases */
         for (m = 0; m < alpha_size; m++)
         {

            /* for all bases */
            for (n = 0; n < alpha_size; n++)
            {
               p_bb = seqmatrix_get_probability(m, pi1 + 1, sm)
                  * seqmatrix_get_probability(n, pj1 - 1, sm);

               /* for all bases */
               for (o = 0; o < alpha_size; o++)
               {
                  /* i1 */
                  p =  p_bp2 * p_bb
                     * seqmatrix_get_probability(bpp, pj1, sm)
                     * seqmatrix_get_probability(o, pj2 + 1, sm);

                  cell_i1 += p * nn_scores_get_G_internal_1x2_loop (row,
                                                                    bpp,
                                                                    m,
                                                                    o,
                                                                    n,
                                                                    bj,
                                                                    bi,
                                                                  this->scores);

                  /* j1 */
                  p =  p_bp2 * p_bb
                     * seqmatrix_get_probability(bpp, pi1, sm)
                     * seqmatrix_get_probability(o, pj2 + 1, sm);

                  cell_j1 += p * nn_scores_get_G_internal_1x2_loop (bpp,
                                                                    row,
                                                                    m,
                                                                    o,
                                                                    n,
                                                                    bj,
                                                                    bi,
                                                                  this->scores);

                  /* i2 */
                  p =  p_bp1 * p_bb
                     * seqmatrix_get_probability(bpp,  pj2, sm)
                     * seqmatrix_get_probability(o, pj2 + 1, sm);

                  cell_i2 += p * nn_scores_get_G_internal_1x2_loop (bi,
                                                                    bj,
                                                                    m,
                                                                    o,
                                                                    n,
                                                                    bpp,
                                                                    row,
                                                                  this->scores);

               /* j2 */
               p =  p_bp1 * p_bb
                  * seqmatrix_get_probability(bpp,  pi2, sm)
                  * seqmatrix_get_probability(o, pj2 + 1, sm);
               
               cell_j2 += p * nn_scores_get_G_internal_1x2_loop (bi,
                                                                 bj,
                                                                 m,
                                                                 o,
                                                                 n,
                                                                 row,
                                                                 bpp,
                                                                 this->scores);
               }
            }
         }
      }
   }
   /* 7 bases are involved in a 1x2 loop*/
   cell_i1 = cell_i1 / 7; /* SB 08-12-12 */
   cell_j1 = cell_j1 / 7; /* SB 08-12-12 */
   cell_i2 = cell_i2 / 7; /* SB 08-12-12 */
   cell_j2 = cell_j2 / 7; /* SB 08-12-12 */
   seqmatrix_add_2_eeff (cell_i1, row, pi1, sm);
   seqmatrix_add_2_eeff (cell_j1, row, pj1, sm);
   seqmatrix_add_2_eeff (cell_i2, row, pi2, sm);
   seqmatrix_add_2_eeff (cell_j2, row, pj2, sm);

   /* design unpaired bases */
   cell_i1 = 0.0f;
   cell_j1 = 0.0f;
   cell_i2 = 0.0f;

   /* for all possible pairs */
   for (k = 0; k < allowed_bp; k++)
   {
      nn_scores_get_allowed_basepair (k, &bi, &bj, this->scores);

      p_bp1 = seqmatrix_get_probability(bi, pi1, sm)
         * seqmatrix_get_probability(bj, pj1, sm);

      /* for all possible pairs */
      for (l = 0; l < allowed_bp; l++)
      {
         nn_scores_get_allowed_basepair (l, &bi2, &bj2, this->scores);

         p_bp2 = p_bp1 * seqmatrix_get_probability(bi2,  pi2, sm)
            * seqmatrix_get_probability(bj2, pj2, sm);

         /* for all bases */
         for (m = 0; m < alpha_size; m++)
         {

            /* for all bases */
            for (n = 0; n < alpha_size; n++)
            {
               /* i1 + 1 */
               p =  p_bp2
                  * seqmatrix_get_probability(n, pj2 + 1, sm)
                  * seqmatrix_get_probability(m, pj1 - 1, sm);

               cell_i1 += p * nn_scores_get_G_internal_1x2_loop (bi,
                                                                 bj,
                                                                 row,
                                                                 n,
                                                                 m,
                                                                 bj2,
                                                                 bi2,
                                                                 this->scores);

               /* j1 - 1 */
               p =  p_bp2
                  * seqmatrix_get_probability(n, pj2 + 1, sm)
                  * seqmatrix_get_probability(m, pi1 + 1, sm);

               cell_j1 += p * nn_scores_get_G_internal_1x2_loop (bi,
                                                                 bj,
                                                                 m,
                                                                 n,
                                                                 row,
                                                                 bj2,
                                                                 bi2,
                                                                 this->scores);

               /* j2 + 1 */
               p =  p_bp2
                  * seqmatrix_get_probability(n, pj1 - 1, sm)
                  * seqmatrix_get_probability(m, pi1 + 1, sm);

               cell_i2 += p * nn_scores_get_G_internal_1x2_loop (bi,
                                                                 bj,
                                                                 m,
                                                                 row,
                                                                 n,
                                                                 bj2,
                                                                 bi2,
                                                                 this->scores);
            }
         }
      }
   }
   /* 7 bases are involved in a 1x2 loop*/
   cell_i1 = cell_i1 / 7; /* SB 08-12-12 */
   cell_j1 = cell_j1 / 7; /* SB 08-12-12 */
   cell_i2 = cell_i2 / 7; /* SB 08-12-12 */
   seqmatrix_add_2_eeff (cell_i1, row, pi1 + 1, sm);
   seqmatrix_add_2_eeff (cell_j1, row, pj1 - 1, sm);
   seqmatrix_add_2_eeff (cell_i2, row, pj2 + 1, sm);
}

static __inline__ void
scmf_rna_opt_calc_int11 (const unsigned long row,
                         unsigned long allowed_bp,
                         unsigned long alpha_size,
                         unsigned long pi1, unsigned long pj1,
                         unsigned long pi2, unsigned long pj2,
                         SeqMatrix* sm,
                         Scmf_Rna_Opt_data* this)
{
   unsigned long k, l, m, n;

   float cell_i1 = 0.0f;
   float cell_j1 = 0.0f;
   float cell_i2 = 0.0f;
   float cell_j2 = 0.0f;
   char bpp, bi, bj, bi2, bj2;
   float p_bp2, p_bp1, p;

   /* for all allowed pairs */
   for (k = 0; this->bp_allowed[row][k] != 0; k++)
   {
      bpp = (char) (this->bp_allowed[row][k] - 1);
      
      /* for all possible pairs */
      for (l = 0; l < allowed_bp; l++)
      {
         nn_scores_get_allowed_basepair (l, &bi, &bj, this->scores);
         
         p_bp2 = seqmatrix_get_probability(bi, pi2, sm)
            * seqmatrix_get_probability(bj, pj2, sm);
         
         p_bp1 = seqmatrix_get_probability(bi, pi1, sm)
            * seqmatrix_get_probability(bj, pj1, sm);
         
         /* for all bases */
         for (m = 0; m < alpha_size; m++)
         {
            /* for all bases */
            for (n = 0; n < alpha_size; n++)
            {
               p = seqmatrix_get_probability(bpp, pj1, sm)
                  * p_bp2
                  * seqmatrix_get_probability(m, pi1 + 1, sm)
                  * seqmatrix_get_probability(n, pj1 - 1, sm);
               
               cell_i1 += p * nn_scores_get_G_internal_1x1_loop (row, bpp,
                                                                 m, n,
                                                                 bi, bj,
                                                                 this->scores);

               p = seqmatrix_get_probability(bpp, pi1, sm)
                  * p_bp2
                  * seqmatrix_get_probability(m, pi1 + 1, sm)
                  * seqmatrix_get_probability(n, pj1 - 1, sm);
               
               cell_j1 += p * nn_scores_get_G_internal_1x1_loop (bpp, row,
                                                                 m, n,
                                                                 bi, bj,
                                                                 this->scores);
               
               p = p_bp1
                  * seqmatrix_get_probability(bpp, pj2, sm)
                  * seqmatrix_get_probability(m, pi1 + 1, sm)
                  * seqmatrix_get_probability(n, pj1 - 1, sm);
               
               cell_i2 += p * nn_scores_get_G_internal_1x1_loop (bi, bj,
                                                                 m, n,
                                                                 row, bpp,
                                                                 this->scores);
               
               p = p_bp1
                  * seqmatrix_get_probability(bpp, pi2, sm)
                  * seqmatrix_get_probability(m, pi1 + 1, sm)
                  * seqmatrix_get_probability(n, pj1 - 1, sm);
               
               cell_j2 += p * nn_scores_get_G_internal_1x1_loop (bi, bj,
                                                                 m, n,
                                                                 bpp, row,
                                                                 this->scores);
            }
         }
      }
   }
   /* each of the 6 bases involved in this loop gets an energy contribution */
   cell_i1 = cell_i1 / 6; /* SB 08-12-12 */
   cell_j1 = cell_j1 / 6; /* SB 08-12-12 */
   cell_i2 = cell_i2 / 6; /* SB 08-12-12 */
   cell_j2 = cell_j2 / 6; /* SB 08-12-12 */
   seqmatrix_add_2_eeff (cell_i1, row, pi1, sm);
   seqmatrix_add_2_eeff (cell_j1, row, pj1, sm);
   seqmatrix_add_2_eeff (cell_i2, row, pi2, sm);
   seqmatrix_add_2_eeff (cell_j2, row, pj2, sm);

   /* design unpaired bases */
   cell_i1 = 0.0f;
   cell_j1 = 0.0f;
   /* for all possible pairs */
   for (k = 0; k < allowed_bp; k++)
   {
      nn_scores_get_allowed_basepair (k, &bi, &bj, this->scores);
      
      p_bp1 = seqmatrix_get_probability(bi, pi1, sm)
         * seqmatrix_get_probability(bj, pj1, sm);
      
      /* for all possible pairs */
      for (l = 0; l < allowed_bp; l++)
      {
         nn_scores_get_allowed_basepair (l, &bi2, &bj2, this->scores);
         
         p_bp2 = p_bp1
            * seqmatrix_get_probability(bi2, pi2, sm)
            * seqmatrix_get_probability(bj2, pj2, sm);
         
         /* for all bases */
         for (m = 0; m < alpha_size; m++)
         {
            cell_i1 += p_bp2 * seqmatrix_get_probability(m, pj1 - 1, sm)
               * nn_scores_get_G_internal_1x1_loop (bi, bj,
                                                    row, m,
                                                    bi2, bj2,
                                                    this->scores);
            
            cell_j1 += p_bp2 * seqmatrix_get_probability(m, pi1 + 1, sm)
               * nn_scores_get_G_internal_1x1_loop (bi, bj,
                                                    m, row,
                                                    bi2, bj2,
                                                    this->scores);
         }
      }
   }
   /* each of the 6 bases involved in this loop gets an energy contribution */
   cell_i1 = cell_i1 / 6; /* SB 08-12-12 */
   cell_j1 = cell_j1 / 6; /* SB 08-12-12 */
   seqmatrix_add_2_eeff (cell_i1, row, pi1 + 1, sm);
   seqmatrix_add_2_eeff (cell_j1, row, pj1 - 1, sm);
}

static void
scmf_rna_opt_calc_internals (const unsigned long row,
                             const unsigned long loop,
                             SeqMatrix* sm,
                             Scmf_Rna_Opt_data* this)
{
   unsigned long pi1, pj1, size1;
   unsigned long pi2, pj2, size2;
   unsigned long allowed_bp = nn_scores_no_allowed_basepairs (this->scores);
   unsigned long alpha_size = alphabet_size (this->sigma);

   /* fetch loop geometry */
   rna_secstruct_get_geometry_internal (&pi1, &pj1, &pi2, &pj2, &size1, &size2,
                                        loop,
                                        this->rna);

   /*mfprintf (stderr, "%lu: i1: %lu j1: %lu size1: %lu "
             "i2: %lu j2: %lu size2: %lu\n", loop,
             i1, j1, size1, i2, j2, size2);*/

   if ((size1 == 1) && (size2 == 1))
   {
      /* 1x1 internal loop */
      scmf_rna_opt_calc_int11 (row, allowed_bp, alpha_size, 
                               pi1, pj1, pi2, pj2, sm, this);
   }
   else if ((size1 == 1) && (size2 == 2))
   {
      /* 1x2 internal loop */
      scmf_rna_opt_calc_int12 (row, allowed_bp, alpha_size, 
                               pi1, pj1, pi2, pj2, sm, this);
   }
   else if ((size1 == 2) && (size2 == 1))
   {
      /* 2x1 internal loop */
      scmf_rna_opt_calc_int12 (row, allowed_bp, alpha_size, 
                               pj2, pi2, pj1, pi1, sm, this);
   }
   else if ((size1 == 2) && (size2 == 2))
   {
      /* 2x2 internal loop */
      scmf_rna_opt_calc_int22 (row, allowed_bp, alpha_size, 
                               pi1, pj1, pi2, pj2, sm, this);
   }
   else
   {
      /* generic internal loop */
      scmf_rna_opt_calc_internal (row, allowed_bp, alpha_size, 
                                  pi1, pj1, pi2, pj2, sm, this);      
   }
}

static float
scmf_rna_opt_get_interaction_energy (char state,
                                     const unsigned long site,
                                     const unsigned long partner,
                                     unsigned long allowed_bp,
                                     const Scmf_Rna_Opt_data* this,
                                     SeqMatrix* sm)
{
   float pe = 0.0f;
   unsigned long abp, bp;
   char bip, bjm, bj;
   char* p5;
   char* p3;
   unsigned long pos_j, pos_i;

   if (site < partner)
   {
      p5 = &state;
      p3 = &bj;
      pos_i = site;
      pos_j = partner;
   }
   else
   {
      p3 = &state;
      p5 = &bj;
      pos_j = site;
      pos_i = partner;
   }

   /* for all possible pairing partners */
   for (abp = 0; this->bp_allowed[(int)state][abp] != 0; abp++)
   {
      bj = (char) (this->bp_allowed[(int)state][abp] - 1);
      
      /* for all stacking bp */
      for (bp = 0; bp < allowed_bp; bp++)
      {
         nn_scores_get_allowed_basepair (bp, &bjm, &bip, this->scores);

         pe += (nn_scores_get_G_stack (*p5, *p3, bjm, bip, this->scores)
                * seqmatrix_get_probability (bj, partner/* pos_j */, sm)
                * seqmatrix_get_probability (bjm, pos_j - 1, sm)
                * seqmatrix_get_probability (bip, pos_i + 1, sm));
      }
   }
   
   return pe;
}

static void
scmf_rna_opt_calc_upstream_cont (const unsigned long state,
                                 const unsigned long j,
                                 const unsigned long allowed_bp,
                                 const Scmf_Rna_Opt_data* this,
                                 const SeqMatrix* sm)
{
   unsigned long abp, bp;
   char bip, bjm, bj;

   for (abp = 0; this->bp_allowed[state][abp] != 0; abp++)
   {
      bj = (char) (this->bp_allowed[state][abp] - 1);

      for (bp = 0; bp < allowed_bp; bp++)
      {
         nn_scores_get_allowed_basepair (bp, &bjm, &bip, this->scores);
         
         this->en_neg[(int)bip] -= (nn_scores_get_G_stack (state, bj, bjm, bip,
                                                           this->scores)
                                    * seqmatrix_get_probability (bjm, (j-1),sm)
                                    * seqmatrix_get_probability (bj, j, sm));
      }
   }
}

static void
scmf_rna_opt_calc_downstream_cont (const unsigned long state,
                                   const unsigned long i,
                                   const unsigned long allowed_bp,
                                   const Scmf_Rna_Opt_data* this,
                                   const SeqMatrix* sm)
{
   unsigned long abp, bp;
   char bi, bip, bjm;

   for (abp = 0; this->bp_allowed[state][abp] != 0; abp++)
   {
      bi = (char)(this->bp_allowed[state][abp] - 1);

      for (bp = 0; bp < allowed_bp; bp++)
      {
         nn_scores_get_allowed_basepair (bp, &bip, &bjm, this->scores);

      this->en_neg_35[0][(int)bjm] += (nn_scores_get_G_stack (bi,state,bjm,bip,
                                                              this->scores)
                                       * seqmatrix_get_probability (bi, i, sm)
                                 * seqmatrix_get_probability (bip, i + 1, sm));
      }
   }
}

static void
scmf_rna_opt_calc_neg_loop (const unsigned long state,
                            const unsigned long n_sites,
                            const unsigned long allowed_bp,
                            const unsigned long alpha_size,
                            const Scmf_Rna_Opt_data* this,
                            SeqMatrix* sm)
{
   /* called for each state seperately */
   unsigned long j, abp, bp, paired_2;
   char bj, bip, bjm;
   float prob;

   /* init upstream direction */
   /* Idea: For the first state, add up all negative interactions. But only
            count energies multiplied with probabilities INDEPENDENT of the
            first position (that is j, j-1 since cur.state is i, i+1 is not
            independent). This gives us a quasi force-field which can be
            iterated to produce neg.design terms for all other sites. This
            turns the naive quadratic approach into linear time. */

   memset(this->en_neg, 0, alpha_size * sizeof (*(this->en_neg)));
   memset(this->en_neg_35[0], 0, alpha_size * sizeof(**(this->en_neg_35)));

   /* walk over all pairing partners of given state */
   for (abp = 0; this->bp_allowed[state][abp] != 0; abp++)
   {
      bj = (char) (this->bp_allowed[state][abp] - 1);
      
      /* stack with all pairs allowed in our alphabet */
      for (bp = 0; bp < allowed_bp; bp++)
      {
         nn_scores_get_allowed_basepair (bp, &bjm, &bip, this->scores);
         
         /* iterate all but the first site to be the j-side in a pair */
         prob = 0.0f;
         for (j = 1; j < n_sites; j++)
         {
            prob += (seqmatrix_get_probability (bjm, (j - 1), sm)
                   * seqmatrix_get_probability (bj,   j,      sm));
         }

         this->en_neg[(int)bip] += (nn_scores_get_G_stack (state, bj, bjm, bip,
                                                          this->scores) * prob);
      }
   }

   /* iterate over all sites */
   /* treat site 0 as special case */
   if (!seqmatrix_is_col_fixed (0, sm))
   {
      prob = 0.0f;
      for (abp = 0; abp < alpha_size; abp++)
      {
         prob += (this->en_neg[abp] * seqmatrix_get_probability (abp, 1, sm));
      }
      
      if ((paired_2 = rna_base_pairs_with (0, this->rna)) != NOT_PAIRED)
      {
         prob -= scmf_rna_opt_get_interaction_energy ((char) state, 0, paired_2,
                                                      allowed_bp,
                                                      this, sm);
      }
   
      /*mfprintf (stdout, "%lu,%lu: %f\n", state, 0L, prob);*/
      seqmatrix_add_2_eeff ((prob/n_sites) * -1.0f * this->neg_scale,state, 0,
                            sm);
   }

   scmf_rna_opt_calc_upstream_cont (state, 1, allowed_bp, this, sm);
   scmf_rna_opt_calc_downstream_cont (state, 0, allowed_bp, this, sm);

   for (j = 1; j < (n_sites - 1); j++)
   {
      if (!seqmatrix_is_col_fixed (j, sm))
      {
         /* create current neg. design term */
         prob = 0.0f;
         
         /* add all contributions of possible stacking pairing parnters */
         for (abp = 0; abp < alpha_size; abp++)
         {
            prob += (this->en_neg[abp] 
                     * seqmatrix_get_probability (abp, j + 1, sm));
            prob += (this->en_neg_35[0][abp]
                     * seqmatrix_get_probability (abp, (j - 1), sm));
         }
         
         /*   remove possible interaction */      
         if ((paired_2 = rna_base_pairs_with (j, this->rna)) != NOT_PAIRED)
         {
            prob -= scmf_rna_opt_get_interaction_energy ((char) state, j,
                                                         paired_2,
                                                         allowed_bp,
                                                         this, sm);
         }
         
         
         /* if we pair, we have stacks (i, i+1, site-1, site) and
            (site, site+1, j-1, j) counted without need */
         
         /*   store */
         /*mfprintf (stdout, "%lu,%lu: %f\n", state, j, prob);*/
         seqmatrix_add_2_eeff ((prob/n_sites) * -1.0f * this->neg_scale, state,
                               j,
                               sm);
      }

      /*   update downstream */
      scmf_rna_opt_calc_downstream_cont (state, j, allowed_bp, this, sm);

      /*   update upstream */
      scmf_rna_opt_calc_upstream_cont (state, j + 1, allowed_bp, this, sm);

      /* exclude one sonderfall (0 and last) */
   }

   if (!seqmatrix_is_col_fixed (j, sm))
   {
      /* treat last site as special case */
      prob = 0.0f;
      for (abp = 0; abp < alpha_size; abp++)
      {
         prob += (this->en_neg_35[0][abp]
                  * seqmatrix_get_probability (abp, (j - 1), sm));
      }
      
      if ((paired_2 = rna_base_pairs_with (j, this->rna)) != NOT_PAIRED)
      {
         prob -= scmf_rna_opt_get_interaction_energy ((char) state, j,
                                                      paired_2, allowed_bp,
                                                      this, sm);
      }
      
      /*mfprintf (stdout, "%lu,%lu: %f\n", state, j, prob);*/
      seqmatrix_add_2_eeff ((prob/n_sites) * -1.0f * this->neg_scale,state, j,
                            sm);
   }
}

static void
scmf_rna_opt_calc_het_term (const unsigned long state,
                            const unsigned long n_sites,
                            const Scmf_Rna_Opt_data* this,
                            SeqMatrix* sm)
{
   unsigned long s_c/* , s_h */;
   float het;
   /* float het_count = 0.0f; */

   /* do it like gundolf - START */
   /***************************************/
   /* first try: additive, whole sequence */
   /***************************************/
   het = 0.0f;
   /* calc first pos */
   for (s_c = 0; s_c < n_sites; s_c++)
   {
      het += seqmatrix_get_probability (state, s_c, sm);
   }

   /* propagate on whole sequence */
   for (s_c = 0; s_c < n_sites; s_c++)
   {
      seqmatrix_add_2_eeff (((het - seqmatrix_get_probability (state, s_c, sm))
                             / n_sites) * this->het_scale, state, s_c, sm);
   }


   /***************************************/
   /*      2nd try: additive, window      */
   /***************************************/

   /* do it like gundolf - END */

/*    /\* for each site *\/ */
/*    for (s_c = 0; s_c < n_sites; s_c++) */
/*    { */
/*       het = 0.0f; */

/*       /\* pay attention to all sites before the current *\/ */
/*       for (s_h = 0; s_h < s_c; s_h++) */
/*       { */
/*          if (rna_base_pairs_with (s_c, this->rna) != s_h) */
/*          { */
/*             het += (seqmatrix_get_probability (state, s_h, sm) */
/*                     * expf(this->het_rate * (s_c - (s_h + 1)))); */
/*             het_count += expf ( (this->het_rate* (s_c - (s_h + 1)))); */
/*          } */
/*       } */

/*       /\* pay attention to all remaining sites *\/ */
/*       for (s_h = s_c + 1; s_h < n_sites; s_h++) */
/*       { */
/*          if (rna_base_pairs_with (s_c, this->rna) != s_h) */
/*          { */
/*             het += (seqmatrix_get_probability (state, s_h, sm) */
/*                         * expf(this->het_rate * (s_h - (s_c + 1)))); */
/*             het_count += (expf (this->het_rate * (s_h - (s_c + 1)))); */
/*          } */
/*       } */

/*       /\* store contribution *\/ */
/*       /\*mfprintf (stdout, "%lu: %f\n", s_c, this->het_scale);*\/ */
/*       seqmatrix_add_2_eeff ((het/\* / het_count *\/) * this->het_scale, state, s_c, sm); */
/*    } */
}

/** @brief SCMF simulation function.
 *
 * This is the substitute for the column iteration function of a SCMF
 * simulation. Instead of iterating the columns of a sequence matrix we iterate
 * over the structural components of a RNA secondary structure.
 *
 * @params[in] sm Sequence matrix.
 * @params[in] t temperature.
 * @params[in] sco Data.
 */
int
scmf_rna_opt_calc_col_nn (SeqMatrix* sm,
                          const float t,
                          void* sco)
{
   /*unsigned long i;*/
   int error = 0;
   Scmf_Rna_Opt_data* this;
   unsigned long n;
   unsigned long i;
   unsigned long n_states;
   unsigned long n_sites;
   unsigned long r, c;
   /* unsigned long pi, c_neg, abp; */
   /* char b_i, b_ip, b_jm, b_j; */
   float cell;
   /* float prob; */
   unsigned long allowed_bp;
   unsigned long alpha_size;

   assert (sm);
   assert (sco);

   this = (Scmf_Rna_Opt_data*) sco;

   allowed_bp = nn_scores_no_allowed_basepairs (this->scores);
   alpha_size = alphabet_size (this->sigma);

   seqmatrix_set_eeff_matrix_zero (sm);

   /* iterate all bases (n_states) */
   n_states = seqmatrix_get_rows (sm);
   n_sites = seqmatrix_get_width (sm);
   r = 0;
   while ((r < n_states) && (!error))
   {
      /*seqmatrix_print_2_stdout (2, sm);*/
      /* process structure components */
      /* external loop */
      scmf_rna_opt_calc_ext_loop (r, sm, this);

      /* stacking pairs */
      n = rna_secstruct_get_noof_stacks (this->rna);
      for (i = 0; i < n; i++)
      {
         scmf_rna_opt_calc_stack (r, i, sm, this);
      }

      /* bulge loops */
      n = rna_secstruct_get_noof_bulges (this->rna);
      for (i = 0; i < n; i++)
      {
         scmf_rna_opt_calc_bulge (r, i, sm, this);
      }

      /* internal loops */
      n = rna_secstruct_get_noof_internals (this->rna);
      for (i = 0; i < n; i++)
      {
         scmf_rna_opt_calc_internals (r, i, sm, this);
      }      

      /* hairpin loops */
      n = rna_secstruct_get_noof_hairpins (this->rna);
      for (i = 0; i < n; i++)
      {
         scmf_rna_opt_calc_hairpin (r, i, sm, this);
      }

      /* multiloops */
      n = rna_secstruct_get_noof_multiloops (this->rna);
      for (i = 0; i < n; i++)
      {
         scmf_rna_opt_calc_multi_loop (r, i, sm, this);
      }

      /* calc. neg. design term, iteratevily */
      scmf_rna_opt_calc_neg_loop (r, n_sites, allowed_bp, alpha_size, this, sm);

      /* heterogenity term */
      scmf_rna_opt_calc_het_term (r, n_sites, this, sm);

      for (c = 0; c < n_sites; c++)
      {
         if (!seqmatrix_is_col_fixed (c, sm))
         {
/*             mprintf ("CM(%lu,%lu) = %.2f\n", r, c, */
/*               seqmatrix_get_eeff (r, c, sm)); */
            cell = expf((-1.0f) * seqmatrix_get_eeff (r, c, sm)
                        / (t * seqmatrix_get_gas_constant(sm)));
            seqmatrix_set_eeff (cell, r, c, sm);
            /*mprintf (" %.2f\n", seqmatrix_get_eeff (r, c, sm));*/
         }         
      }

      r++;
   }
   /* mfprintf (stderr, "\n"); */


   return error;
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
