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
 *  @file libcrbbrot/scmf_rna_opt.h
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

#include <libcrbrna/crbrna.h>
#include "seqmatrix.h"


#ifdef __cplusplus
extern "C" {
#endif

#ifndef SCMF_RNA_OPT_H
#define SCMF_RNA_OPT_H


typedef struct Scmf_Rna_Opt_data Scmf_Rna_Opt_data;

Scmf_Rna_Opt_data*
scmf_rna_opt_data_new (const char*, const int);

#define SCMF_RNA_OPT_DATA_NEW scmf_rna_opt_data_new (__FILE__, __LINE__)

void
scmf_rna_opt_data_delete (Scmf_Rna_Opt_data*);

Scmf_Rna_Opt_data*
scmf_rna_opt_data_new_init (const char*,
                            const unsigned long,
                            const char*,
                            const unsigned long,
                            float,
                            unsigned int,
                            const char*, const int);

#define SCMF_RNA_OPT_DATA_NEW_INIT(S, L, A, G, H, I)                      \
   scmf_rna_opt_data_new_init (S, L, A, G, H, I, __FILE__, __LINE__)

/* int */
/* scmf_rna_opt_data_init_negative_design_energies (void*, */
/*                                                  SeqMatrix*); */

int
scmf_rna_opt_data_secstruct_init (Scmf_Rna_Opt_data*);

int
scmf_rna_opt_data_update_neg_design_energy (void*, unsigned long, SeqMatrix*);

int
scmf_rna_opt_data_init_negative_design_energies_alt (void*,
                                                 SeqMatrix*);

int
scmf_rna_opt_data_transform_row_2_base (const unsigned long,
                                        const unsigned long,
                                        void*);

void
scmf_rna_opt_data_set_scales (float, float, Scmf_Rna_Opt_data*);

void
scmf_rna_opt_data_set_scores (void*, Scmf_Rna_Opt_data*);

void
scmf_rna_opt_data_set_bp_allowed (char**, Scmf_Rna_Opt_data*);

void
scmf_rna_opt_data_set_het_window (const long, Scmf_Rna_Opt_data*);

Alphabet*
scmf_rna_opt_data_get_alphabet (Scmf_Rna_Opt_data*);

char*
scmf_rna_opt_data_get_seq (Scmf_Rna_Opt_data*);

unsigned long
scmf_rna_opt_data_get_rna_size (Scmf_Rna_Opt_data*);

float
scmf_rna_opt_calc_nussinov (const unsigned long, const unsigned long,
                           void*,
                           SeqMatrix*);

float
scmf_rna_opt_calc_simplenn (const unsigned long, const unsigned long,
                     void*,
                     SeqMatrix*);

int
scmf_rna_opt_calc_col_nn (SeqMatrix*, const float, void*);

#endif /* SCMF_RNA_OPT_H */

#ifdef __cplusplus
}
#endif
