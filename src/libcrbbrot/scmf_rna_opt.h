/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
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
scmf_rna_opt_data_delete_nn (Scmf_Rna_Opt_data*);

void
scmf_rna_opt_data_delete_nussi (Scmf_Rna_Opt_data*);

void
scmf_rna_opt_data_set (NN_scores*, Alphabet*, Scmf_Rna_Opt_data*);

Scmf_Rna_Opt_data*
scmf_rna_opt_data_new_nn (const unsigned long, const char*, const int);

#define SCMF_RNA_OPT_DATA_NEW_NN(L) \
   scmf_rna_opt_data_new_nn (L, __FILE__, __LINE__)

Scmf_Rna_Opt_data*
scmf_rna_opt_data_new_nussi (const unsigned long, const char*, const int);

#define SCMF_RNA_OPT_DATA_NEW_NUSSI(L) \
   scmf_rna_opt_data_new_nussi (L, __FILE__, __LINE__)

int
scmf_rna_opt_data_transform_row_2_base (const unsigned long,
                                        const unsigned long,
                                        void*);

/*Alphabet*
  scmf_rna_opt_data_get_alphabet (Scmf_Rna_Opt_data*);*/

char*
scmf_rna_opt_data_get_seq (Scmf_Rna_Opt_data*);

float
scmf_rna_opt_calc_nussinov (const unsigned long, const unsigned long,
                           void*,
                           SeqMatrix*);

float
scmf_rna_opt_calc_nn (const unsigned long, const unsigned long,
                     void*,
                     SeqMatrix*);
   
#endif /* SCMF_RNA_OPT_H */

#ifdef __cplusplus
}
#endif
