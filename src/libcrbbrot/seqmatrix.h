/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbrot/seqmatrix.h
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


#include <libcrbrna/crbrna.h>


#ifdef __cplusplus
extern "C" {
#endif

#ifndef SEQMATRIX_H
#define SEQMATRIX_H

enum seqmatrix_retvals{
   ERR_SM_ALLOC = 1,      /* (re)allocation problems */
   ERR_SM_PRINT,          /* problems on proper printing */
};

typedef struct SeqMatrix SeqMatrix;


/**********************   Constructors and destructors   **********************/
SeqMatrix*
seqmatrix_new (const char*, const int);

#define SEQMATRIX_NEW seqmatrix_new (__FILE__, __LINE__)

void
seqmatrix_delete (SeqMatrix*);


/*********************************   Access   *********************************/

bool
seqmatrix_is_col_fixed (const unsigned long, const SeqMatrix*);

unsigned long
seqmatrix_get_width (const SeqMatrix*);

/********************************   Altering   ********************************/

int
seqmatrix_init (const unsigned long*,
                const unsigned long,
                SeqMatrix*);

void
seqmatrix_fix_col (const unsigned long,
                   const unsigned long,
                   SeqMatrix*);

int
sequence_matrix_simulate_scmf (const unsigned long,
                               const float,
                               SeqMatrix*,
                               float**);

int
sequence_matrix_simulate_scmf_nn (const unsigned long,
                                  const float,
                                  SeqMatrix*,
                                  const NN_scores*, const Alphabet*);

/*********************************   Output   *********************************/

void
seqmatrix_fprintf_sequence (FILE*, SeqMatrix*);

void
seqmatrix_printf_sequence (SeqMatrix*);

int
seqmatrix_collate_is (const float,
                      const unsigned long,
                      const float,
                      float**,
                      SeqMatrix*,
                      const Alphabet*);

int
seqmatrix_collate_is_nn (const float,
                         const unsigned long,
                         const float,
                         const NN_scores*,
                         SeqMatrix*,
                         const Alphabet*);

int
seqmatrix_collate_mv (SeqMatrix*, const Alphabet*);

void
seqmatrix_fprintf (FILE*, const int, const SeqMatrix*);

void
seqmatrix_print_2_stdout (const int, const SeqMatrix*);

void
seqmatrix_print_2_stderr (const int, const SeqMatrix*);

#endif /* SEQMATRIX_H */

#ifdef __cplusplus
}
#endif
