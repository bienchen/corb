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

/*__inline__*/ bool
seqmatrix_is_row_fixed (unsigned long, SeqMatrix*);


/********************************   Altering   ********************************/

int
seqmatrix_init (const unsigned long*,
                const unsigned long,
                const PresetArray*,
                SeqMatrix*);

int
sequence_matrix_simulate_scmf (const unsigned long,
                               const float,
                               SeqMatrix*,
                               float**);

/*********************************   Output   *********************************/

void
seqmatrix_fprintf_sequence (FILE*, SeqMatrix*);

void
seqmatrix_printf_sequence (SeqMatrix*);

int
seqmatrix_collate_is (float, unsigned long, float, float**, SeqMatrix*);

int
seqmatrix_collate_mv (SeqMatrix*);

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
