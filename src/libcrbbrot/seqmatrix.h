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

bool
seqmatrix_is_col_fixed (const unsigned long, const SeqMatrix*);

unsigned long
seqmatrix_get_width (const SeqMatrix*);

unsigned long
seqmatrix_get_rows (const SeqMatrix*);

float
seqmatrix_get_probability (const unsigned long, const unsigned long,
                            const SeqMatrix*);

/********************************   Altering   ********************************/

int
seqmatrix_init (const unsigned long*,
                const unsigned long,
                const unsigned long,
                SeqMatrix*,
                const char*, const int);

#define SEQMATRIX_INIT(A, B, C, SM) \
   seqmatrix_init (A, B, C, SM, __FILE__, __LINE__)

void
seqmatrix_set_cell (const float, const unsigned long, const unsigned long,
                    SeqMatrix*);

void
seqmatrix_set_gas_constant (const float, SeqMatrix*);

void
seqmatrix_set_func_calc_eeff_col (int (*calc_eeff_col) (SeqMatrix*,
                                                        const float,
                                                        void*),
                                                        SeqMatrix*);

void
seqmatrix_set_func_calc_cell_energy (float (*calc_cell_energy)
                                     (const unsigned long, const unsigned long,
                                      void*,
                                      SeqMatrix*),
                                     SeqMatrix*);

void
seqmatrix_set_func_calc_eeff_row (int (*calc_eeff_row) (const unsigned long,
                                                        SeqMatrix*,
                                                        const float,
                                                        void*),
                                  SeqMatrix*);

void
seqmatrix_set_transform_row (int (*transform_row) (const unsigned long,
                                                   const unsigned long,
                                                   void*),
                             SeqMatrix*);

void
seqmatrix_fix_col (const unsigned long,
                   const unsigned long,
                   SeqMatrix*);

int
seqmatrix_simulate_scmf (const unsigned long,
                         const float,
                         const float,
                         const float,
                         const float,
                         const float,
                         SeqMatrix*,
                         void*);

/*********************************   Output   *********************************/

int
seqmatrix_collate_is (const float,
                      const unsigned long,
                      const float,
                      const float,
                      const float,
                      const float,
                      const float,
                      SeqMatrix*,
                      void*);

int
seqmatrix_collate_mv (SeqMatrix*, void*);

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
