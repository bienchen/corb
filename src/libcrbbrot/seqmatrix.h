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
   ERR_SM_WRITE,          /* problems on proper writing to a file */
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

float
seqmatrix_get_eeff (const unsigned long, const unsigned long,
                    const SeqMatrix*);

float
seqmatrix_get_gas_constant (const SeqMatrix*);

/********************************   Altering   ********************************/

int
seqmatrix_init (const unsigned long,
                const unsigned long,
                SeqMatrix*,
                const char*, const int);

#define SEQMATRIX_INIT(A, B, SM)                         \
   seqmatrix_init (A, B, SM, __FILE__, __LINE__)


void
seqmatrix_set_eeff_matrix_zero (SeqMatrix*);

void
seqmatrix_set_eeff (const float,
                    const unsigned long,
                    const unsigned long,
                    SeqMatrix*);

void
seqmatrix_add_2_eeff (const float,
                                const unsigned long,
                                const unsigned long,
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
seqmatrix_set_pre_col_iter_hook(int (*pre_col_iter_hook) (void*, SeqMatrix*),
                                SeqMatrix* sm);

void
seqmatrix_set_fixed_site_hook (int (*fixed_site_hook) (void*, unsigned long,
                                                       SeqMatrix*),
                               SeqMatrix* sm);

void
seqmatrix_set_fixing_site_hook (int (*fixing_site_hook) (void*, unsigned long,
                                                         SeqMatrix*),
                                SeqMatrix* sm);

void
seqmatrix_set_get_seq_string (char* (*get_seq_string) (void*), SeqMatrix* sm);

int
seqmatrix_fix_col (const unsigned long,
                   const unsigned long,
                   void*,
                   SeqMatrix*);

int
seqmatrix_simulate_scmf (unsigned long,
                         const float,
                         const float,
                         const float,
                         const float,
                         const float,
                         const float,
                         /*const float,*/
                         const float,
                         GFile*,
                         GFile*,
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
                      const float,
                      /* const float, */
                      SeqMatrix*,
                      void*);

int
seqmatrix_collate_mv (const SeqMatrix*, void*);

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
