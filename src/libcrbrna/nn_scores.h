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
 *  @file libcrbrna/nn_scores.h
 *
 *  @brief Neares neighbour Model for evaluating RNA secondary structures
 *
 *  Module: nn_scores
 *
 *  Library: crbrna
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-07-03
 *
 *
 *  Revision History:
 *         - 2008Jul03 bienert: created
 *
 */


#ifdef __cplusplus
extern "C" {
#endif

#ifndef NN_SCORES_H
#define NN_SCORES_H

#include "alphabet.h"

/* indeces for multiloop arrays */
enum stem_indeces {
   P5_Strand = 0,
   P3_Strand,
   No_Of_Strands
};

enum dangle_indeces {
   P5_Dangle = 0,
   P3_Dangle,
   Ne_Dangle,
   No_Of_Dangles
};

typedef struct NN_scores NN_scores;
   

/**********************   Constructors and destructors   **********************/

NN_scores*
nn_scores_new (const char*, const int);

#define NN_SCORES_NEW nn_scores_new (__FILE__, __LINE__)

NN_scores*
nn_scores_new_init (float, Alphabet*, const char*, const int);

#define NN_SCORES_NEW_INIT(A, B) nn_scores_new_init (A, B, __FILE__, __LINE__)

void
nn_scores_delete (NN_scores*);


/*********************************   Access   *********************************/

unsigned long
nn_scores_get_size_tetra_loop (const NN_scores*);

unsigned long
nn_scores_get_size_tetra_loop_full (const NN_scores*);

unsigned long
nn_scores_get_no_of_tetra_loops (const NN_scores*);

const char*
nn_scores_get_tetra_loop (const unsigned long, const NN_scores*);

void
nn_scores_get_allowed_basepair (unsigned, char*, char*, const NN_scores*);

float
nn_scores_get_G_non_gc_penalty_for_bp (const int, const int,
                                       const NN_scores*);

float
nn_scores_get_G_dangle5 (const int, const int, const int, const NN_scores*);

float
nn_scores_get_G_dangle3 (const int, const int, const int, const NN_scores*);

float
nn_scores_get_G_extloop_multiloop (const char*,
                                   const unsigned long,
                                   const unsigned long,
                                   unsigned long (*stems)[No_Of_Strands],
                                   const unsigned long,
                                   unsigned long (*dangle5)[No_Of_Dangles],
                                   const unsigned long,
                                   unsigned long (*dangle3)[No_Of_Dangles],
                                   const bool,
                                   const NN_scores*);

float
nn_scores_get_G_tetra_loop (const char*,
                            const unsigned long,
                            const NN_scores*);
float
nn_scores_get_G_hairpin_mismatch (const int,
                                  const int,
                                  const int,
                                  const int,
                                  const unsigned long,
                                  const NN_scores*);

float
nn_scores_get_G_hairpin_loop (const char*,
                              const unsigned long,
                              const unsigned long,
                              const unsigned long,
                              const NN_scores*);

float
nn_scores_get_G_bulge_stack (const int, const int, const int, const int,
                             const unsigned long,
                             const NN_scores*);

float
nn_scores_get_G_bulge_loop (const int, const int,
                            const int, const int,
                            const unsigned long,
                            const NN_scores*);

float
nn_scores_get_G_stack (const int, const int, const int, const int,
                       const NN_scores*);

float
nn_scores_get_G_internal_2x2_loop (const int,
                                   const int,
                                   const int,
                                   const int,
                                   const int,
                                   const int,
                                   const int,
                                   const int,
                                   const NN_scores*);

float
nn_scores_get_G_internal_1x2_loop (const int,
                                   const int,
                                   const int,
                                   const int,
                                   const int,
                                   const int,
                                   const int,
                                   const NN_scores*);

float
nn_scores_get_G_internal_1x1_loop (const int,
                                   const int,
                                   const int,
                                   const int,
                                   const int,
                                   const int,
                                   const NN_scores*);

float
nn_scores_get_G_mismatch_interior (const int,
                                   const int,
                                   const int,
                                   const int,
                                   const NN_scores*);

float
nn_scores_get_G_internal_loop (const char*,
                               const unsigned long,
                               const unsigned long,
                               const unsigned long,
                               const unsigned long,
                               const unsigned long,
                               const unsigned long,
                               const NN_scores*);

float
nn_scores_get_G_mm_stack (const char, const char, const char, const char,
                          const NN_scores*);

/*********************************    Size    *********************************/

unsigned long
nn_scores_no_allowed_basepairs (const NN_scores*);

/*********************************   Output   *********************************/

void
nn_scores_fprintf_bp_allowed (FILE*, const NN_scores*, const Alphabet*);

void
nn_scores_fprintf_bp_idx (FILE*, const NN_scores*, const Alphabet*);

void
nn_scores_fprintf_G_stack (FILE*, const NN_scores*, const Alphabet*);

void
nn_scores_fprintf_mm_G_stack (FILE*, const NN_scores*, const Alphabet*);

void
nn_scores_fprintf_G_hairpin_loop (FILE*, const NN_scores*);

void
nn_scores_fprintf_G_hairpin_loop (FILE*, const NN_scores*);

void
nn_scores_fprintf_G_mismatch_hairpin (FILE*, const NN_scores*, const Alphabet*);

void
nn_scores_fprintf_G_bulge_loop (FILE*, const NN_scores*);

void
nn_scores_fprintf_non_gc_penalty_for_bp(FILE*,
                                        const NN_scores*,
                                        const Alphabet*);

void
nn_scores_fprintf_tetra_loop(FILE*, const NN_scores*, const Alphabet*);

void
nn_scores_fprintf_G_dangle5(FILE*, const NN_scores*, const Alphabet*);

void
nn_scores_fprintf_G_dangle3(FILE*, const NN_scores*, const Alphabet*);

void
nn_scores_fprintf_G_mismatch_interior (FILE*,
                                       const NN_scores*,
                                       const Alphabet*);

void
nn_scores_fprintf_G_internal_loop (FILE*, const NN_scores*);

void
nn_scores_fprintf_G_int11 (FILE*, const NN_scores*, const Alphabet*);

void
nn_scores_fprintf_G_int21 (FILE*, const NN_scores*, const Alphabet*);

void
nn_scores_fprintf_G_int22 (FILE*, const NN_scores*, const Alphabet*);

/******************************   Miscellaneous   *****************************/

unsigned long
nn_scores_bp_2_idx (const char, const char, const NN_scores*);

bool
nn_scores_is_allowed_basepair (const char, const char, void*);

#endif /* NN_SCORES_H */

#ifdef __cplusplus
}
#endif
