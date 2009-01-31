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
 *  @file libcrbrna/rna.h
 *
 *  @brief RNA datastructure
 *
 *  Module: rna
 *
 *  Library: crbrna
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-08-19
 *
 *
 *  Revision History:
 *         - 2008Aug19 bienert: created
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <libcrbbasic/crbbasic.h>
#include "nn_scores.h"
#include "alphabet.h"

#ifndef RNA_H
#define RNA_H


/* unpaired nucleotide */
#define NOT_PAIRED ULONG_UNDEF
#define VIENNA_OPEN  '('
#define VIENNA_CLOSE ')'
#define VIENNA_UNPAIRED '.'

/* error numbers */
enum rna_retvals{
   ERR_RNA_ALLOC = 1,      /* (re)allocation problems */
   ERR_RNA_VIENNA_FORMAT,  /* Wrong format of Vienna structure string */
   ERR_RNA_VIENNA_MMC,     /* Mismatched closing pairing partner found */
   ERR_RNA_VIENNA_MMO,     /* Mismatched opening pairing partner found */
   ERR_RNA_NO_BASE,        /* Given symbol is not a valid base */
   ERR_RNA_FILE_OPEN,      /* Problems on opening a file */
   ERR_RNA_CT_NAN,         /* Ct column does not contain a number */
   ERR_RNA_CT_NN,          /* Ct column does not contain a single letter base */
   ERR_RNA_CT_SM,          /* Ct line semantic problem */
   ERR_RNA_EXT_NR,         /* File extension not recognised*/
};

typedef struct Rna Rna;


/**********************   Constructors and destructors   **********************/

Rna*
rna_new (const char*, const int);

#define RNA_NEW rna_new (__FILE__, __LINE__)

void
rna_delete (Rna*);
   
int
rna_allocate_pairlist (const unsigned long, Rna*, const char*, const int);

#define RNA_ALLOCATE_PAIRLIST(A, B) \
   rna_allocate_pairlist (A, B, __FILE__, __LINE__)

int
rna_realloc_pairlist (const unsigned long, Rna*, const char*, const int);

#define RNA_REALLOC_PAIRLIST(A, B) \
   rna_realloc_pairlist (A, B, __FILE__, __LINE__)

int
rna_init_pairlist_vienna (const char*, const unsigned long, Rna*,
                          const char*, const int);

#define RNA_INIT_PAIRLIST_VIENNA(A, B, C) \
   rna_init_pairlist_vienna (A, B, C, __FILE__, __LINE__);

int
rna_alloc_sequence (const unsigned long, Rna*, const char*, const int);

int
rna_realloc_sequence (const unsigned long, Rna*, const char*, const int);

int
rna_init_sequence (const char*, const unsigned long, Alphabet*, Rna*,
                   const char*, const int);

#define RNA_INIT_SEQUENCE(SEQ, LENGTH, SIGMA, RNA)  \
   rna_init_sequence (SEQ, LENGTH, SIGMA, RNA, __FILE__, __LINE__)

int
rna_init_sequence_structure (const char*, const char*, const unsigned long,
                             Alphabet*, Rna*,
                             const char*, const int);

#define RNA_INIT_SEQUENCE_STRUCTURE(A, B, C, D, E) \
   rna_init_sequence_structure (A, B, C, D, E, __FILE__, __LINE__)


#define RNA_ALLOC_SEQUENCE(A, B) \
   rna_alloc_sequence (A, B, __FILE__, __LINE__)

#define RNA_REALLOC_SEQUENCE (A, B) \
   rna_realloc_sequence (A, B, __FILE__, __LINE__)

int
rna_read_from_file_ct (Rna*, GFile*);

int
rna_read_from_file (Rna*, const char*, const unsigned long, 
                    const char*, const unsigned long, const GFileType);

int
rna_secstruct_init (Rna*, const char*, const int);

#define RNA_SECSTRUCT_INIT(A) \
   rna_secstruct_init (A, __FILE__, __LINE__)

/********************************   Altering   ********************************/

void
rna_set_sequence_base (const char, const unsigned long, Rna*);

void
rna_set_pairing_base (const unsigned long, const unsigned long, Rna*);

void
rna_set_sequence (const char*, const unsigned long, Rna*);

int
rna_transform_sequence_2_no (const Alphabet*, Rna*);

int
rna_transform_sequence_2_bases (const Alphabet*, Rna*);

/*********************************   Access   *********************************/

unsigned long
rna_get_size (const Rna*);

const unsigned long*
rna_get_pairlist (const Rna*);

char*
rna_get_sequence (const Rna*);

const Str*
rna_get_info (const Rna*);

unsigned long
rna_base_pairs_with (const unsigned long, const Rna*);

char
rna_get_sequence_base (const unsigned long, const Rna*);

unsigned long
rna_validate_basepairs (bool (*validate_basepair) (const char, const char,
                                                   void*), void*, const Rna*);

unsigned long
rna_validate_basepairs_nn_scores (NN_scores*, const Rna*);

int
rna_secstruct_calculate_DG (const NN_scores*, const Rna*);

unsigned long
rna_secstruct_get_noof_stacks (const Rna*);

void
rna_secstruct_get_i_geometry_stack (unsigned long*, unsigned long*,
                                    const unsigned long,
                                    const Rna*);

unsigned long
rna_secstruct_get_noof_hairpins (const Rna*);

void
rna_secstruct_get_geometry_hairpin (unsigned long*,
                                    unsigned long*,
                                    unsigned long*,
                                    const unsigned long,
                                    const Rna*);

unsigned long
rna_secstruct_get_i_start_hairpin (const unsigned long, const Rna*);

unsigned long
rna_secstruct_get_i_end_hairpin (const unsigned long, const Rna*);

unsigned long
rna_secstruct_get_i_size_hairpin (const unsigned long, const Rna*);

unsigned long
rna_secstruct_get_noof_stems_extloop (const Rna*);

unsigned long
rna_secstruct_get_i_5p_stem_extloop (const unsigned long, const Rna*);

unsigned long
rna_secstruct_get_i_3p_stem_extloop (const unsigned long, const Rna*);

void
rna_secstruct_get_i_stem_extloop (unsigned long*, unsigned long*,
                                  const unsigned long, const Rna*);

unsigned long
rna_secstruct_get_noof_5pdangles_extloop (const Rna*);

unsigned long
rna_secstruct_get_noof_3pdangles_extloop (const Rna*);

void
rna_secstruct_get_i_5pdangle_extloop (unsigned long*,
                                      unsigned long*,
                                      unsigned long*,
                                      const unsigned long,
                                      const Rna*);

void
rna_secstruct_get_i_3pdangle_extloop (unsigned long*,
                                      unsigned long*,
                                      unsigned long*,
                                      const unsigned long,
                                      const Rna*);

unsigned long
rna_secstruct_get_noof_multiloops (const Rna*);

unsigned long
rna_secstruct_get_i_noof_stems_multiloop (unsigned long, const Rna*);

void
rna_secstruct_get_i_stem_multiloop (unsigned long*, unsigned long*,
                                    const unsigned long,
                                    const unsigned long,
                                    const Rna*);

unsigned long
rna_secstruct_get_i_noof_5pdangles_multiloop (unsigned long, const Rna*);

void
rna_secstruct_get_i_5pdangle_multiloop (unsigned long*,
                                        unsigned long*,
                                        unsigned long*,
                                        const unsigned long,
                                        const unsigned long,
                                        const Rna*);

unsigned long
rna_secstruct_get_i_noof_3pdangles_multiloop (unsigned long, const Rna*);

void
rna_secstruct_get_i_3pdangle_multiloop (unsigned long*,
                                        unsigned long*,
                                        unsigned long*,
                                        const unsigned long,
                                        const unsigned long,
                                        const Rna*);

unsigned long
rna_secstruct_get_noof_bulges (const Rna*);

void
rna_secstruct_get_geometry_bulge (unsigned long*, unsigned long*,
                                  unsigned long*, unsigned long*,
                                  unsigned long*,
                                  const unsigned long,
                                  const Rna*);

unsigned long
rna_secstruct_get_noof_internals (const Rna*);

void
rna_secstruct_get_geometry_internal (unsigned long*,
                                     unsigned long*,
                                     unsigned long*,
                                     unsigned long*,
                                     unsigned long*,
                                     unsigned long*,
                                     const unsigned long,
                                     const Rna*);

#endif /* RNA_H */

#ifdef __cplusplus
}
#endif
