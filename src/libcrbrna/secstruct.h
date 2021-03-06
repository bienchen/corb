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
 *  @file libcrbrna/secstruct.h
 *
 *  @brief Component container for RNA 2D structures
 *
 *  Module: secstruct
 *
 *  Library: crbrna
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-09-16
 *
 *
 *  Revision History:
 *         - 2008Sep16 bienert: created
 *
 */

#ifdef __cplusplus
extern "C" {
#endif


#ifndef SECSTRUCT_H
#define SECSTRUCT_H

/*#include "rna.h"*/
#include "nn_scores.h"

/* error numbers */
enum {
   ERR_SCSTRCT_ALLOC = 1,      /* (re)allocation problems */
   /* ERR_SCSTRCT_UKWN_FEAT, */      /* unknwon structural feature */
};

/* structural features */
typedef enum {
   SCSTRCT_VOID = 0,            /* base has no class assigned */
   SCSTRCT_HAIRPIN,             /* hairpin loop */
   SCSTRCT_STACK,               /* stacked basepair */
   SCSTRCT_BULGE,               /* bulge loop */
   SCSTRCT_INTERNAL,            /* internal loop */
   SCSTRCT_MULTI,               /* multi loop */
   SCSTRCT_EXTERNAL,            /* external loop */
   SCSTRCT_MTO,                 /* base is in M ore T han O ne class */
   N_FEATURES                   /* No. of elements */
} SecStructFtrs;

typedef struct SecStruct SecStruct;


/**********************   Constructors and destructors   **********************/
SecStruct*
secstruct_new (const char*, const int);

#define SECSTRUCT_NEW secstruct_new (__FILE__, __LINE__)

void
secstruct_delete (SecStruct*);

/********************************   Altering   ********************************/
int
secstruct_find_interactions (const unsigned long*,
                             const unsigned long,
                             SecStruct*);

void
secstruct_delete_stack (const unsigned long, SecStruct*);

void
secstruct_delete_element (const unsigned long,
                          const SecStructFtrs,
                          SecStruct*);

/*********************************   Access   *********************************/

SecStructFtrs
secstruct_get_feature_at_pos (const unsigned long,
                                unsigned long*,
                                const SecStruct*);

SecStructFtrs
secstruct_get_feature_multi_1st (const unsigned long,
                                   unsigned long*,
                                   const SecStruct*);

SecStructFtrs
secstruct_get_feature_multi_2nd (const unsigned long,
                                   unsigned long*,
                                   const SecStruct*);

/* hairpin loops */
unsigned long
secstruct_get_noof_hairpins (const SecStruct*);

void
secstruct_get_geometry_hairpin (unsigned long*, unsigned long*, unsigned long*,
                                const unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_start_hairpin (const unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_end_hairpin (const unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_size_hairpin (const unsigned long, const SecStruct*);

/* stacks */
unsigned long
secstruct_get_noof_stacks (const SecStruct*);

unsigned long
secstruct_get_i_5p_stack (const unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_3p_stack (const unsigned long, const SecStruct*);

void
secstruct_get_i_geometry_stack (unsigned long*, unsigned long*,
                                const unsigned long,
                                const SecStruct*);

/* bulge loops */
unsigned long
secstruct_get_noof_bulges (const SecStruct*);

unsigned long
secstruct_get_i_start_bulge (const unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_end_bulge (const unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_size_bulge (const unsigned long, const SecStruct*);

void
secstruct_get_geometry_bulge (unsigned long*, unsigned long*,
                              unsigned long*, unsigned long*,
                              unsigned long*,
                              const unsigned long,
                              const SecStruct*);

/* internal loops */
unsigned long
secstruct_get_noof_internals (const SecStruct*);

void
secstruct_get_geometry_internal (unsigned long*,
                                 unsigned long*,
                                 unsigned long*,
                                 unsigned long*,
                                 unsigned long*,
                                 unsigned long*,
                                 const unsigned long,
                                 const SecStruct*);

/* multiloops */
unsigned long
secstruct_get_noof_multiloops (const SecStruct*);

unsigned long
secstruct_get_i_noof_unpaired_multiloop (unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_noof_stems_multiloop (unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_5p_stem_multiloop (const unsigned long, const unsigned long,
                                   const SecStruct*);

unsigned long
secstruct_get_i_3p_stem_multiloop (const unsigned long, const unsigned long,
                                   const SecStruct*);

void
secstruct_get_i_stem_multiloop (unsigned long*, unsigned long*,
                                const unsigned long,
                                const unsigned long,
                                const SecStruct*);

unsigned long
secstruct_get_i_noof_5pdangles_multiloop (unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_5p_5pdangle_multiloop (const unsigned long, const unsigned long,
                                       const SecStruct*);

unsigned long
secstruct_get_i_3p_5pdangle_multiloop (const unsigned long, const unsigned long,
                                       const SecStruct*);

unsigned long
secstruct_get_i_dangle_5pdangle_multiloop (const unsigned long,
                                           const unsigned long,
                                           const SecStruct*);

void
secstruct_get_i_5pdangle_multiloop (unsigned long*,
                                    unsigned long*,
                                    unsigned long*,
                                    const unsigned long,
                                    const unsigned long,
                                    const SecStruct*);

unsigned long
secstruct_get_i_noof_3pdangles_multiloop (unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_5p_3pdangle_multiloop (const unsigned long, const unsigned long,
                                       const SecStruct*);

unsigned long
secstruct_get_i_3p_3pdangle_multiloop (const unsigned long, const unsigned long,
                                       const SecStruct*);

unsigned long
secstruct_get_i_dangle_3pdangle_multiloop (const unsigned long,
                                           const unsigned long,
                                           const SecStruct*);

void
secstruct_get_i_3pdangle_multiloop (unsigned long*,
                                    unsigned long*,
                                    unsigned long*,
                                    const unsigned long,
                                    const unsigned long,
                                    const SecStruct*);

/* external loops */
unsigned long
secstruct_get_i_noof_unpaired_extloop (const SecStruct*);

unsigned long
secstruct_get_noof_stems_extloop (const SecStruct*);

unsigned long
secstruct_get_i_5p_stem_extloop (const unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_3p_stem_extloop (const unsigned long, const SecStruct*);

void
secstruct_get_i_stem_extloop (unsigned long*, unsigned long*,
                              const unsigned long, const SecStruct*);

unsigned long
secstruct_get_noof_5pdangles_extloop (const SecStruct*);

unsigned long
secstruct_get_noof_3pdangles_extloop (const SecStruct*);

unsigned long
secstruct_get_i_5p_3pdangle_extloop (const unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_3p_3pdangle_extloop (const unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_5p_5pdangle_extloop (const unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_3p_5pdangle_extloop (const unsigned long, const SecStruct*);

void
secstruct_get_i_5pdangle_extloop (unsigned long*,
                                  unsigned long*,
                                  unsigned long*,
                                  const unsigned long,
                                  const SecStruct*);

void
secstruct_get_i_3pdangle_extloop (unsigned long*,
                                  unsigned long*,
                                  unsigned long*,
                                  const unsigned long,
                                  const SecStruct*);

int
secstruct_calculate_DG (const char*, const NN_scores*, const SecStruct*);

/*********************************   Output   *********************************/

void
secstruct_fprintf_i_stack (FILE*, const unsigned long, const SecStruct*);

void
secstruct_fprintf_stacks (FILE*, const SecStruct*);

void
secstruct_fprintf_i_hairpin (FILE*, const unsigned long, const SecStruct*);

void
secstruct_fprintf_hairpins (FILE*, const SecStruct*);

void
secstruct_fprintf_i_bulge (FILE*, const unsigned long, const SecStruct*);

void
secstruct_fprintf_bulges (FILE*, const SecStruct*);

void
secstruct_fprintf_i_internal (FILE*, const unsigned long, const SecStruct*);

void
secstruct_fprintf_internals (FILE*, const SecStruct*);

void
secstruct_fprintf_external (FILE*, const SecStruct*);

void
secstruct_fprintf_i_multiloop (FILE*, const unsigned long, const SecStruct*);

void
secstruct_fprintf_multiloops (FILE*, const SecStruct*);

void
secstruct_fprintf_seqpos_map (FILE*, const SecStruct*);

#endif /* SECSTRUCT_H */

#ifdef __cplusplus
}
#endif
