/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
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

#include "rna.h"

#ifndef SECSTRUCT_H
#define SECSTRUCT_H

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


/*********************************   Access   *********************************/

/* hairpin loops */
unsigned long
secstruct_get_noof_hairpins (const SecStruct*);

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

/* bulge loops */
unsigned long
secstruct_get_noof_bulges (const SecStruct*);

unsigned long
secstruct_get_i_start_bulge (const unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_end_bulge (const unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_size_bulge (const unsigned long, const SecStruct*);

/* internal loops */
unsigned long
secstruct_get_noof_internals (const SecStruct*);

unsigned long
secstruct_get_i_start_internal (const unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_end_internal (const unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_size_internal (const unsigned long, const SecStruct*);

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

unsigned long
secstruct_get_i_noof_5pdangles_multiloop (unsigned long,
                                          const SecStruct*);

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

/* external loops */
unsigned long
secstruct_get_i_noof_unpaired_extloop (const SecStruct*);

unsigned long
secstruct_get_i_noof_stems_extloop (const SecStruct*);

unsigned long
secstruct_get_i_5p_stem_extloop (const unsigned long, const SecStruct*);

unsigned long
secstruct_get_i_3p_stem_extloop (const unsigned long, const SecStruct*);

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

#endif /* SECSTRUCT_H */

#ifdef __cplusplus
}
#endif
