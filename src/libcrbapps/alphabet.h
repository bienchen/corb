/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbapps/alphabet.h
 *
 *  @brief RNA alphabet
 *
 *  Module: alphabet
 *
 *  Library: libcrbapps
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-06-18
 *
 *
 *  Revision History:
 *         - 2008Jun18 bienert: created
 *
 */


#ifdef __cplusplus
extern "C" {
#endif

#ifndef ALPHABET_H
#define ALPHABET_H

char
transform_base_2_number (const char);

char
transform_number_2_base (const char);
 
float**
create_scoring_matrix (void);
  
#endif /* ALPHABET_H */

#ifdef __cplusplus
}
#endif
