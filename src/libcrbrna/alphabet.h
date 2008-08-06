/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbrna/alphabet.h
 *
 *  @brief RNA alphabet
 *
 *  Module: alphabet
 *
 *  Library: libcrbrna
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-06-30
 *
 *
 *  Revision History:
 *         - 2008Jun30 bienert: created
 *
 */


#ifdef __cplusplus
extern "C" {
#endif

#ifndef ALPHABET_H
#define ALPHABET_H


/* alphabet strings */
#define RNA_ALPHABET "AUGCaugc"
#define IUPAC_NA_STRICT_ALPHABET "ACGTURYMKWSBDHVN"
#define IUPAC_NA_ALPHABET "ACGTURYMKWSBDHVNacgturymkwsbdhvn"


typedef struct Alphabet Alphabet;


/**********************   Constructors and destructors   **********************/

Alphabet*
alphabet_new (const char*, const int);

#define ALPHABET_NEW alphabet_new (__FILE__, __LINE__)

Alphabet*
alphabet_new_pair (const char*, const char*, const unsigned long,
                   const char*, const int);

#define ALPHABET_NEW_PAIR(U, L, S) \
        alphabet_new_pair (U, L, S, __FILE__, __LINE__)

Alphabet*
alphabet_new_single (const char*, const unsigned long, const char*, const int);

#define ALPHABET_NEW_SINGLE(A,S) \
        alphabet_new_single (A, S, __FILE__, __LINE__)

void
alphabet_delete (Alphabet*);

/********************************   Altering   ********************************/
/*********************************   Access   *********************************/
/******************* Size *******************/

unsigned long
alphabet_size (const Alphabet*);

/******************* Searching *******************/
/******************* Comparison *******************/

bool
alphabet_is_standard_rna (Alphabet*);

char
alphabet_base_2_no (const char, const Alphabet*);

char
alphabet_no_2_base (const char, const Alphabet*);

float**
create_scoring_matrix (const Alphabet*);
  
#endif /* ALPHABET_H */

#ifdef __cplusplus
}
#endif
