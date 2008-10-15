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
#define RNA_ALPHABET "ACGUacgu"
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
