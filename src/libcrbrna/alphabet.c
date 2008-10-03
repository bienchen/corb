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
 *  @file libcrbrna/alphabet.c
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
 */


#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <libcrbbasic/crbbasic.h>
#include "alphabet.h"


struct Alphabet {
      char* upper_case;
      char* lower_case;
      char idx[CHAR_MAX];
      unsigned long size;
};


/**********************   Constructors and destructors   **********************/

/** @brief Create a new alphabet object.
 *
 * The constructor for @c Alphabet objects. If compiled with enabled
 * memory checking, @c file and @c line should point to the position where the
 * function was called. Both parameters are automatically set by using the
 * macro @c ALPHABET_NEW.\n
 * Returns @c NULL on error.
 *
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
Alphabet*
alphabet_new (const char* file, const int line)
{
   unsigned char i;

   /* allocate 1 object */
   Alphabet* this = XOBJ_MALLOC(sizeof (Alphabet), file, line);
   
   if (this != NULL)
   {
      this->upper_case = NULL;
      this->lower_case = NULL;
      this->size       = 0;

      for (i = 0; i < CHAR_MAX; i++)
      {
         this->idx[i] = CHAR_UNDEF;
      }
   }

   return this;
}

/** @brief Create a new alphabet object from two c strings.
 *
 * A constructor for an initialised @c Alphabet object. If compiled with enabled
 * memory checking, @c file and @c line should point to the position where the
 * function was called. Both parameters are automatically set by using the
 * macro @c ALPHABET_NEW_PAIR. @c upper and @lower have to be of same size.
 * Symbols at equal positions are taken as upper and lower case pendants of the
 * same symbol in the alphabet.\n
 * Returns @c NULL on error.
 *
 * @param[in] upper alphabet in upper case letters.
 * @param[in] lower alphabet in lower case letters.
 * @param[in] size no. of elements of the alphabet.
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
Alphabet*
alphabet_new_pair (const char* upper,
                   const char* lower,
                   const unsigned long size,
                   const char* file, const int line)
{
   unsigned long i;
   Alphabet* this = alphabet_new (file, line);

   if ((this != NULL) && (size != 0))
   {
      /* set size */
      this->size = size;

      /* copy strings */
      this->upper_case = XMALLOC (this->size * sizeof (this->upper_case[0]));
      if (this->upper_case == NULL)
      {
         alphabet_delete (this);
         return NULL;
      }
      this->lower_case = XMALLOC (this->size * sizeof (this->upper_case[0]));   
      if (this->upper_case == NULL)
      {
         alphabet_delete (this);
         return NULL;
      }

      for (i = 0; i < this->size; i++)
      {
         this->upper_case[i] = upper[i];
         this->lower_case[i] = lower[i];
         this->idx[(int) this->upper_case[i]] = i;
         this->idx[(int) this->lower_case[i]] = i;
      }
   }
   
   return this;
}

/** @brief Create a new alphabet object from a single c string.
 *
 * A constructor for an initialised @c Alphabet object. If compiled with enabled
 * memory checking, @c file and @c line should point to the position where the
 * function was called. Both parameters are automatically set by using the
 * macro @c ALPHABET_NEW_SINGLE. @c alphabet must contain an even number of
 * elements. The string has to consist of corresponding pairs (e.g. ACGUacgu)\n
 * Returns @c NULL on error.
 *
 * @param[in] alphabet alphabet string.
 * @param[in] size no. of elements of the alphabet.
 * @param[in] format format of the alphabet, either 
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
Alphabet*
alphabet_new_single (const char* alphabet,
                     const unsigned long size,
                     const char* file, const int line)
{
   Alphabet* this;

   this = alphabet_new_pair (alphabet,
                             alphabet + size,
                             size,
                             file, line);

   return this;
}

/** @brief Delete an alphabet.
 *
 * The destructor for @c Alphabet objects.
 *
 * @param[in] this object to be freed.
 */
void
alphabet_delete (Alphabet* this)
{
   if (this != NULL)
   {
     XFREE (this->upper_case);
     XFREE (this->lower_case);
     XFREE (this);
   }
}


/********************************   Altering   ********************************/
/*********************************   Access   *********************************/
/******************* Size *******************/

unsigned long
alphabet_size (const Alphabet* this)
{
   assert (this != NULL);

   return this->size;
}


/******************* Searching *******************/
/******************* Comparison *******************/

/** @brief Check if an alphabet is the standard 'AUGC' RNA alphabet.
 *
 * Checks if an alphabet contains only the standrad RNA nucleotides 'AUGC'.\n
 * Returns @c true or @c false.
 *
 * @param[in] sigma alphabet to be checked.
 */
bool
alphabet_is_standard_rna (Alphabet* sigma)
{
   int i;
   char* rna = RNA_ALPHABET;
   int a;

   assert (sigma != NULL);
   assert (sigma->idx != NULL);
   assert (sigma->upper_case != NULL);
   assert (sigma->lower_case != NULL);

   /* simplest check: Size of alphabet */
   if (sigma->size != (strlen (RNA_ALPHABET) / 2))
   {
      THROW_ERROR_MSG ("Alphabet contains non-standard RNA nucleotides: %s%s",
                       sigma->upper_case, sigma->lower_case);
      return false;
   }

   /* check if an index exists for each nucleotide*/
   for (i = 0; (unsigned) i < (sigma->size * 2); i++)
   {
      a = rna[i];
      if (sigma->idx[a] == CHAR_UNDEF)
      {
         THROW_ERROR_MSG ("Alphabet does not contain standard RNA nucleotide "
                          "'%c': %s%s", rna[i], sigma->upper_case,
                          sigma->lower_case);
         return false;         
      }
   }

   return true;
}



char
alphabet_base_2_no (const char base,
                    const Alphabet* sigma)
{
   /* RNA: "AaUuGgCc" HP: "LlRr" */
   char no;

   assert (sigma != NULL);
   assert (sigma->idx != NULL);
   
   no = sigma->idx[(int) base];

   if (no == CHAR_UNDEF)
   {
      THROW_ERROR_MSG ("Not a valid RNA nucleotide identifier: %c", base);
   }

   return no;
}

char
alphabet_no_2_base (const char base,
                    const Alphabet* sigma __attribute__((unused)))
{
   assert (sigma != NULL);
   assert (sigma->upper_case != NULL);

   if (sigma->size > (unsigned) base)
   {
      return sigma->upper_case[(int) base];
   }

   THROW_ERROR_MSG ("Not a valid RNA nucleotide identifier: %d",
                    (int)base);

   return CHAR_UNDEF;
}

float**
create_scoring_matrix (const Alphabet* sigma)
{
   float** matrix = (float**) XMALLOC_2D (4, 4, sizeof (float));

   if (matrix == NULL)
   {
      return matrix;
   }

   matrix[(int) alphabet_base_2_no ('A', sigma)]
         [(int) alphabet_base_2_no ('A', sigma)] =  0.0f;
   matrix[(int) alphabet_base_2_no ('A', sigma)]
         [(int) alphabet_base_2_no ('U', sigma)] = -2.0f;
   matrix[(int) alphabet_base_2_no ('A', sigma)]
         [(int) alphabet_base_2_no ('G', sigma)] =  0.0f;
   matrix[(int) alphabet_base_2_no ('A', sigma)]
         [(int) alphabet_base_2_no ('C', sigma)] =  0.0f;
   matrix[(int) alphabet_base_2_no ('U', sigma)]
      [(int) alphabet_base_2_no ('A', sigma)] = -2.00f/*-1.99f*//*-1.99f*/;
   matrix[(int) alphabet_base_2_no ('U', sigma)]
         [(int) alphabet_base_2_no ('U', sigma)] =  0.0f;
   matrix[(int) alphabet_base_2_no ('U', sigma)]
         [(int) alphabet_base_2_no ('G', sigma)] = -1.5f;
   matrix[(int) alphabet_base_2_no ('U', sigma)]
         [(int) alphabet_base_2_no ('C', sigma)] =  0.0f;
   matrix[(int) alphabet_base_2_no ('G', sigma)]
         [(int) alphabet_base_2_no ('A', sigma)] =  0.0f;
   matrix[(int) alphabet_base_2_no ('G', sigma)]
      [(int) alphabet_base_2_no ('U', sigma)] = -1.50f/*-1.49f*//*-1.49f*/;
   matrix[(int) alphabet_base_2_no ('G', sigma)]
         [(int) alphabet_base_2_no ('G', sigma)] =  0.0f;
   matrix[(int) alphabet_base_2_no ('G', sigma)]
         [(int) alphabet_base_2_no ('C', sigma)] = -3.0f;
   matrix[(int) alphabet_base_2_no ('C', sigma)]
         [(int) alphabet_base_2_no ('A', sigma)] =  0.0f;
   matrix[(int) alphabet_base_2_no ('C', sigma)]
         [(int) alphabet_base_2_no ('U', sigma)] =  0.0f;
   matrix[(int) alphabet_base_2_no ('C', sigma)]
      [(int) alphabet_base_2_no ('G', sigma)] = -3.00f/*-2.99f*//*-2.99f*/;
   matrix[(int) alphabet_base_2_no ('C', sigma)]
         [(int) alphabet_base_2_no ('C', sigma)] =  0.0f;

   /* add random nos to 3'-5' values */
   srand(30459);
   /* U - A */
   matrix[(int) alphabet_base_2_no ('U', sigma)]
         [(int) alphabet_base_2_no ('A', sigma)] = rand();
   matrix[(int) alphabet_base_2_no ('U', sigma)]
         [(int) alphabet_base_2_no ('A', sigma)] /=
       pow (10, (floor (log10 (matrix[(int) alphabet_base_2_no ('U', sigma)]
                                     [(int) alphabet_base_2_no ('A', sigma)]) 
                        + 1.0f)) + 2);
   matrix[(int) alphabet_base_2_no ('U', sigma)]
         [(int) alphabet_base_2_no ('A', sigma)] +=
                                  matrix[(int) alphabet_base_2_no ('A', sigma)]
                                  [(int) alphabet_base_2_no ('U', sigma)];

   /* G - U */
  matrix[(int) alphabet_base_2_no ('G', sigma)]
         [(int) alphabet_base_2_no ('U', sigma)] = rand();
   matrix[(int) alphabet_base_2_no ('G', sigma)]
         [(int) alphabet_base_2_no ('U', sigma)] /=
       pow (10, (floor (log10 (matrix[(int) alphabet_base_2_no ('G', sigma)]
                                     [(int) alphabet_base_2_no ('U', sigma)]) 
                        + 1.0f)) + 2);
   matrix[(int) alphabet_base_2_no ('G', sigma)]
         [(int) alphabet_base_2_no ('U', sigma)] +=
                                  matrix[(int) alphabet_base_2_no ('U', sigma)]
                                  [(int) alphabet_base_2_no ('G', sigma)];

   /* C - G */   
   matrix[(int) alphabet_base_2_no ('C', sigma)]
         [(int) alphabet_base_2_no ('G', sigma)] = rand();
   matrix[(int) alphabet_base_2_no ('C', sigma)]
         [(int) alphabet_base_2_no ('G', sigma)] /=
       pow (10, (floor (log10 (matrix[(int) alphabet_base_2_no ('C', sigma)]
                                     [(int) alphabet_base_2_no ('G', sigma)]) 
                        + 1.0f)) + 2);
   matrix[(int) alphabet_base_2_no ('C', sigma)]
         [(int) alphabet_base_2_no ('G', sigma)] +=
                                  matrix[(int) alphabet_base_2_no ('G', sigma)]
                                  [(int) alphabet_base_2_no ('C', sigma)];

   return matrix;
}
