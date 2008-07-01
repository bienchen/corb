/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
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
 *
 */


#include <config.h>
/* #include <stdio.h> */
#include <stdlib.h>
#include <math.h>
/* #include <assert.h> */
#include <libcrbbasic/crbbasic.h>
#include "alphabet.h"


struct Alphabet {
      char* upper_case;
      char* lower_case;
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
   /* allocate 1 object */
   Alphabet* this = XOBJ_MALLOC(sizeof (Alphabet), file, line);
   
   if (this != NULL)
   {
      this->upper_case = NULL;
      this->lower_case = NULL;
      this->size       = 0;
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
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
/*Alphabet*
alphabet_new_pair (const char* upper,
                   const char* lower,
                   unsigned long size,
                   const char* file, const int line)
{
  Alphabet* this = alphabet_new (file, line);

  if ((this != NULL) && (size != 0))
  {*/
     /* copy strings */
     /* set size */
/*}

  return this;
  }*/

/* alpha_new_single */

/********************************   Altering   ********************************/
/*********************************   Access   *********************************/
/******************* Size *******************/
/******************* Searching *******************/
/******************* Comparison *******************/

char
transform_base_2_number (const char base)
{
   /* RNA: "AaUuGgCc" HP: "LlRr" */
   char* alphabet="AaUuGgCc";   /* RrYyMmKkSsWwHhBbVvDdNn */
   size_t alpha_size = strlen (alphabet);
   size_t i;
   unsigned long shift = 0;

   for (i = 0; i < alpha_size; i++)
   {
      if (base == alphabet[i])
      {
         return i - shift;
      }
      if (i % 2 == 0)
      {
         shift++;
      }
   }
   
   THROW_ERROR_MSG ("Not a valid RNA nucleotide identifier: %c",
                    base);

   return CHAR_UNDEF;
}

char
transform_number_2_base (const char base)
{
   /* RNA: "AUGC" HP: "LR" */
   char* alphabet="AUGC";   /* RrYyMmKkSsWwHhBbVvDdNn */
   size_t alpha_size = strlen (alphabet);
   
   if ((size_t) base < alpha_size)
   {
      return alphabet[(size_t) base];
   }

   THROW_ERROR_MSG ("Not a valid RNA nucleotide identifier: %d",
                    (int)base);

   return CHAR_UNDEF;
}

float**
create_scoring_matrix (void)
{
   float** matrix = (float**) XMALLOC_2D (4, 4, sizeof (float));

   if (matrix == NULL)
   {
      return matrix;
   }

   matrix[(int) transform_base_2_number ('A')]
         [(int) transform_base_2_number ('A')] =  0.0f;
   matrix[(int) transform_base_2_number ('A')]
         [(int) transform_base_2_number ('U')] = -2.0f;
   matrix[(int) transform_base_2_number ('A')]
         [(int) transform_base_2_number ('G')] =  0.0f;
   matrix[(int) transform_base_2_number ('A')]
         [(int) transform_base_2_number ('C')] =  0.0f;
   matrix[(int) transform_base_2_number ('U')]
      [(int) transform_base_2_number ('A')] = -2.00f/*-1.99f*//*-1.99f*/;
   matrix[(int) transform_base_2_number ('U')]
         [(int) transform_base_2_number ('U')] =  0.0f;
   matrix[(int) transform_base_2_number ('U')]
         [(int) transform_base_2_number ('G')] = -1.5f;
   matrix[(int) transform_base_2_number ('U')]
         [(int) transform_base_2_number ('C')] =  0.0f;
   matrix[(int) transform_base_2_number ('G')]
         [(int) transform_base_2_number ('A')] =  0.0f;
   matrix[(int) transform_base_2_number ('G')]
      [(int) transform_base_2_number ('U')] = -1.50f/*-1.49f*//*-1.49f*/;
   matrix[(int) transform_base_2_number ('G')]
         [(int) transform_base_2_number ('G')] =  0.0f;
   matrix[(int) transform_base_2_number ('G')]
         [(int) transform_base_2_number ('C')] = -3.0f;
   matrix[(int) transform_base_2_number ('C')]
         [(int) transform_base_2_number ('A')] =  0.0f;
   matrix[(int) transform_base_2_number ('C')]
         [(int) transform_base_2_number ('U')] =  0.0f;
   matrix[(int) transform_base_2_number ('C')]
      [(int) transform_base_2_number ('G')] = -3.00f/*-2.99f*//*-2.99f*/;
   matrix[(int) transform_base_2_number ('C')]
         [(int) transform_base_2_number ('C')] =  0.0f;

   /* add random numbers to 3'-5' values */
   srand(30459);
   /* U - A */
   matrix[(int) transform_base_2_number ('U')]
         [(int) transform_base_2_number ('A')] = rand();
   matrix[(int) transform_base_2_number ('U')]
         [(int) transform_base_2_number ('A')] /=
       pow (10, (floor (log10 (matrix[(int) transform_base_2_number ('U')]
                                     [(int) transform_base_2_number ('A')]) 
                        + 1.0f)) + 2);
   matrix[(int) transform_base_2_number ('U')]
         [(int) transform_base_2_number ('A')] +=
                                  matrix[(int) transform_base_2_number ('A')]
                                  [(int) transform_base_2_number ('U')];

   /* G - U */
  matrix[(int) transform_base_2_number ('G')]
         [(int) transform_base_2_number ('U')] = rand();
   matrix[(int) transform_base_2_number ('G')]
         [(int) transform_base_2_number ('U')] /=
       pow (10, (floor (log10 (matrix[(int) transform_base_2_number ('G')]
                                     [(int) transform_base_2_number ('U')]) 
                        + 1.0f)) + 2);
   matrix[(int) transform_base_2_number ('G')]
         [(int) transform_base_2_number ('U')] +=
                                  matrix[(int) transform_base_2_number ('U')]
                                  [(int) transform_base_2_number ('G')];

   /* C - G */   
   matrix[(int) transform_base_2_number ('C')]
         [(int) transform_base_2_number ('G')] = rand();
   matrix[(int) transform_base_2_number ('C')]
         [(int) transform_base_2_number ('G')] /=
       pow (10, (floor (log10 (matrix[(int) transform_base_2_number ('C')]
                                     [(int) transform_base_2_number ('G')]) 
                        + 1.0f)) + 2);
   matrix[(int) transform_base_2_number ('C')]
         [(int) transform_base_2_number ('G')] +=
                                  matrix[(int) transform_base_2_number ('G')]
                                  [(int) transform_base_2_number ('C')];

   return matrix;
}
