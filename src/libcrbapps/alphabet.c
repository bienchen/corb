/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbapps/alphabet.c
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


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <libcrbbasic/crbbasic.h>
#include "alphabet.h"


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
