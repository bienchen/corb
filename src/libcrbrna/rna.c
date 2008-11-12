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
 *  @file libcrbrna/rna.c
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


#include <config.h>
#include <limits.h>
#include <stddef.h>
#include <libcrbbasic/crbbasic.h>
/*#include "alphabet.h"*/
#include "secstruct.h"
/*#include "nn_scores.h"*/
#include "rna.h"


typedef struct {
      ArrayUlong tetra_loop;
} struct_comp;


struct Rna {
      char*          seq;       /* the nucleotide sequence */
      /* char* vienna; */       /* vienna string */
      unsigned long* pairs;     /* base pairs */
      unsigned long  size;      /* size of the RNA (sequence & 2D structure) */
       SecStruct*    structure; /* decomposed structure */
};


/**********************   Constructors and destructors   **********************/

/** @brief Init a struct_comp structure.
 *
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
/* static __inline__ struct_comp */
/* struct_comp_new (void/\* const char* file, const int line *\/) */
/* { */
/*    struct_comp this; */

/*    ARRAY_ULONG_INIT(this.tetra_loop, 0); */

/*    return this; */
/* } */

/** @brief Delete the components of a struct_comp structure.
 *
 * @param[in] this Object to be freed.
 */
/* static __inline__ void */
/* struct_comp_delete (struct_comp* this) */
/* { */
/*    if (this != NULL) */
/*    { */
/*      ARRAY_DELETE(this->tetra_loop); */
/*    } */
/* } */

/** @brief Create a new rna object.
 *
 * The constructor for @c Rna objects. If compiled with enabled
 * memory checking, @c file and @c line should point to the position where the
 * function was called. Both parameters are automatically set by using the
 * macro @c RNA_NEW.\n
 * Returns @c NULL on error.
 *
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
Rna*
rna_new (const char* file, const int line)
{
   Rna* this = XOBJ_MALLOC(sizeof (*this), file, line);
   
   if (this != NULL)
   {
      this->size      = 0;
      this->seq       = NULL;
      /* this->vienna    = NULL; */
      this->pairs     = NULL;
      this->structure = NULL;
   }

   return this;
}

/** @brief Delete a Rna object.
 *
 * The destructor for @c Rna objects.
 *
 * @param[in] this Object to be freed.
 */
void
rna_delete (Rna* this)
{

   if (this != NULL)
   {
     XFREE(this->seq);
     /* XFREE(this->vienna); */
     XFREE(this->pairs);
     secstruct_delete (this->structure);
     XFREE(this);
   }
}

/** @brief Allocate memory for the sequence component of an @c Rna object.
 *
 * Allocate the memory for a sequence stored within an @c Rna data object. If
 * compiled with enabled memory checking, @c file and @c line should point to
 * the position where the function was called. Both parameters are 
 * automatically set by using the macro @c RNA_ALLOCATE_SEQUENCE.\n
 * Returns 0 on success, @c ERR_RNA_ALLOC on memory problems.
 *
 * @param[in] size Size of the sequence/ structure of the RNA.
 * @param[in] this Rna data object.
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
int
rna_alloc_sequence (const unsigned long size, Rna* this,
                    const char* file, const int line)
{
   assert (this);
   assert (this->seq == NULL);
   assert ((this->size == 0) || (this->size == size));

   this->size = size;

   this->seq = XOBJ_CALLOC ((this->size + 1), sizeof (this->seq[0]),
                            file, line);
   if (this->seq == NULL)
   {
      return ERR_RNA_ALLOC;
   }

   return 0;
}

/** @brief Store a copy of an RNA sequence in an Rna data object.
 *
 * Initialise the sequence component of an Rna object. The data source is a
 * character string over a certain alphabet. The sequence is transforme into
 * the internal RNA representation and therefore stored as a copy in the
 * object. If compiled with enabled memory checking, @c file and @c line
 * should point to the position where the function was called. Both parameters
 * are automatically set by using the macro @c RNA_INIT_SEQUENCE.\n
 * Returns 0 on success, @c ERR_RNA_ALLOC on memory problems. On problems
 * parsing the sequence string, the following error codes are returned:\n
 * @list
 * @c ERR_RNA_NO_BASE A base in the sequence is not in the given alphabet.
 *
 * @params[in] seq Sequence string.
 * @params[in] length Length of the sequence.
 * @params[in] sigma Alphabet.
 * @params[in] this Rna object to store structure.
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
int
rna_init_sequence (const char* seq,
                   const unsigned long length,
                   Alphabet* sigma,
                   Rna* this,
                   const char* file, const int line)
{
   int error = 0;
   unsigned long i;

   assert (this);
   assert (this->seq == NULL);
   assert ((this->size == 0) || (this->size == length));

   /* allocate sequence */
   error = rna_alloc_sequence (length, this, file, line);

   /* store & validate sequence */
   i = 0;
   while ((!error) && (i < length))
   {
      this->seq[i] = alphabet_base_2_no (seq[i], sigma);

      if (this->seq[i] == CHAR_UNDEF)
      {
         error =  ERR_RNA_NO_BASE;
      }
      i++;
   }

   return error;
}

/** @brief Allocate memory for the pair list component of an @c Rna object.
 *
 * Allocate the memory for a list of pairs stored within an @c Rna data object.
 * If compiled with enabled memory checking, @c file and @c line should point to
 * the position where the function was called. Both parameters are 
 * automatically set by using the macro @c RNA_ALLOCATE_PAIRLIST.\n
 * Returns 0 on success, @c ERR_RNA_ALLOC on memory problems.
 *
 * @param[in] size Size of the sequence/ structure of the RNA.
 * @param[in] this Rna data object.
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
int
rna_allocate_pairlist (const unsigned long size, Rna* this,
                       const char* file, const int line)
{
   assert (this);
   assert (this->pairs == NULL);
   assert ((this->size == 0) || (this->size == size));  
   
   this->pairs = XOBJ_MALLOC(sizeof (this->pairs[0]) * size, file, line);
   memset (this->pairs, INT_MAX, sizeof (this->pairs[0]) * size);   

   if (this->pairs == NULL)
   {
      return ERR_RNA_ALLOC;
   }

   return 0;
}

/** @brief Read the list of base pairs from a Vienna string.
 *
 * Initialise the base pairs component of an Rna object. The data source is a
 * structure string in Vienna notation which defines size of the and pairs in
 * the list. Position counting starts by 1. Non-paired positions will get a
 * value of @c NOT_PAIRED. Additionally the Vienna string itself is copied to
 * the data object. If compiled with enabled memory checking, @c file and
 * @c line should point to the position where the function was called. Both
 * parameters are automatically set by using the macro
 * @c RNA_INIT_PAIRLIST_VIENNA.\n
 * Returns 0 on success, @c ERR_RNA_ALLOC on memory problems. For problems
 * parsing the Vienna string, the following error codes are resrved:\n
 * @list
 * @c ERR_RNA_VIENNA_FORMAT if the string contains unknown symbols
 * @c ERR_RNA_VIENNA_MMC    if there are more closing then opening paring
 *                          partners
 * @c ERR_RNA_VIENNA_MMO    if there are more opening then closing paring
 *                          partners
 *
 * @params[in] vienna Structure string.
 * @params[in] length Length of @c vienna.
 * @params[in] this Rna object to store structure.
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
int
rna_init_pairlist_vienna (const char* vienna,
                          const unsigned long length,
                          Rna* this,
                          const char* file, const int line)
{
   int error = 0;
   unsigned long i;
   unsigned long open_pair;
   unsigned long* p_stack = NULL;

   assert (this);
   /* assert (this->vienna == NULL); */
   assert (this->pairs == NULL);
   assert ((this->size == 0) || (this->size == length));

   /* allocate list */
   error = rna_allocate_pairlist (length, this, file, line);

   if (!error)
   {
      this->size = length;
      p_stack = XCALLOC (this->size, sizeof (p_stack[0]));
      if (p_stack == NULL)
      {
         error = ERR_RNA_ALLOC;
      }
   }

   /* transform structure string to base pair list */
   i = 0;
   open_pair = 0;
   while ((!error) && (i < this->size))
   {
      switch (vienna[i])
      {
         case VIENNA_OPEN :
            p_stack[open_pair] = i;
            open_pair++;
            break;     
            
         case VIENNA_CLOSE :
            if (open_pair == 0)
            {
               THROW_ERROR_MSG ("Mismatched nucleotide (closing base pair "
                                "partner) found at position %lu in structure "
                                "string (\'%s\')", i, vienna);
               error = ERR_RNA_VIENNA_MMC;               
            }
            else
            {
               open_pair--;
               this->pairs[i] = p_stack[open_pair];
               this->pairs[p_stack[open_pair]] = i;
            }
            break;
            
         case VIENNA_UNPAIRED : ;
            break;     
            
         default  : 
            THROW_ERROR_MSG ("Non-valid symbol found in structure string "
                             "(\'%s\'). Allowed characters: \'%c\', \'%c\', "
                             "\'%c\'. Found \'%c\' at position %lu.",
                             vienna,
                             VIENNA_OPEN,
                             VIENNA_CLOSE,
                             VIENNA_UNPAIRED,
                             vienna[i], i);
            error = ERR_RNA_VIENNA_FORMAT;
            break;
      }

      i++;
   }

   if (open_pair != 0)
   {
      THROW_ERROR_MSG ("Mismatched nucleotide (opening base pair "
                       "partner) found in structure string (\'%s\')", vienna);
      error = ERR_RNA_VIENNA_MMO;
   }
   
   XFREE (p_stack);

   if (error)
   {
      XFREE (this->pairs);
      this->pairs = NULL;
   }

   return error;
}

/** @brief Initialise an Rna data object with a sequence and a structure.
 *
 * Store sequence and structure in an Rna data object. Does the same as using
 * @c rna_init_sequence together with @c rna_init_pairlist_vienna. Since there
 * is only one @c length parameter, we assume that sequence and structure are
 * of equal length. If compiled with enabled memory checking, @c file and
 * @c line should point to the position where the function was called. Both
 * parameters are automatically set by using the macro
 * @c RNA_INIT_SEQUENCE_STRUCTURE.\n
 * Returns 0 on success or an error code from the utilised functions.
 *
 * @param[in] seq Sequence.
 * @param[in] vienna Structure string in Vienna notation.
 * @param[in] length Length of the RNA.
 * @param[in] sigma Alphabet.
 * @param[in] this Rna data object.
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
int
rna_init_sequence_structure (const char* seq,
                             const char* vienna,
                             const unsigned long length,
                             Alphabet* sigma,
                             Rna* this,
                             const char* file, const int line)
{
   int error = 0;

   /* allocate and store sequence */
   error = rna_init_sequence (seq, length, sigma, this, file, line);

   /* allocate and store structure */
   if (!error)
   {
      error = rna_init_pairlist_vienna (vienna, length, this, file, line);
   }

   return error;
}

/** @brief Decompose an already stored secondary structue.
 *
 * Initialise the secondary structure component of an Rna data object. Instead
 * of a simple list of pairs, this creates a set of structural components the
 * RNA structure consists of. If compiled with enabled memory checking,
 * @c file and @c line should point to the position where the function was
 * called. Both parameters are automatically set by using the macro
 * @c RNA_SECSTRUCT_INIT.\n
 * Returns 0 on success, @c ERR_RNA_ALLOC if memory allocation failed.
 *
 * @param[in] this Rna object.
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
int
rna_secstruct_init (Rna* this, const char* file, const int line)
{
   int error = 0;

   assert (this);
   assert (this->pairs);
   assert (this->structure == NULL);

   this->structure = secstruct_new (file, line);

   if (this->structure == NULL)
   {
      error = ERR_RNA_ALLOC;
   }
   else
   {
      error = secstruct_find_interactions (this->pairs, this->size,
                                           this->structure);
   }

/*    if (!error) */
/*    { */
/*       mprintf ("Stacked base piars:\n"); */
/*       secstruct_fprintf_stacks (stdout, this->structure); */
/*       mprintf ("\nHairpin loops:\n"); */
/*       secstruct_fprintf_hairpins (stdout, this->structure); */
/*       mprintf ("\nBulge loops:\n"); */
/*       secstruct_fprintf_bulges (stdout, this->structure); */
/*       mprintf ("\nInternal loops:\n"); */
/*       secstruct_fprintf_internals (stdout, this->structure); */
/*       mprintf ("\nExternal loop:\n"); */
/*       secstruct_fprintf_external (stdout, this->structure); */
/*       mprintf ("\nMulti loops:\n"); */
/*       secstruct_fprintf_multiloops (stdout, this->structure); */
/*    } */

   return error;
}

/********************************   Altering   ********************************/

/** @brief Set a certain base at a certain position in an Rna sequence.
 *
 * Set a given base at a certain position in the sequence component of an Rna
 * object. Position counting starts at 0.\n
 *
 * @param[in] base The base type to set.
 * @param[in] pos Position where to place base.
 * @param[in] this Rna data object.
 */
void
rna_set_sequence_base (const char base, const unsigned long pos, Rna* this)
{
   assert (this);
   assert (this->seq);
   assert (this->size > pos);

   this->seq[pos] = base;

}

/** @brief Copy a given sequence into an Rna object.
 *
 * Set the sequence component of an Rna object. Since the input sequence is
 * copied, memory for it has to be allocated before calling. The sequence is
 * stored verbatim, this means that it is not transformed into numbers or vice
 * versa automatically. If @c len exceeds the allocated memory, only a fitting
 * prefix of the sequence is copied. We do not check whether sequence consists
 * only of symbols from a certain RNA alphabet.
 *
 * @param[in] sequence Sequence to be copied.
 * @param[in] len Length of the sequence.
 * @param[in] this Rna object.
 */
void
rna_set_sequence (const char* sequence, const unsigned long len, Rna* this)
{
   assert (this);
   assert (sequence);

   memcpy (this->seq, sequence, (len < this->size ? len : this->size));
}

/** @brief Transform Rna sequence to numbers.
 *
 * Transform the sequence component of an Rna object to numbers.\n
 * Returns 0 on success, @c ERR_RNA_NO_BASE if the sequence contains invalid
 * bases. In case of error, the sequence is unchanged.
 *
 * @param[in] sigma Alphabet to use for transformation.
 * @param[in] this Rna data object.
 */
int
rna_transform_sequence_2_no (const Alphabet* sigma, Rna* this)
{
   unsigned long i;
   char tmp;

   assert (sigma);
   assert (this);
   assert (this->seq);

   for (i = 0; i < this->size; i++)
   {
      tmp = alphabet_base_2_no (this->seq[i], sigma);

      if (tmp == CHAR_UNDEF)
      {
         if (i > 0)
         {
            for (; i > 0; i--)
            {
               this->seq[i - 1] = alphabet_no_2_base (this->seq[i - 1], sigma);
            }
         }
         return ERR_RNA_NO_BASE;
      }

      this->seq[i] = tmp;
   }

   return 0;
}

/** @brief Transform Rna sequence from numbers to bases.
 *
 * Transform the sequence component of an Rna object to bases.\n
 * Returns 0 on success, @c ERR_RNA_NO_BASE if the sequence contains invalid
 * bases. In case of error, the sequence is unchanged.
 *
 * @param[in] sigma Alphabet to use for transformation.
 * @param[in] this Rna data object.
 */
int
rna_transform_sequence_2_bases (const Alphabet* sigma, Rna* this)
{
   unsigned long i;
   char tmp;

   assert (sigma);
   assert (this);
   assert (this->seq);

   for (i = 0; i < this->size; i++)
   {
      tmp = alphabet_no_2_base (this->seq[i], sigma);

      if (tmp == CHAR_UNDEF)
      {
         if (i > 0)
         {
            for (; i > 0; i--)
            {
               this->seq[i - 1] = alphabet_base_2_no (this->seq[i - 1], sigma);
            }
         }
         return ERR_RNA_NO_BASE;
      }

      this->seq[i] = tmp;
   }

   return 0;
}


/*********************************   Access   *********************************/

/** @brief Get the no. of nucleotides of an Rna object.
 *
 * Returns the no. of nucleotides of a given Rna object.
 *
 * @param[in] this Rna data object.
 */
unsigned long
rna_get_size (const Rna* this)
{
   assert (this);

   return this->size;
}

/** @brief Get the list of pairs of an Rna object.
 *
 * Returns the pairing list of an Rna object. If no list was set before,
 * returns @c NULL.
 *
 * @params[in] this Rna data object.
 */
unsigned long*
rna_get_pairlist (const Rna* this)
{
   assert (this);

   return this->pairs;
}

/** @brief Get the sequence of an Rna object.
 *
 * Returns the sequence component of an Rna object as a cstring. If no sequence
 * was set before, returns @c NULL.
 *
 * @params[in] this Rna data object.
 */
char*
rna_get_sequence (const Rna* this)
{
   assert (this);

   return this->seq;
}

/** @brief Get a base from a certain position in an Rna sequence.
 *
 * Get a base from a certain position in the sequence component of an Rna
 * object. Position counting starts at 0.\n
 *
 * @param[in] pos Position where to place base.
 * @param[in] this Rna data object.
 */
char
rna_get_sequence_base (const unsigned long pos, const Rna* this)
{
   assert (this);
   assert (this->seq);
   assert (this->size > pos);

   return this->seq[pos];

}

/** @brief Get the pairing partner of a base.
 *
 * Returns the position of a pairing partner of a position in the RNA sequence.
 * Counting starts at 1, @c NOT_PAIRED is returned for unpaired positions.\n
 *
 * @params[in] pos Position to be evaluated.
 * @params[in] this Rna data object.
 */
unsigned long
rna_base_pairs_with (const unsigned long pos, const Rna* this)
{
   assert (this);
   assert (this->pairs);
   assert (this->size > pos);

   return this->pairs[pos];
}

/** @brief Check all base pairs of a structure to be valid.
 *
 * For a given secondary structure and sequence, test all base pairs to be
 * allowed. This should be used to verify if a RNA secondary structure can be
 * evaluated with a certain scoring function.\n
 * Returns the size of the RNA if all base pairs are allowed or the index of
 * the 5' base of the first base pair not allowed.
 *
 * @param[in] this RNA structure to be checked.
 */
unsigned long
rna_validate_basepairs (bool (*validate_basepair) (const char, const char,
                                                   void*),
                        void* scores, const Rna* this)
{
   unsigned long k;

   assert (validate_basepair);
   assert (this);
   assert (this->seq);
   assert (this->pairs);

   for (k = 0; k < this->size; k++)
   {
      if (this->pairs[k] != NOT_PAIRED)
      {
         if (validate_basepair (this->seq[k], this->seq[this->pairs[k]], scores)
             == false)
         {
            return k;
         }  
      }
   }

   return this->size;
}

unsigned long
rna_validate_basepairs_nn_scores (NN_scores* scores, const Rna* this)
{
   return rna_validate_basepairs (nn_scores_is_allowed_basepair, scores, this);
}

int
rna_secstruct_calculate_DG (const NN_scores* scores, const Rna* this)
{
   assert (this);

   return secstruct_calculate_DG (this->seq, scores, this->structure);
}
