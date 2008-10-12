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
#include "rna.h"


typedef struct {
      ArrayUlong tetra_loop;
} struct_comp;


struct Rna {
      char*          seq;       /* the nucleotide sequence */
      /* char* vienna; */       /* vienna string */
      unsigned long* pairs;     /* base pairs */
      unsigned long  size;      /* size of the RNA (sequence & 2D structure) */
      struct_comp    structure; /* decomposed structure */
};


/**********************   Constructors and destructors   **********************/

/** @brief Init a struct_comp structure.
 *
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
static __inline__ struct_comp
struct_comp_new (void/* const char* file, const int line */)
{
   struct_comp this;

   ARRAY_ULONG_INIT(this.tetra_loop, 0);

   return this;
}

/** @brief Delete the components of a struct_comp structure.
 *
 * @param[in] this Object to be freed.
 */
static __inline__ void
struct_comp_delete (struct_comp* this)
{
   if (this != NULL)
   {
     ARRAY_DELETE(this->tetra_loop);
   }
}

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
   Rna* this = XOBJ_MALLOC(sizeof (Rna), file, line);
   
   if (this != NULL)
   {
      this->size      = 0;
      this->seq       = NULL;
      /* this->vienna    = NULL; */
      this->pairs     = NULL;
      this->structure = struct_comp_new();
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
     struct_comp_delete (&this->structure);
     XFREE(this);
   }
}

/** @brief Allocate memory for the sequence component of an @c Rna object.
 *
 * Allocate the memory for a sequence stored within an @c Rna dat aobject. If
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
   this->pairs = XOBJ_MALLOC(sizeof (this->pairs[0]) * length, file, line);
   memset (this->pairs, INT_MAX, sizeof (this->pairs[0]) * length);
   if (this->pairs == NULL)
   {
      error = ERR_RNA_ALLOC;
   }

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
