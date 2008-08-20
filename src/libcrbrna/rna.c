/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
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
#include <stddef.h>
#include <libcrbbasic/crbbasic.h>
#include "rna.h"


struct Rna {
      char* seq;                /* the nucleotide sequence */
      char* vienna;             /* vienna string */
      unsigned long* pairs;     /* base pairs */
      unsigned long size;       /* size of the RNA (sequence & 2D structure) */
};


/**********************   Constructors and destructors   **********************/

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
      this->size   = 0;
      this->seq    = NULL;
      this->vienna = NULL;
      this->pairs  = NULL;
   }

   return this;
}

/** @brief Delete a Rna.
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
     XFREE(this->vienna);
     XFREE(this->pairs);
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
rna_allocate_sequence (const unsigned long size, Rna* this,
                       const char* file, const int line)
{
   assert (this);
   assert (this->seq == NULL);
   assert ((this->size == 0) || (this->size == size));

   this->size = size;

   this->seq = XOBJ_MALLOC (this->size * sizeof (this->seq[0]), file, line);
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
   assert (this->vienna == NULL);
   assert (this->pairs == NULL);
   assert ((this->size == 0) || (this->size == length));

   /* allocate list */
   /* !!! by using calloc we assume, that NOT_PAIRED equals zero !!! */
   this->pairs = XOBJ_CALLOC (length, sizeof (this->pairs[0]), file, line);
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
               this->pairs[i] = p_stack[open_pair] + 1;
               this->pairs[p_stack[open_pair]] = i + 1;
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


/*********************************   Access   *********************************/

/** @brief Get the pairing partner of a base.
 *
 * Returns the position of a pairing partner of a position in the RNA sequence.
 * Counting starts at 1, @c NOT_PAIRED is returned for unpaired positions.\n
 *
 * @params[in] pos Position to be evaluated.
 * @params[in] this Rna data object.
 */
unsigned long
rna_base_pairs_with (unsigned long pos, Rna* this)
{
   assert (this);
   assert (this->pairs);
   assert (this->size > pos);

   return this->pairs[pos];
}
