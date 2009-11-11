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
 *  @file libcrbrna/secstruct.c
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
 *         - 2009Oct19 bienert: Added new mapping, sequence position ->
 *                              secondary structure element, makes it possible
 *                              to query the structural feature a certain
 *                              position in the sequence belong to.
 *
 */


#include <config.h>
#include <stdlib.h> /* only as long as EXIT_FAILURE is in use */
#include <math.h>
#include <libcrbbasic/crbbasic.h>
#include "rna.h"
#include "secstruct.h"
#include "nn_scores.h"

/* structure to store hairpins */
typedef struct {
      unsigned long i;          /* 5' end of loop */
      unsigned long j;          /* 3' end of loop */
      unsigned long size;       /* loop size */
} HairpinLoop;

ARRAY_CREATE_CLASS(HairpinLoop);

/* structure to store stacks */
typedef struct {
/* TODO: Is it sensible to store both pairs here? E.g. for working on a duplex?
         5' ------------------- 3' in this case, the last stack is not followed
             |||||||||||||||||     by a closing bp of some loop!!!
         3' ------------------- 5'
*/
      unsigned long i, j;
} StackLoop;

ARRAY_CREATE_CLASS(StackLoop);

/* structure to store bulge loops */
typedef struct {
      unsigned long i1, j1, i2, j2;
      unsigned long size;
} BulgeLoop;

ARRAY_CREATE_CLASS(BulgeLoop);

/* structure to store internal loops */
typedef struct {
      unsigned long i1, j1, i2, j2;
      unsigned long size1, size2;
} IntLoop;

ARRAY_CREATE_CLASS(IntLoop);

/* store multiloops */
typedef struct {
      unsigned long unpaired;
      unsigned long nstems, ndangle5, ndangle3;
      unsigned long (*stems)[No_Of_Strands];
      unsigned long (*dangle5)[No_Of_Dangles];
      unsigned long (*dangle3)[No_Of_Dangles];
} MultiLoop;

ARRAY_CREATE_CLASS(MultiLoop);

typedef struct {
   SecStructFtrs type;
   unsigned long idx;
} StructElem;

ARRAY_CREATE_CLASS(StructElem);

/* data structure to store decomposed secondary structures */
struct SecStruct {
   ArrayHairpinLoop hairpin_loop;
   ArrayStackLoop stack;
   ArrayBulgeLoop bulge_loop;
   ArrayIntLoop internal_loop;
   ArrayMultiLoop multi_loop;
   MultiLoop ext_loop;
   /* bienert 19-10-09: sequence -> structure mapping */
   StructElem* seq_struct_map;
   ArrayStructElem multi_map;
   unsigned long size;

   /* operations: delete?
      del: take type & idx, set all seq.idxs to not in secstruct on multiples?
      del: just take idx, individual delete functions
      two del funcs: fst uses 2nd */
};


/**********************   Constructors and destructors   **********************/

/** @brief Create a new secstruct object.
 *
 * The constructor for @c SecStruct objects. If compiled with enabled
 * memory checking, @c file and @c line should point to the position where the
 * function was called. Both parameters are automatically set by using the
 * macro @c SECSTRUCT_NEW.\n
 * Returns @c NULL on error.
 *
 * @param[in] size No. of bases in the structure.
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
SecStruct*
secstruct_new (const char* file, const int line)
{
   SecStruct* this = XOBJ_MALLOC(sizeof (*this), file, line);

   if (this != NULL)
   {
      ARRAY_SET_VOID (this->hairpin_loop);
      ARRAY_SET_VOID (this->stack);
      ARRAY_SET_VOID (this->bulge_loop);
      ARRAY_SET_VOID (this->internal_loop);
      ARRAY_SET_VOID (this->multi_loop);
      ARRAY_SET_VOID (this->multi_map);

      ARRAY_INIT (this->hairpin_loop, 1, HairpinLoop);
      if (ARRAY_IS_NULL (this->hairpin_loop))
      {
         secstruct_delete (this);
         return NULL;
      }

      ARRAY_INIT (this->stack, 1, StackLoop);
      if (ARRAY_IS_NULL (this->stack))
      {
         secstruct_delete (this);
         return NULL;
      }

      ARRAY_INIT (this->bulge_loop, 1, BulgeLoop);
      if (ARRAY_IS_NULL (this->bulge_loop))
      {
         secstruct_delete (this);
         return NULL;
      }

      ARRAY_INIT (this->internal_loop, 1, IntLoop);
      if (ARRAY_IS_NULL (this->internal_loop))
      {
         secstruct_delete (this);
         return NULL;
      }

      ARRAY_INIT (this->multi_loop, 1, MultiLoop);
      if (ARRAY_IS_NULL (this->multi_loop))
      {
         secstruct_delete (this);
         return NULL;
      }

      this->seq_struct_map = NULL;
      this->size = 0;

      ARRAY_INIT (this->multi_map, 1, StructElem);
      if (ARRAY_IS_NULL (this->multi_map))
      {
         secstruct_delete (this);
         return NULL;
      }
   }

   return this;
}

static int
multiloop_allocate (unsigned long size, MultiLoop* ml)
{
   assert (ml);

   ml->nstems = ml->ndangle5 = ml->ndangle3 = 0;

   ml->stems = XMALLOC (size * sizeof (*ml->stems));
   if (ml->stems == NULL)
   {
      ml->dangle5 = NULL;
      ml->dangle3 = NULL;
      return 1;
   }

   ml->dangle5 = XMALLOC (size * sizeof (*ml->dangle5));
   if (ml->dangle5 == NULL)
   {
      ml->dangle3 = NULL;
      return 1;
   }

   ml->dangle3 = XMALLOC (size * sizeof (*ml->dangle3));
   if (ml->dangle3 == NULL)
   {
      return 1;
   }

   return 0;
}

static void
multiloop_delete (MultiLoop* ml)
{
   assert (ml);

   XFREE (ml->stems);
   XFREE (ml->dangle5);
   XFREE (ml->dangle3);
}

/** @brief Delete a secstruct object.
 *
 * The destructor for @c SecStruct objects.
 *
 * @param[in] this Object to be freed.
 */
void
secstruct_delete (SecStruct* this)
{
   unsigned long i;

   if (this != NULL)
   {
      ARRAY_DELETE (this->hairpin_loop);
      ARRAY_DELETE (this->stack);
      ARRAY_DELETE (this->bulge_loop);
      ARRAY_DELETE (this->internal_loop);

      XFREE(this->seq_struct_map);
      this->size = 0;

      ARRAY_DELETE (this->multi_map);

      for (i = 0; i < ARRAY_CURRENT (this->multi_loop); i++)
      {
         multiloop_delete(&(ARRAY_ACCESS(this->multi_loop, i)));
      }
      ARRAY_DELETE (this->multi_loop);

      multiloop_delete (&this->ext_loop);

      XFREE(this);
   }
}


/********************************   Altering   ********************************/

static int
s_seq_struct_map_multi_pair (const unsigned long i, const unsigned long j,
                             SecStructFtrs type,
                             SecStruct* this)
{
   StructElem loop;

   if (this->seq_struct_map[i].type == SCSTRCT_VOID)
   {
      this->seq_struct_map[i].type = type;
      this->seq_struct_map[i].idx = ARRAY_CURRENT(this->multi_loop);
   }
   else if (   this->seq_struct_map[i].type != type 
            || this->seq_struct_map[i].idx != ARRAY_CURRENT(this->multi_loop))
   {
      assert (this->seq_struct_map[i].type != SCSTRCT_MTO);

      loop.type = this->seq_struct_map[i].type;
      loop.idx = this->seq_struct_map[i].idx;
      
      this->seq_struct_map[i].type = SCSTRCT_MTO;
      this->seq_struct_map[i].idx = ARRAY_CURRENT(this->multi_map);
      
      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });

      loop.type = type;
      loop.idx = ARRAY_CURRENT(this->multi_loop);

      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });
   }

   if (this->seq_struct_map[j].type == SCSTRCT_VOID)
   {
      this->seq_struct_map[j].type = type;
      this->seq_struct_map[j].idx = ARRAY_CURRENT(this->multi_loop);
   }
   else if (   this->seq_struct_map[j].type != type 
            || this->seq_struct_map[j].idx != ARRAY_CURRENT(this->multi_loop))
   {
      assert (this->seq_struct_map[j].type != SCSTRCT_MTO);

      loop.type = this->seq_struct_map[j].type;
      loop.idx = this->seq_struct_map[j].idx;
      
      this->seq_struct_map[j].type = SCSTRCT_MTO;
      this->seq_struct_map[j].idx = ARRAY_CURRENT(this->multi_map);
      
      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });

      loop.type = type;
      loop.idx = ARRAY_CURRENT(this->multi_loop);

      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });
   }

   return 0;
}

static int
s_seq_struct_map_multi_dangle (unsigned long b, SecStructFtrs type,
                               SecStruct* this)
{
   StructElem loop;

   if (this->seq_struct_map[b].type == SCSTRCT_VOID)
   {
      this->seq_struct_map[b].type = type;
      this->seq_struct_map[b].idx = ARRAY_CURRENT(this->multi_loop);
   }
   else if ( this->seq_struct_map[b].type != type 
           ||this->seq_struct_map[b].idx != ARRAY_CURRENT(this->multi_loop))
   {
      assert (this->seq_struct_map[b].type != SCSTRCT_MTO);

      loop.type = this->seq_struct_map[b].type;
      loop.idx = this->seq_struct_map[b].idx;
      
      this->seq_struct_map[b].type = SCSTRCT_MTO;
      this->seq_struct_map[b].idx = ARRAY_CURRENT(this->multi_map);
      
      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });
      
      loop.type = type;
      loop.idx = ARRAY_CURRENT(this->multi_loop);
      
      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });
   }
   return 0;
}

static int
multiloop_find (unsigned long i, unsigned long j, const unsigned long* pairs,
                unsigned long size, MultiLoop* ml, SecStructFtrs type,
                SecStruct* strct)
{
   int error;

   assert (pairs);
   assert (ml);

   ml->unpaired = 0;

   while (i <= j)
   {
      if (pairs[i] == NOT_PAIRED)
      {
         ml->unpaired++;
      }
      else if (i < pairs[i])
      {
         /* add stem */
         ml->stems[ml->nstems][P5_Strand] = i;
         ml->stems[ml->nstems][P3_Strand] = pairs[i];
         ml->nstems++;

         strct->seq_struct_map[i].type = type;
         strct->seq_struct_map[i].idx = ARRAY_CURRENT(strct->multi_loop);
         strct->seq_struct_map[pairs[i]].type = type;
         strct->seq_struct_map[pairs[i]].idx = ARRAY_CURRENT(strct->multi_loop);

         error = s_seq_struct_map_multi_pair (i, pairs[i], type, strct);
         if (error)
            return error;

         if (i > 0)
         {
            /* add dangle5 */
            ml->dangle5[ml->ndangle5][P5_Dangle] = i;
            ml->dangle5[ml->ndangle5][P3_Dangle] = pairs[i];
            ml->dangle5[ml->ndangle5][Ne_Dangle] = i - 1;
            ml->ndangle5++;

            error = s_seq_struct_map_multi_dangle (i - 1, type, strct);
            if (error)
               return error;
         }
        
         if (pairs[i] < size - 1)
         {
            /* add dangle3 */
            ml->dangle3[ml->ndangle3][P5_Dangle] = i;
            ml->dangle3[ml->ndangle3][P3_Dangle] = pairs[i];
            ml->dangle3[ml->ndangle3][Ne_Dangle] = pairs[i] + 1;
            ml->ndangle3++;

            error = s_seq_struct_map_multi_dangle (pairs[i] + 1, type, strct);
            if (error)
               return error;
         }

         i = pairs[i];
      }
      else
      {
         THROW_ERROR_MSG ("Error not supposed");
         exit (EXIT_FAILURE);
      }
      i++;
   }

   return 0;
}

static int
s_seq_struct_map_multi (MultiLoop* ml, SecStructFtrs type, SecStruct* this)
{
   unsigned long lp, k, l;

   /* mark unpaired bases */
   lp = ml->stems[ml->nstems - 1][P5_Strand] + 1;

   for (k = 0; k < ml->nstems; k++)
   {
      for (l = lp; l < ml->stems[k][P5_Strand]; l++)
      {
         if (this->seq_struct_map[l].type == SCSTRCT_VOID)
         {
            this->seq_struct_map[l].type = type;
            this->seq_struct_map[l].idx = ARRAY_CURRENT(this->multi_loop);
         }
         assert (this->seq_struct_map[l].type == type
             && this->seq_struct_map[l].idx == ARRAY_CURRENT(this->multi_loop));
      }
      lp = ml->stems[k][P3_Strand] + 1;
   }

   return 0;
}

/* i: seqpos, j: pairing partner */
static int
s_seq_struct_map_hairpin (unsigned long i, unsigned long j, SecStruct* this)
{
   unsigned long k;
   StructElem loop;

   if (this->seq_struct_map[i].type == SCSTRCT_VOID)
   {
      this->seq_struct_map[i].type = SCSTRCT_HAIRPIN;
      this->seq_struct_map[i].idx = ARRAY_CURRENT(this->hairpin_loop);
      this->seq_struct_map[j].type = SCSTRCT_HAIRPIN;
      this->seq_struct_map[j].idx = ARRAY_CURRENT(this->hairpin_loop);
   }
   else
   {
      assert (this->seq_struct_map[i].type != SCSTRCT_MTO);

      loop.type = this->seq_struct_map[i].type;
      loop.idx = this->seq_struct_map[i].idx;

      this->seq_struct_map[i].type = SCSTRCT_MTO;
      this->seq_struct_map[i].idx = ARRAY_CURRENT(this->multi_map);
      this->seq_struct_map[j].type = SCSTRCT_MTO;
      this->seq_struct_map[j].idx = ARRAY_CURRENT(this->multi_map);

      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });

      loop.type = SCSTRCT_HAIRPIN;
      loop.idx = ARRAY_CURRENT(this->hairpin_loop);

      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });
   }

   /* set unpaired bases */
   for (k = i + 1; k < j; k++)
   {
      assert (this->seq_struct_map[k].type == SCSTRCT_VOID);

      this->seq_struct_map[k].type = SCSTRCT_HAIRPIN;
      this->seq_struct_map[k].idx = ARRAY_CURRENT(this->hairpin_loop);
   }

   return 0;
}

static int
s_seq_struct_map_stack (unsigned long i, unsigned long j, SecStruct* this)
{
   StructElem stack;

   if (this->seq_struct_map[i].type == SCSTRCT_VOID)
   {
      this->seq_struct_map[i].type = SCSTRCT_STACK;
      this->seq_struct_map[i].idx = ARRAY_CURRENT(this->stack);
      this->seq_struct_map[j].type = SCSTRCT_STACK;
      this->seq_struct_map[j].idx = ARRAY_CURRENT(this->stack);
   }
   else
   {
      assert (this->seq_struct_map[i].type != SCSTRCT_MTO);

      stack.type = this->seq_struct_map[i].type;
      stack.idx = this->seq_struct_map[i].idx;

      this->seq_struct_map[i].type = SCSTRCT_MTO;
      this->seq_struct_map[i].idx = ARRAY_CURRENT(this->multi_map);
      this->seq_struct_map[j].type = SCSTRCT_MTO;
      this->seq_struct_map[j].idx = ARRAY_CURRENT(this->multi_map);

      ARRAY_PUSH(this->multi_map, stack, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });

      stack.type = SCSTRCT_STACK;
      stack.idx = ARRAY_CURRENT(this->stack);

      ARRAY_PUSH(this->multi_map, stack, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });
   }

   return 0;
}

static int
s_seq_struct_map_bulge (unsigned long pi1, unsigned long pj1,
                        unsigned long pi2, unsigned long pj2,
                        unsigned long size1, unsigned long size2,
                        SecStruct* this)
{
   unsigned long j;
   StructElem loop;

   /* fst base pair */
   if (this->seq_struct_map[pi1].type == SCSTRCT_VOID)
   {
      this->seq_struct_map[pi1].type = SCSTRCT_BULGE;
      this->seq_struct_map[pi1].idx =
         ARRAY_CURRENT(this->bulge_loop);
      this->seq_struct_map[pj1].type = SCSTRCT_BULGE;
      this->seq_struct_map[pj1].idx =
         ARRAY_CURRENT(this->bulge_loop);
   }
   else
   {
      assert (this->seq_struct_map[pi1].type != SCSTRCT_MTO);

      loop.type = this->seq_struct_map[pi1].type;
      loop.idx = this->seq_struct_map[pi1].idx;

      this->seq_struct_map[pi1].type = SCSTRCT_MTO;
      this->seq_struct_map[pi1].idx = ARRAY_CURRENT(this->multi_map);
      this->seq_struct_map[pj1].type = SCSTRCT_MTO;
      this->seq_struct_map[pj1].idx = ARRAY_CURRENT(this->multi_map);

      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });

      loop.type = SCSTRCT_BULGE;
      loop.idx = ARRAY_CURRENT(this->bulge_loop);

      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });
   }

   /* snd base pair */
   if (this->seq_struct_map[pi2].type == SCSTRCT_VOID)
   {
      this->seq_struct_map[pi2].type = SCSTRCT_BULGE;
      this->seq_struct_map[pi2].idx = ARRAY_CURRENT(this->bulge_loop);
      this->seq_struct_map[pj2].type = SCSTRCT_BULGE;
      this->seq_struct_map[pj2].idx = ARRAY_CURRENT(this->bulge_loop);
   }
   else
   {
      assert (this->seq_struct_map[pi2].type != SCSTRCT_MTO);

      loop.type = this->seq_struct_map[pi2].type;
      loop.idx = this->seq_struct_map[pi2].idx;

      this->seq_struct_map[pi2].type = SCSTRCT_MTO;
      this->seq_struct_map[pi2].idx = ARRAY_CURRENT(this->multi_map);
      this->seq_struct_map[pj2].type = SCSTRCT_MTO;
      this->seq_struct_map[pj2].idx = ARRAY_CURRENT(this->multi_map);

      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });

      loop.type = SCSTRCT_BULGE;
      loop.idx = ARRAY_CURRENT(this->bulge_loop);

      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });
   }

   if (size1 < size2)
   {
      pi1 = pj2;
      pi2 = pj1;
   }
   
   for (j = pi1 + 1; j < pi2; j++)
   {
      assert (this->seq_struct_map[j].type == SCSTRCT_VOID);

      this->seq_struct_map[j].type = SCSTRCT_BULGE;
      this->seq_struct_map[j].idx =
         ARRAY_CURRENT(this->bulge_loop);
   }

   return 0;
}

static int
s_seq_struct_map_internal (unsigned long pi1, unsigned long pj1,
                           unsigned long pi2, unsigned long pj2,
                           SecStruct* this)
{
   unsigned long i;
   StructElem loop;

   /* fst base pair */
   if (this->seq_struct_map[pi1].type == SCSTRCT_VOID)
   {
      this->seq_struct_map[pi1].type = SCSTRCT_INTERNAL;
      this->seq_struct_map[pi1].idx = ARRAY_CURRENT(this->internal_loop);
      this->seq_struct_map[pj1].type = SCSTRCT_INTERNAL;
      this->seq_struct_map[pj1].idx = ARRAY_CURRENT(this->internal_loop);
   }
   else
   {
      assert (this->seq_struct_map[pi1].type != SCSTRCT_MTO);

      loop.type = this->seq_struct_map[pi1].type;
      loop.idx = this->seq_struct_map[pi1].idx;

      this->seq_struct_map[pi1].type = SCSTRCT_MTO;
      this->seq_struct_map[pi1].idx = ARRAY_CURRENT(this->multi_map);
      this->seq_struct_map[pj1].type = SCSTRCT_MTO;
      this->seq_struct_map[pj1].idx = ARRAY_CURRENT(this->multi_map);

      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });

      loop.type = SCSTRCT_INTERNAL;
      loop.idx = ARRAY_CURRENT(this->internal_loop);

      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });
   }

   /* snd base pair */
   if (this->seq_struct_map[pi2].type == SCSTRCT_VOID)
   {
      this->seq_struct_map[pi2].type = SCSTRCT_INTERNAL;
      this->seq_struct_map[pi2].idx = ARRAY_CURRENT(this->internal_loop);
      this->seq_struct_map[pj2].type = SCSTRCT_INTERNAL;
      this->seq_struct_map[pj2].idx = ARRAY_CURRENT(this->internal_loop);
   }
   else
   {
      assert (this->seq_struct_map[pi2].type != SCSTRCT_MTO);

      loop.type = this->seq_struct_map[pi2].type;
      loop.idx = this->seq_struct_map[pi2].idx;

      this->seq_struct_map[pi2].type = SCSTRCT_MTO;
      this->seq_struct_map[pi2].idx = ARRAY_CURRENT(this->multi_map);
      this->seq_struct_map[pj2].type = SCSTRCT_MTO;
      this->seq_struct_map[pj2].idx = ARRAY_CURRENT(this->multi_map);

      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });

      loop.type = SCSTRCT_INTERNAL;
      loop.idx = ARRAY_CURRENT(this->internal_loop);

      ARRAY_PUSH(this->multi_map, loop, StructElem,
                 { return ERR_SCSTRCT_ALLOC; });
   }

   for (i = pi1 + 1; i < pi2; i++)
   {
      assert (this->seq_struct_map[i].type == SCSTRCT_VOID);

      this->seq_struct_map[i].type = SCSTRCT_INTERNAL;
      this->seq_struct_map[i].idx = ARRAY_CURRENT(this->internal_loop);
   }
   
   for (i = pj2 + 1; i < pj1; i++)
   {
      assert (this->seq_struct_map[i].type == SCSTRCT_VOID);

      this->seq_struct_map[i].type = SCSTRCT_INTERNAL;
      this->seq_struct_map[i].idx = ARRAY_CURRENT(this->internal_loop);
   }

   return 0;
}

/** @brief Decompose an RNA structure stored.
 *
 * Analyses a RNA structure and stores information about its components in the
 * SecStruct object.\n
 * Returns 0 on success, ... else.
 *
 * @param[in] pairs List of pairing partners, where i pairs pairs[i], an
 *                  unpaired position k is indicated by pairs[k] = NOT_PAIRED
 * @param[in] size Length of the structure
 * @param[in] this SecStruct data object.
 */
int
secstruct_find_interactions (const unsigned long* pairs,
                             const unsigned long size,
                             SecStruct* this)
{
   unsigned long i, j, k, p, q, size1, size2, pi1, pj1, pi2, pj2;
   unsigned long remains = size;
   int error = 0;

   assert (this);
   assert (this->seq_struct_map == NULL);

   /* allocate seq -> structure map */
   this->size = size;
   this->seq_struct_map =
      XMALLOC(sizeof (this->seq_struct_map[0]) * this->size);
   if (this->seq_struct_map == NULL)
   {
      return ERR_SCSTRCT_ALLOC;
   }
   memset (this->seq_struct_map, 0, /* needs SCSTRCT_VOID == 0 */
           sizeof (this->seq_struct_map[0]) * this->size); 

   /* exterior loop */
   error = multiloop_allocate (remains, &this->ext_loop);
   if ((size > 0) && (!error))
   {
      error = multiloop_find (0, (size - 1),
                              pairs,
                              size,
                              &this->ext_loop,
                              SCSTRCT_EXTERNAL, this);

      if (this->ext_loop.nstems > 0)
      {
         for (j = 0; j < this->ext_loop.stems[0][P5_Strand]; j++)
         {
            this->seq_struct_map[j].type = SCSTRCT_EXTERNAL;
            this->seq_struct_map[j].idx = 0;
         }

         pi1 = this->ext_loop.stems[0][P3_Strand] + 1;
         for (k = 1; k < this->ext_loop.nstems; k++)
         {
            for (j = pi1; j < this->ext_loop.stems[k][P5_Strand]; j++)
            {
               if (this->seq_struct_map[j].type == SCSTRCT_VOID)
               {
                  this->seq_struct_map[j].type = SCSTRCT_EXTERNAL;
                  this->seq_struct_map[j].idx = 0;
               }
               assert (this->seq_struct_map[j].type == SCSTRCT_EXTERNAL);
            }
            pi1 = this->ext_loop.stems[k][P3_Strand] + 1;
         }

         for (j = pi1; j < this->size; j++)
         {
            this->seq_struct_map[j].type = SCSTRCT_EXTERNAL;
            this->seq_struct_map[j].idx = 0;
         }
      }
   }

   /* hairpins, interior loops (including stacking basepairs) and multiloops */
   i = 0;
   while ((i < size) && (!error)) 
   {
      /* continue until we have found an opening basepairing
         i.e. we're at i of basepair (i, pairs[i]) */
      if ((pairs[i] == NOT_PAIRED) || (i > pairs[i]))
      {
         i++;
      }
      else
      {
         /* now search inwards from basepair (i, pairs[i] - 1) for the next two
            basepairs (p, pairs[p]) and ( (pairs[q], q) or (q, pairs[q]) )
            TODO: i think all 4 combinations are possible:
            (p,pairs[p])  (pairs[p],p)  (pairs[q],q)  (q,pairs[q]) */
         
         p = i + 1;
         q = pairs[i] - 1; /* holds for "0 pairs 1" -> pairs[0] = 2 */
         while ((pairs[p] == NOT_PAIRED) && (p < pairs[i]))
            p++;
         while ((pairs[q] == NOT_PAIRED) && (q > i))
            q--;

         /* now we can tell what loop type we have */
         if (q < p)
         {
            /* hairpin loop - pairs have run past each other */
            HairpinLoop hl = { .i = i, .j = pairs[i],
                               .size = pairs[i] - i - 1 };

            /* seq -> structure mapping */
            error = s_seq_struct_map_hairpin (i, pairs[i], this);

            /* store element */
            if (!error)
               ARRAY_PUSH(this->hairpin_loop, hl, HairpinLoop,
                          { error = ERR_SCSTRCT_ALLOC; });
         }
         else
         {
            /* TODO: have all cases really been taken care of here?
               i think so, but need to recheck and document here why */
            if (q == pairs[p])
            {
               /* stacking pair, bulge loop or internal loop - pairs are
                  the same */
               size1 = p - i - 1;
               size2 = pairs[i] - q - 1;
               pi1 = i;
               pj1 = pairs[i];
               pi2 = p;
               pj2 = pairs[p];
               
               if ((size1 == 0) && (size2 == 0))
               {
                  /* stacking pair of base pairs */
                  StackLoop st = { .i = pi1, .j = pj1};

                  /* seq -> structure mapping */
                  error = s_seq_struct_map_stack (pi1, pj1, this);

                  /* store element */
                  if (!error)
                     ARRAY_PUSH(this->stack, st, StackLoop,
                                { error = ERR_SCSTRCT_ALLOC; });
               }
               else if ((size1 == 0) || (size2 == 0))
               {
                  /* bulge loop */
                  BulgeLoop bg = { .i1 = pi1, .j1 = pj1, .i2 = pi2, .j2 = pj2,
                                   .size = (size1 > size2 ? size1 : size2) };

                  /* seq -> structure mapping */
                  error = s_seq_struct_map_bulge (pi1, pj1, pi2, pj2,
                                                  size1, size2, this);

                  /* store element */
                  if (!error)
                     ARRAY_PUSH(this->bulge_loop, bg, BulgeLoop,
                                { error = ERR_SCSTRCT_ALLOC; });
               }
               else
               {
                  /* generic internal loop */
                  IntLoop in = { .i1 = pi1, .j1 = pj1, .i2 = pi2, .j2 = pj2,
                                 .size1 = size1, .size2 = size2};

                  /* seq -> structure mapping */
                  error = s_seq_struct_map_internal (pi1, pj1, pi2, pj2, this);

                  /* store element */
                  if (!error)
                     ARRAY_PUSH(this->internal_loop, in, IntLoop,
                                { error = ERR_SCSTRCT_ALLOC; });
               }
            }
            else
            {
               MultiLoop ml;
               /* multiloop - pairs are different */
               error = multiloop_allocate (remains, &ml);
               if (!error)
               {
                  /* multiloop_find (p, q, pairs, size, &ml); */
                  error = multiloop_find (i + 1,
                                          pairs[i] - 1,
                                          pairs,
                                          size,
                                          &ml,
                                          SCSTRCT_MULTI, 
                                          this);

                  /* TODO: perhaps move this into nn_multiloop_xfind
                     assumes there is enough space in arrays */
                  /* add the basepair initiating the multiloop to the stems */
                  ml.stems[ml.nstems][P5_Strand] = i;
                  ml.stems[ml.nstems][P3_Strand] = pairs[i];
                  ml.nstems++;

                  /* add extra dangles for outer basepair of multiloop
                     TODO: this is the way it is done in vienna rna, but is
                     this really right ??? */
                  ml.dangle5[ml.ndangle5][P5_Dangle] = pairs[i];
                  ml.dangle5[ml.ndangle5][P3_Dangle] = i;
                  ml.dangle5[ml.ndangle5][Ne_Dangle] = pairs[i] - 1;
                  ml.ndangle5++;

                  ml.dangle3[ml.ndangle3][P5_Dangle] = pairs[i];
                  ml.dangle3[ml.ndangle3][P3_Dangle] = i;
                  ml.dangle3[ml.ndangle3][Ne_Dangle] = i + 1;
                  ml.ndangle3++;

                  if (!error)
                  {
                     error = s_seq_struct_map_multi_pair (i, pairs[i],
                                                          SCSTRCT_MULTI, this);
                  }
                  if (!error)
                  {
                     error = s_seq_struct_map_multi_dangle (pairs[i] - 1,
                                                            SCSTRCT_MULTI,
                                                            this);
                  }
                  if (!error)
                  {
                     error = s_seq_struct_map_multi_dangle (i + 1,
                                                            SCSTRCT_MULTI, 
                                                            this);
                  }
                  if (!error)
                  {
                     error = s_seq_struct_map_multi (&ml, SCSTRCT_MULTI, this);
                  }
                  if (!error)
                  { 
                     ARRAY_PUSH (this->multi_loop, ml, MultiLoop,
                                 { error = ERR_SCSTRCT_ALLOC; });
                  }
               }
            }
         }
         
         /* move to next paired base */
         i = p;
      }
   }

   /*mprintf ("Stacked base pairs:\n");
   secstruct_fprintf_stacks (stdout, this);
   mprintf ("\nHairpin loops:\n");
   secstruct_fprintf_hairpins (stdout, this);
   mprintf ("\nBulge loops:\n");
   secstruct_fprintf_bulges (stdout, this);
   mprintf ("\nInternal loops:\n");
   secstruct_fprintf_internals (stdout, this);
   mprintf ("\nExternal loop:\n");
   secstruct_fprintf_external (stdout, this);
   mprintf ("\nMulti loops:\n");
   secstruct_fprintf_multiloops (stdout, this);*/

   return error;
}

/** @brief Delete a stacked base pair of a secondary structure
 *
 * Removes a base pair from a structure. @c index points to the stack as stored
 * in @c this, not to a sequence position. It may be obtained using
 * @c secstruct_get_feature_at_pos(). Please note that we do not check for a
 * valid @c index, here.
 *
 * @param[in] index Internal index of stack to be deleted.
 * @param[in] this SecStruct data object.
 */
void
secstruct_delete_stack (const unsigned long index, SecStruct* this)
{
   unsigned long i, j; /* seq.pos. of i and j component */

   assert (this);
   assert (ARRAY_NOT_NULL (this->stack));
   assert (ARRAY_CURRENT (this->stack) > index);
   assert (this->seq_struct_map);

   /* update seq.map for item to be deleted */
   i = ARRAY_ACCESS (this->stack, index).i;
   j = ARRAY_ACCESS (this->stack, index).j;

   if (this->seq_struct_map[i].type == SCSTRCT_STACK)
   {
      this->seq_struct_map[i].type = SCSTRCT_VOID;
   }
   else
   {
      if (   ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).type ==
             SCSTRCT_STACK
          && ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).idx ==
             index)
      {
         this->seq_struct_map[i].type =
            ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx + 1).type;
         this->seq_struct_map[i].idx =
            ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx + 1).idx;
      }
      else
      {
         this->seq_struct_map[i].type =
            ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).type;
         this->seq_struct_map[i].idx =
            ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).idx;
      }
   }

   if (this->seq_struct_map[j].type == SCSTRCT_STACK)
   {
      this->seq_struct_map[j].type = SCSTRCT_VOID;
   }
   else
   {
      if (   ARRAY_ACCESS(this->multi_map, this->seq_struct_map[j].idx).type ==
             SCSTRCT_STACK
          && ARRAY_ACCESS(this->multi_map, this->seq_struct_map[j].idx).idx ==
             index)
      {
         this->seq_struct_map[j].type =
            ARRAY_ACCESS(this->multi_map, this->seq_struct_map[j].idx + 1).type;
         this->seq_struct_map[j].idx =
            ARRAY_ACCESS(this->multi_map, this->seq_struct_map[j].idx + 1).idx;
      }
      else
      {
         this->seq_struct_map[j].type =
            ARRAY_ACCESS(this->multi_map, this->seq_struct_map[j].idx).type;
         this->seq_struct_map[j].idx =
            ARRAY_ACCESS(this->multi_map, this->seq_struct_map[j].idx).idx;
      }
   }

   /* swap entries & shrink array */
   ARRAY_CURRENT_DECREMENT(this->stack);
   ARRAY_SET(this->stack,
             ARRAY_ACCESS(this->stack, ARRAY_CURRENT (this->stack)),
             index);

   /* update seq.map for swaped entry */
   i = ARRAY_ACCESS (this->stack, index).i;
   j = ARRAY_ACCESS (this->stack, index).j;

   /* we have moved the last entry of the stack-array to the pos.2b.deleted.
      Now we have to update the seq.pos map:
      - if the pos. is only bound to a single feature, it must be a stack so we
        just set the new pos. in the array
      - if pos. points to a MTO site, we have to find the stack element and
        update the index, there*/
   if (this->seq_struct_map[i].type == SCSTRCT_STACK)
   {
      this->seq_struct_map[i].idx = index;
   }
   else
   {
      if (   ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).type ==
             SCSTRCT_STACK
          && ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).idx ==
             ARRAY_CURRENT (this->stack))
      {
         ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).idx = index;
      }
      else
      {
         ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx + 1).idx =
            index;
      }
   }

   if (this->seq_struct_map[j].type == SCSTRCT_STACK)
   {
      this->seq_struct_map[j].idx = index;
   }
   else
   {
      if (   ARRAY_ACCESS(this->multi_map, this->seq_struct_map[j].idx).type ==
             SCSTRCT_STACK
          && ARRAY_ACCESS(this->multi_map, this->seq_struct_map[j].idx).idx ==
             ARRAY_CURRENT (this->stack))
      {
         ARRAY_ACCESS(this->multi_map, this->seq_struct_map[j].idx).idx = index;
      }
      else
      {
         ARRAY_ACCESS(this->multi_map, this->seq_struct_map[j].idx + 1).idx =
            index;
      }
   }
}

/** @brief Delete a structural feature of a secondary structure 
 *
 * If the type and the index of a feature in a structure is known, this
 * function may be used for deletion. Index and type can be queried using
 * @c secstruct_get_feature_at_pos(). Please note that this function does not
 * check whether a valid @c index or @c feature is provided.
 * - needs ...
 * - returns ...
 */
void
secstruct_delete_element (const unsigned long index,
                          const SecStructFtrs feature,
                          SecStruct* this)
{
   assert (this);
   assert (feature < N_FEATURES);

   switch (feature)
   {
      case SCSTRCT_VOID:
         /* nothing to be done for non-classified position */
         return;
         break;
         /* SCSTRCT_HAIRPIN */
      case SCSTRCT_STACK:
         secstruct_delete_stack (index, this);
         break;
         /*  SCSTRCT_BULGE  
           SCSTRCT_INTERNAL
           SCSTRCT_MULTI  
           SCSTRCT_EXTERNAL
           SCSTRCT_MTO    
           N_FEATURES*/   
      default:
         THROW_ERROR_MSG ("Error not supposed");
         exit (EXIT_FAILURE);
   }
}

/*********************************   Access   *********************************/

/** @brief Get the type of secondary structure for a certain sequence position
 *
 * Retrieve the sort of structural element, a certain sequence position belongs
 * to. Additionally, an index usable with an appropriate @c get function is
 * provided.\n
 * Possible types are taken from the @c SecStructFtrs @c enum.\n
 * If @c SCSTRCT_MTO is returned, the complementary sequence position belongs
 * to two distinct classes. The retrived index might then be used with
 * @c secstruct_get_feature_multi_1st and
 * @c secstruct_get_feature_multi_2nd to retrieve the true classes and
 * indices.\n
 * No errors are defined for return values. If a position in the sequence is not assigned to a structural feature, @c SCSTRCT_VOID is returned.
 *
 * @param[in]  i Sequence position.
 * @param[out] index Address of the returned element.
 * @param[in]  this Secondary structure.
 */
SecStructFtrs
secstruct_get_feature_at_pos (const unsigned long i,
                                unsigned long* index,
                                const SecStruct* this)
{
   assert (index);
   assert (this);
   assert (this->seq_struct_map);
   assert (ARRAY_NOT_NULL (this->multi_map));
   assert (i < this->size);

   index[0] = this->seq_struct_map[i].idx;

   return this->seq_struct_map[i].type;
}

/** @brief Get the sturcutral feature for a certain position classified as MTO.
 *
 * For the use after @c secstruct_get_feature_at_pos (), if @c SCSTRCT_MTO
 * was returned. Delivers the first of possibly two sturctural features.
 * Additionally, an index usable with an appropriate @c get function is
 * provided.\n
 * Possible types are taken from the @c SecStructFtrs @c enum.
 * @c SCSTRCT_MTO is not a valid return value for this function.\n
 * Please note: @c i is again the sequence position. This is for your
 * convenience since it enables safe use of the same variable for @c index as
 * for the call of @c secstruct_get_feature_at_pos () and you might just
 * recycle @c i as well from a previous call.\n
 * No errors are defined for return values.
 *
 * @param[in]  i Sequence position, already verified to be multi-classified.
 * @param[out] index Address of the returned element.
 * @param[in]  this Secondary structure.
 */
SecStructFtrs
secstruct_get_feature_multi_1st (const unsigned long i,
                                   unsigned long* index,
                                   const SecStruct* this)
{
   assert (index);
   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_map));
   assert (i < this->size);
   assert (this->seq_struct_map[i].type == SCSTRCT_MTO);

   index[0] = ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).idx;

   return ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).type;
}

/** @brief Get the sturcutral feature for a certain position classified as MTO.
 *
 * For the use after @c secstruct_get_feature_at_pos (), if @c SCSTRCT_MTO
 * was returned. Delivers the second of possibly two sturctural features.
 * Additionally, an index usable with an appropriate @c get function is
 * provided.\n
 * Possible types are taken from the @c SecStructFtrs @c enum.
 * @c SCSTRCT_MTO is is not a valid return value for this function.\n
 * Please note: @c i is again the sequence position. This is for your
 * convenience since it enables safe use of the same variable for @c index as
 * for the call of @c secstruct_get_feature_at_pos () and you might just
 * recycle @c i as well from a previous call.\n
 * No errors are defined for return values.
 *
 * @param[in]  i Sequence position, already verified to be multi-classified.
 * @param[out] index Address of the returned element.
 * @param[in]  this Secondary structure.
 */
SecStructFtrs
secstruct_get_feature_multi_2nd (const unsigned long i,
                                   unsigned long* index,
                                   const SecStruct* this)
{
   assert (index);
   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_map));
   assert (i < this->size);
   assert (this->seq_struct_map[i].type == SCSTRCT_MTO);

   index[0] = ARRAY_ACCESS(this->multi_map,this->seq_struct_map[i].idx + 1).idx;

   return ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx + 1).type;
}

/** @brief Get the no. of hairpin loops of a 2D structure
 *
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_noof_hairpins (const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->hairpin_loop));

   return ARRAY_CURRENT (this->hairpin_loop);
}

/** @brief Retrieve hairpin geometry in one go.
 *
 * @param[in/out] start Start position in sequence.
 * @param[in/out] end   End position in sequence.
 * @param[in/out] size  Size of the loop.
 * @param[in]     i    No. of hairpin.
 * @param[in]     this Secondary structure.
 */
void
secstruct_get_geometry_hairpin (unsigned long* start,
                                unsigned long* end,
                                unsigned long* size,
                                const unsigned long i,
                                const SecStruct* this)
{
   assert(this);
   assert (ARRAY_NOT_NULL (this->hairpin_loop));
   assert (ARRAY_CURRENT (this->hairpin_loop) > i);
   assert(start);
   assert(end);
   assert(size);

   *start = ARRAY_ACCESS(this->hairpin_loop, i).i;
   *end = ARRAY_ACCESS(this->hairpin_loop, i).j;
   *size = ARRAY_ACCESS(this->hairpin_loop, i).size;
}

/** @brief get the start base of the ith hairpin loop of a 2D structure
 *
 * @param[in] i Index of loop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_start_hairpin (const unsigned long i, const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->hairpin_loop));
   assert (ARRAY_CURRENT (this->hairpin_loop) > i);

   return ARRAY_ACCESS(this->hairpin_loop, i).i;
}

/** @brief get the last base of the ith hairpin loop of a 2D structure
 *
 * @param[in] i Index of loop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_end_hairpin (const unsigned long i, const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->hairpin_loop));
   assert (ARRAY_CURRENT (this->hairpin_loop) > i);

   return ARRAY_ACCESS(this->hairpin_loop, i).j;
}

/** @brief get the no. of unpaired bases of the ith hairpin loop
 *
 * @param[in] i Index of loop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_size_hairpin (const unsigned long i, const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->hairpin_loop));
   assert (ARRAY_CURRENT (this->hairpin_loop) > i);

   return ARRAY_ACCESS(this->hairpin_loop, i).size;
}

/** @brief get the no. of stacks of a 2D structure
 *
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_noof_stacks (const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->stack));

   return ARRAY_CURRENT (this->stack);
}

/** @brief get the 5' base of the ith stack of a 2D structure
 *
 * @param[in] i Index of stack.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_5p_stack (const unsigned long i, const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->stack));
   assert (ARRAY_CURRENT (this->stack) > i);

   return ARRAY_ACCESS(this->stack, i).i;
}

/** @brief get the 3' base of the ith stack of a 2D structure
 *
 * @param[in] i Index of stack.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_3p_stack (const unsigned long i, const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->stack));
   assert (ARRAY_CURRENT (this->stack) > i);

   return ARRAY_ACCESS(this->stack, i).j;
}

/** @brief Get the 5' and 3' base of a stem base pair 
 *
 * @param[in] i     5' base container.
 * @param[in] j     3' base container.
 * @param[in] stack no. of stack.
 * @param[in] this  Secondary structure.
 */
void
secstruct_get_i_geometry_stack (unsigned long* i, unsigned long* j,
                                const unsigned long stack,
                                const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->stack));
   assert (ARRAY_CURRENT (this->stack) > stack);
   assert (i);
   assert (j);

   *i = ARRAY_ACCESS(this->stack, stack).i;
   *j = ARRAY_ACCESS(this->stack, stack).j;
}

/** @brief get the no. of bulge loops of a 2D structure
 *
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_noof_bulges (const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->bulge_loop));

   return ARRAY_CURRENT (this->bulge_loop);
}

/** @brief get the start base of the ith bulge loop of a 2D structure
 *
 * @param[in] i Index of loop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_start_bulge (const unsigned long i, const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->bulge_loop));
   assert (ARRAY_CURRENT (this->bulge_loop) > i);

   return ARRAY_ACCESS(this->bulge_loop, i).i1;
}

/** @brief get the last base of the ith bulge loop of a 2D structure
 *
 * @param[in] i Index of loop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_end_bulge (const unsigned long i, const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->bulge_loop));
   assert (ARRAY_CURRENT (this->bulge_loop) > i);

   return ARRAY_ACCESS(this->bulge_loop, i).j1;
}

/** @brief get the no. of unpaired bases of the ith bulge loop
 *
 * @param[in] i Index of loop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_size_bulge (const unsigned long i, const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->bulge_loop));
   assert (ARRAY_CURRENT (this->bulge_loop) > i);

   return ARRAY_ACCESS(this->bulge_loop, i).size;
}

/** @brief Retrieve bulge loop geometry in one go.
 *
 * @param[in/out] i1   Base i of opening pair.
 * @param[in/out] j1   Base j of opening pair.
 * @param[in/out] i2   Base i of closing pair.
 * @param[in/out] j2   Base j of closing pair
 * @param[in/out] size Size of the loop.
 * @param[in]     i    No. of bulge.
 * @param[in]     this Secondary structure.
 */
void
secstruct_get_geometry_bulge (unsigned long* pi1,
                              unsigned long* pj1,
                              unsigned long* pi2,
                              unsigned long* pj2,
                              unsigned long* size,
                              const unsigned long i,
                              const SecStruct* this)
{
   assert(this);
   assert (ARRAY_NOT_NULL (this->bulge_loop));
   assert (ARRAY_CURRENT (this->bulge_loop) > i);
   assert(pi1);
   assert(pj1);
   assert(pi2);
   assert(pj2);
   assert(size);

   *pi1 = ARRAY_ACCESS(this->bulge_loop, i).i1;
   *pj1 = ARRAY_ACCESS(this->bulge_loop, i).j1;
   *pi2 = ARRAY_ACCESS(this->bulge_loop, i).i2;
   *pj2 = ARRAY_ACCESS(this->bulge_loop, i).j2;
   *size = ARRAY_ACCESS(this->bulge_loop, i).size;
}

/** @brief get the no. of internal loops of a 2D structure
 *
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_noof_internals (const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->internal_loop));

   return ARRAY_CURRENT (this->internal_loop);
}

/** @brief get the geometry of a certain internal loop.
 *
 * @param[out] i1    Base i of opening pair.
 * @param[out] j1    Base j of opening pair.
 * @param[out] i2    Base i of closing pair.
 * @param[out] j2    Base j of closing pair.
 * @param[out] size1 Size of first loop.
 * @param[out] size2 Size of 2nd loop.
 * @param[in] i      Index of loop.
 * @param[in] this   Secondary structure.
 */
void
secstruct_get_geometry_internal (unsigned long* pi1,
                                 unsigned long* pj1,
                                 unsigned long* pi2,
                                 unsigned long* pj2,
                                 unsigned long* size1,
                                 unsigned long* size2,
                                 const unsigned long i,
                                 const SecStruct* this)
{
   assert(this);
   assert (ARRAY_NOT_NULL (this->internal_loop));
   assert (ARRAY_CURRENT (this->internal_loop) > i);
   assert(pi1);
   assert(pj1);
   assert(size1);
   assert(pi2);
   assert(pj2);
   assert(size2);

   *pi1    = ARRAY_ACCESS(this->internal_loop, i).i1;
   *pj1    = ARRAY_ACCESS(this->internal_loop, i).j1;
   *size1 = ARRAY_ACCESS(this->internal_loop, i).size1;
   *pi2    = ARRAY_ACCESS(this->internal_loop, i).i2;
   *pj2    = ARRAY_ACCESS(this->internal_loop, i).j2;
   *size2 = ARRAY_ACCESS(this->internal_loop, i).size2;
}

/** @brief get the no. of multiloops of a 2D structure
 *
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_noof_multiloops (const SecStruct* this)
{
   /*unsigned long i;*/

   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_loop));

   return ARRAY_CURRENT (this->multi_loop);
}

/** @brief get the no. of unpaired bases of the ith multiloop  of a 2D structure
 *
 * @param[in] i Index of multiloop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_noof_unpaired_multiloop (unsigned long i, const SecStruct* this)
{
   /*unsigned long i;*/

   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > i);

   return (ARRAY_ACCESS (this->multi_loop, i).unpaired);
}

/** @brief Get the no. of stems involved in a certain multiloop. 
 *
 * Get the no. of base pairs attached to the ith multiloop substructure of
 * an RNA secondary structure.
 *
 * @param[in] i Index of multiloop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_noof_stems_multiloop (unsigned long i, const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > i);

   return (ARRAY_ACCESS (this->multi_loop, i).nstems);
}

/** @brief get the 5' base of the ith stem of the jth multiloop of a structure.
 *
 * @param[in] i Index of stem.
 * @param[in] j Index of multiloop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_5p_stem_multiloop (const unsigned long i,
                                   const unsigned long j,
                                   const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > j);
   assert (ARRAY_ACCESS (this->multi_loop, j).nstems > i);

   return (ARRAY_ACCESS (this->multi_loop, j).stems[i][P5_Strand]);
}

/** @brief get the 3' base of the ith stem of the jth multiloop of a structure.
 *
 * @param[in] i Index of stem.
 * @param[in] j Index of multiloop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_3p_stem_multiloop (const unsigned long i,
                                   const unsigned long j,
                                   const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > j);
   assert (ARRAY_ACCESS (this->multi_loop, j).nstems > i);

   return (ARRAY_ACCESS (this->multi_loop, j).stems[i][P3_Strand]);
}

/** @brief Get the 5' and 3' base position of a stem in an multiloop.
 *
 * @param[out] p5 5' psoition.
 * @param[out] p3 3' position.
 * @param[in]  i    Stem.
 * @param[in]  j    Multiloop.
 * @param[in]  this Secondary structure.
 */
void
secstruct_get_i_stem_multiloop (unsigned long* p5, unsigned long* p3,
                                const unsigned long i,
                                const unsigned long j,
                                const SecStruct* this)
{
   assert(this);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > j);
   assert (ARRAY_ACCESS (this->multi_loop, j).nstems > i);
   assert(p5);
   assert(p3);

   *p5 = ARRAY_ACCESS (this->multi_loop, j).stems[i][P5_Strand];
   *p3 = ARRAY_ACCESS (this->multi_loop, j).stems[i][P3_Strand];
}

/** @brief Get the no. of 5' dangles (???dangling ends???) involved in a certain multiloop. 
 *
 * @param[in] i Index of multiloop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_noof_5pdangles_multiloop (unsigned long i,
                                          const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > i);

   return (ARRAY_ACCESS (this->multi_loop, i).ndangle5);
}

/** @brief get the 5' base of base pair adjacent to the ith 5' dangle of the
    jth multiloop of a structure.
 *
 * @param[in] i Index of 5' dangle.
 * @param[in] j Index of multiloop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_5p_5pdangle_multiloop (const unsigned long i,
                                       const unsigned long j,
                                       const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > j);
   assert (ARRAY_ACCESS (this->multi_loop, j).ndangle5 > i);

   return (ARRAY_ACCESS (this->multi_loop, j).dangle5[i][P5_Dangle]);
}

/** @brief get the 3' base of base pair adjacent to the ith 5' dangle of the
    jth multiloop of a structure.
 *
 * @param[in] i Index of 5' dangle.
 * @param[in] j Index of multiloop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_3p_5pdangle_multiloop (const unsigned long i,
                                       const unsigned long j,
                                       const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > j);
   assert (ARRAY_ACCESS (this->multi_loop, j).ndangle5 > i);

   return (ARRAY_ACCESS (this->multi_loop, j).dangle5[i][P3_Dangle]);
}

/** @brief get the 5' dangling base of the ith 5' dangle of the jth multiloop
    of a structure.
 *
 * @param[in] i Index of 5' dangle.
 * @param[in] j Index of multiloop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_dangle_5pdangle_multiloop (const unsigned long i,
                                           const unsigned long j,
                                           const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > j);
   assert (ARRAY_ACCESS (this->multi_loop, j).ndangle5 > i);

   return (ARRAY_ACCESS (this->multi_loop, j).dangle5[i][Ne_Dangle]);
}

/** @brief Get the ith 5', 3' and free base of a 5' dangling end.
 *
 * @param[out] p5 Container for the 5' base position.
 * @param[out] p3 Container for the 3' base position.
 * @param[out] fb Container for the free base.
 * @param[in]  i  ith 5' dangling end.
 * @param[in]  j  jth multiloop.
 * @param[in]  this Secondary structure.
 */
void
secstruct_get_i_5pdangle_multiloop (unsigned long* p5,
                                    unsigned long* p3,
                                    unsigned long* fb,
                                    const unsigned long i,
                                    const unsigned long j,
                                    const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > j);
   assert (ARRAY_ACCESS (this->multi_loop, j).dangle5 != NULL);
   assert (ARRAY_ACCESS (this->multi_loop, j).ndangle5 > i);
   assert (p5);
   assert (p3);
   assert (fb);

   *p5 = ARRAY_ACCESS (this->multi_loop, j).dangle5[i][P5_Dangle];
   *p3 = ARRAY_ACCESS (this->multi_loop, j).dangle5[i][P3_Dangle];
   *fb = ARRAY_ACCESS (this->multi_loop, j).dangle5[i][Ne_Dangle];
}

/** @brief Get the no. of 3' dangles (???dangling ends???) involved in a certain multiloop. 
 *
 * @param[in] i Index of multiloop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_noof_3pdangles_multiloop (unsigned long i,
                                          const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > i);

   return (ARRAY_ACCESS (this->multi_loop, i).ndangle3);
}

/** @brief get the 3' base of base pair adjacent to the ith 3' dangle of the
    jth multiloop of a structure.
 *
 * @param[in] i Index of 3' dangle.
 * @param[in] j Index of multiloop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_5p_3pdangle_multiloop (const unsigned long i,
                                       const unsigned long j,
                                       const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > j);
   assert (ARRAY_ACCESS (this->multi_loop, j).ndangle3 > i);

   return (ARRAY_ACCESS (this->multi_loop, j).dangle3[i][P5_Dangle]);
}

/** @brief get the 3' base of base pair adjacent to the ith 3' dangle of the
    jth multiloop of a structure.
 *
 * @param[in] i Index of 3' dangle.
 * @param[in] j Index of multiloop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_3p_3pdangle_multiloop (const unsigned long i,
                                       const unsigned long j,
                                       const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > j);
   assert (ARRAY_ACCESS (this->multi_loop, j).ndangle3 > i);

   return (ARRAY_ACCESS (this->multi_loop, j).dangle3[i][P3_Dangle]);
}

/** @brief get the 3' dangling base of the ith 3' dangle of the jth multiloop
    of a structure.
 *
 * @param[in] i Index of 3' dangle.
 * @param[in] j Index of multiloop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_dangle_3pdangle_multiloop (const unsigned long i,
                                           const unsigned long j,
                                           const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > j);
   assert (ARRAY_ACCESS (this->multi_loop, j).ndangle3 > i);

   return (ARRAY_ACCESS (this->multi_loop, j).dangle3[i][Ne_Dangle]);
}

/** @brief Get the ith 5', 3' and free base of a 3' dangling end.
 *
 * @param[out] p5 Container for the 5' base position.
 * @param[out] p3 Container for the 3' base position.
 * @param[out] fb Container for the free base.
 * @param[in]  i  ith 3' dangling end.
 * @param[in]  j  jth multiloop.
 * @param[in]  this Secondary structure.
 */
void
secstruct_get_i_3pdangle_multiloop (unsigned long* p5,
                                    unsigned long* p3,
                                    unsigned long* fb,
                                    const unsigned long i,
                                    const unsigned long j,
                                    const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > j);
   assert (ARRAY_ACCESS (this->multi_loop, j).dangle3 != NULL);
   assert (ARRAY_ACCESS (this->multi_loop, j).ndangle3 > i);
   assert (p5);
   assert (p3);
   assert (fb);

   *p5 = ARRAY_ACCESS (this->multi_loop, j).dangle3[i][P5_Dangle];
   *p3 = ARRAY_ACCESS (this->multi_loop, j).dangle3[i][P3_Dangle];
   *fb = ARRAY_ACCESS (this->multi_loop, j).dangle3[i][Ne_Dangle];
}

/** @brief get the no. of unpaired bases of an external loop of a 2D structure
 *
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_noof_unpaired_extloop (const SecStruct* this)
{
   assert (this);

   return (this->ext_loop.unpaired);
}

/** @brief Get the no. of stems involved in an external loop. 
 *
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_noof_stems_extloop (const SecStruct* this)
{
   assert (this);

   return (this->ext_loop.nstems);
}

/** @brief get the 5' base of the ith stem of the external loop of a structure.
 *
 * @param[in] i Index of stem.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_5p_stem_extloop (const unsigned long i,
                                 const SecStruct* this)
{
   assert (this);
   assert (this->ext_loop.stems);
   assert (this->ext_loop.nstems > i);

   return (this->ext_loop.stems[i][P5_Strand]);
}

/** @brief get the 3' base of the ith stem of the external loop of a structure.
 *
 * @param[in] i Index of stem.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_3p_stem_extloop (const unsigned long i,
                                 const SecStruct* this)
{
   assert (this);
   assert (this->ext_loop.stems);
   assert (this->ext_loop.nstems > i);

   return (this->ext_loop.stems[i][P3_Strand]);
}

/** @brief Get the 5' and 3' base position of a stem in an external loop.
 *
 * @param[out] p5 5' psoition.
 * @param[out] p3 3' position.
 * @param[in]  i    Stem.
 * @param[in]  this Secondary structure.
 */
void
secstruct_get_i_stem_extloop (unsigned long* p5,
                              unsigned long* p3,
                              const unsigned long i,
                              const SecStruct* this)
{
   assert(this);
   assert (this->ext_loop.stems);
   assert (this->ext_loop.nstems > i);
   assert(p5);
   assert(p3);

   *p5 = this->ext_loop.stems[i][P5_Strand];
   *p3 = this->ext_loop.stems[i][P3_Strand];
}


/** @brief Get the no. of 5' dangles involved in an external loop. 
 *
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_noof_5pdangles_extloop (const SecStruct* this)
{
   assert (this);

   return (this->ext_loop.ndangle5);
}

/** @brief Get the no. of 3' dangles involved in an external loop. 
 *
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_noof_3pdangles_extloop (const SecStruct* this)
{
   assert (this);

   return (this->ext_loop.ndangle3);
}

/** @brief get the 5' base of base pair adjacent to the ith 3' dangle of the
    external loop of a structure.
 *
 * @param[in] i Index of 3' dangle.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_5p_3pdangle_extloop (const unsigned long i,
                                     const SecStruct* this)
{
   assert (this);
   assert (this->ext_loop.dangle3);
   assert (this->ext_loop.ndangle3 > i);

   return (this->ext_loop.dangle3[i][P5_Dangle]);
}

/** @brief get the 3' base of base pair adjacent to the ith 3' dangle of the
    external loop of a structure.
 *
 * @param[in] i Index of 3' dangle.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_3p_3pdangle_extloop (const unsigned long i,
                                     const SecStruct* this)
{
   assert (this);
   assert (this->ext_loop.dangle3);
   assert (this->ext_loop.ndangle3 > i);

   return (this->ext_loop.dangle3[i][P3_Dangle]);
}

/** @brief get the 5' base of base pair adjacent to the ith 5' dangle of the
    external loop of a structure.
 *
 * @param[in] i Index of 3' dangle.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_5p_5pdangle_extloop (const unsigned long i,
                                     const SecStruct* this)
{
   assert (this);
   assert (this->ext_loop.dangle5);
   assert (this->ext_loop.ndangle5 > i);

   return (this->ext_loop.dangle5[i][P5_Dangle]);
}

/** @brief get the 3' base of base pair adjacent to the ith 5' dangle of the
    external loop of a structure.
 *
 * @param[in] i Index of 5' dangle.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_3p_5pdangle_extloop (const unsigned long i,
                                     const SecStruct* this)
{
   assert (this);
   assert (this->ext_loop.dangle5);
   assert (this->ext_loop.ndangle5 > i);

   return (this->ext_loop.dangle5[i][P3_Dangle]);
}

/** @brief Get the ith 5', 3' and free base of a 5' dangling end.
 *
 * @param[out] p5 Container for the 5' base position.
 * @param[out] p3 Container for the 3' base position.
 * @param[out] fb Container for the free base.
 * @param[in]  i  ith 5' dangling end.
 * @param[in]  this Secondary structure.
 */
void
secstruct_get_i_5pdangle_extloop (unsigned long* p5,
                                   unsigned long* p3,
                                   unsigned long* fb,
                                   const unsigned long i,
                                   const SecStruct* this)
{
   assert (this);
   assert (this->ext_loop.dangle5);
   assert (this->ext_loop.ndangle5 > i);
   assert (p5);
   assert (p3);
   assert (fb);

   *p5 = this->ext_loop.dangle5[i][P5_Dangle];
   *p3 = this->ext_loop.dangle5[i][P3_Dangle];
   *fb = this->ext_loop.dangle5[i][Ne_Dangle];
}

/** @brief Get the ith 5', 3' and free base of a 3' dangling end.
 *
 * @param[out] p5 Container for the 5' base position.
 * @param[out] p3 Container for the 3' base position.
 * @param[out] fb Container for the free base.
 * @param[in]  i  ith 3' dangling end.
 * @param[in]  this Secondary structure.
 */
void
secstruct_get_i_3pdangle_extloop (unsigned long* p5,
                                  unsigned long* p3,
                                  unsigned long* fb,
                                  const unsigned long i,
                                  const SecStruct* this)
{
   assert (this);
   assert (this->ext_loop.dangle3);
   assert (this->ext_loop.ndangle3 > i);
   assert (p5);
   assert (p3);
   assert (fb);

   *p5 = this->ext_loop.dangle3[i][P5_Dangle];
   *p3 = this->ext_loop.dangle3[i][P3_Dangle];
   *fb = this->ext_loop.dangle3[i][Ne_Dangle];
}

/*********************************    Misc    *********************************/

/** @brief Calculate the free energy of a structure.
 *
 * Evaluate the free energy of a secondary structure using the Nearest
 * Neighbour model. This function does not check for base pairs allowed by the
 * energy model. So feeding a structure with non canonical base pairs to it
 * will result in undefined behaviour.\n
 * Returns the free energy.
 * 
 * 
 */
int
secstruct_calculate_DG (const char* seq, const NN_scores* scores,
                        const SecStruct* this)
{
   unsigned long k;
   int G, Ge = 0, Gs = 0, Gb = 0, Gh = 0, Gi = 0, Gm = 0;
   char i, j, jm1, ip1;

   assert (seq);
   assert (scores);
   assert (this);

   /* external loop */
   Ge += nn_scores_get_G_extloop_multiloop (seq,
                                            this->ext_loop.unpaired,
                                            this->ext_loop.nstems,
                                            this->ext_loop.stems,
                                            this->ext_loop.ndangle5,
                                            this->ext_loop.dangle5,
                                            this->ext_loop.ndangle3,
                                            this->ext_loop.dangle3,
                                            false,
                                            scores);

   /* stacking pairs */
      /* 5' - ii+1
              jj-1 - 3' */
   for (k = 0; k < ARRAY_CURRENT (this->stack); k++)
   {
      /* Why can we access j-1 and i+1 without any checking here?
         If i and j pair, i is at least 0 and j 1. Hence sequence length is at
         least 2. */
      /* Why do we know that j-1 and i+1 form a pair? If they were unpaired,
         i,j would be the starting base pair of a loop and therefore not in the
         stack array. */
      i   = seq[ARRAY_ACCESS (this->stack, k).i];
      j   = seq[ARRAY_ACCESS (this->stack, k).j];
      jm1 = seq[ARRAY_ACCESS (this->stack, k).j - 1];
      ip1 = seq[ARRAY_ACCESS (this->stack, k).i + 1];

      Gs += nn_scores_get_G_stack (i, j, jm1, ip1, scores);
   }

   /* bulge loops */
   for (k = 0; k < ARRAY_CURRENT (this->bulge_loop); k++)
   {
      Gb += nn_scores_get_G_bulge_loop (
                                     seq[ARRAY_ACCESS (this->bulge_loop, k).i1],
                                     seq[ARRAY_ACCESS (this->bulge_loop, k).j1],
                                     seq[ARRAY_ACCESS (this->bulge_loop, k).i2],
                                     seq[ARRAY_ACCESS (this->bulge_loop, k).j2],
                                     ARRAY_ACCESS (this->bulge_loop, k).size,
                                        scores);
   }

   /* internal loops */
   for (k = 0; k < ARRAY_CURRENT (this->internal_loop); k++)
   {
      Gi += nn_scores_get_G_internal_loop (seq,
                                    ARRAY_ACCESS (this->internal_loop, k).size1,
                                    ARRAY_ACCESS (this->internal_loop, k).size2,
                                    ARRAY_ACCESS (this->internal_loop, k).i1,
                                    ARRAY_ACCESS (this->internal_loop, k).j1,
                                    ARRAY_ACCESS (this->internal_loop, k).i2,
                                    ARRAY_ACCESS (this->internal_loop, k).j2,
                                           scores);
   }

   /* hairpins */
   for (k = 0; k < ARRAY_CURRENT (this->hairpin_loop); k++)
   {
      Gh += nn_scores_get_G_hairpin_loop (seq,
                                         ARRAY_ACCESS (this->hairpin_loop, k).i,
                                         ARRAY_ACCESS (this->hairpin_loop, k).j,
                                      ARRAY_ACCESS (this->hairpin_loop, k).size,
                                          scores);
   }

   /* multiloops */
   for (k = 0; k < ARRAY_CURRENT (this->multi_loop); k++)
   {
      Gm += nn_scores_get_G_extloop_multiloop (seq,
                                    ARRAY_ACCESS (this->multi_loop, k).unpaired,
                                    ARRAY_ACCESS (this->multi_loop, k).nstems,
                                    ARRAY_ACCESS (this->multi_loop, k).stems,
                                    ARRAY_ACCESS (this->multi_loop, k).ndangle5,
                                    ARRAY_ACCESS (this->multi_loop, k).dangle5,
                                    ARRAY_ACCESS (this->multi_loop, k).ndangle3,
                                    ARRAY_ACCESS (this->multi_loop, k).dangle3,
                                               true,
                                               scores);
   }

/*    mfprintf (stdout, "external: %d\n" */
/*                      "stack:    %d\n" */
/*                      "bulge:    %d\n" */
/*                      "hairpin:  %d\n" */
/*                      "internal: %d\n" */
/*                      "multi:    %d\n", */
/*              Ge, Gs, Gb, Gh, Gi, Gm); */

   G = Ge + Gs + Gb + Gh + Gi + Gm;

   return G;
}


/*********************************   Output   *********************************/

/** @brief Print a certain stack of a secondary structure to a stream.
 *
 * Format is "k - l". The trick is just, that "k" and "l" have the same column
 * width.\n
 * @params[in] stream Output stream to write to.
 * @params[in] i Stacking pair to print.
 * @params[in] this secondary structure.
 */
void
secstruct_fprintf_i_stack (FILE* stream,
                           const unsigned long i,
                           const SecStruct* this)
{
   unsigned long p5, p3;

   assert (this != NULL);
   assert (ARRAY_NOT_NULL (this->stack));
   assert (ARRAY_CURRENT (this->stack) > i);

   p5 = ARRAY_ACCESS (this->stack, i).i;
   p3 = ARRAY_ACCESS (this->stack, i).j;

   mfprintf (stream, "%lu - %lu", p5, p3);
}

/** @brief Print the list of stacks of a secondary structure to a stream.
 *
 * Format is "index of pair: i - j".\n
 * @params[in] stream Output stream to write to.
 * @params[in] this secondary structure.
 */
void
secstruct_fprintf_stacks (FILE* stream, const SecStruct* this)
{
   unsigned long i;
   unsigned long tmp;
   int rprec, rpreci;
   unsigned long pline_width = 2;
   char* string;
   char* string_start;

   assert (this != NULL);
   assert (ARRAY_NOT_NULL (this->stack));

   /* find largest pair */
   for (i = 0; i < ARRAY_CURRENT (this->stack); i++)
   {
      rprec = 0;
      tmp = ARRAY_ACCESS (this->stack, i).i;
      if (tmp < ARRAY_ACCESS (this->stack, i).j)
      {
         tmp = ARRAY_ACCESS (this->stack, i).j;
      }

      /* get no. of digits */
      if (tmp > 0)
      {
         rprec += floor (log10 (tmp) + 1);
      }
      else
      {
         rprec += 1;
      }
      
      if ((unsigned) rprec > pline_width) 
      {
         pline_width = rprec;
      }
   }
   rprec = pline_width;   
   if (i > 0)
   {
      rpreci = floor (log10 (i) + 1);
   }
   else
   {
      rpreci = 1;
   }

   /* add up components of a line */
   pline_width += rprec;
   pline_width += rpreci;
   pline_width += 5;            /*:\s\s-\s*/
   pline_width += 1;            /* + \n */

   /* allocate buffer */
   string = (char*) XMALLOC (sizeof (*string)
                            * ((pline_width * ARRAY_CURRENT (this->stack))+ 1));
   if (string == NULL)
   {
      return;
   }
   string_start = string;  

   /* print */
   for (i = 0; i < ARRAY_CURRENT (this->stack); i++)
   {
      msprintf (string, "%*lu: %*lu - %*lu\n", rpreci, i, rprec,
                ARRAY_ACCESS (this->stack, i).i, rprec,
                ARRAY_ACCESS (this->stack, i).j);
 
      string += pline_width;      

   }

   string[0] = '\0';
   mfprintf (stream, "%s", string_start);  

   XFREE (string_start);
}

/** @brief Print a certain hairpin loop of a secondary structure to a stream.
 *
 * Format is "k - l (size)". The trick is just, that "k" and "l" have the same
 * column width.\n
 * @params[in] stream Output stream to write to.
 * @params[in] i Hairpin loop to print.
 * @params[in] this secondary structure.
 */
void
secstruct_fprintf_i_hairpin (FILE* stream,
                             const unsigned long i,
                             const SecStruct* this)
{
   unsigned long start, end, size;

   assert (this != NULL);
   assert (ARRAY_NOT_NULL (this->hairpin_loop));
   assert (ARRAY_CURRENT (this->hairpin_loop) > i);

   start = ARRAY_ACCESS (this->hairpin_loop, i).i;
   end   = ARRAY_ACCESS (this->hairpin_loop, i).j;
   size  = ARRAY_ACCESS (this->hairpin_loop, i).size;

   mfprintf (stream, "%lu - %lu (%lu)", start, end, size);
}

/** @brief Print the list of hairpins of a secondary structure to a stream.
 *
 * Format is "index of pair: i - j (size)".\n
 * @params[in] stream Output stream to write to.
 * @params[in] this secondary structure.
 */
void
secstruct_fprintf_hairpins (FILE* stream, const SecStruct* this)
{
   unsigned long i;
   unsigned long tmp, tmps;
   int rprec, rpreci, rprecs = 0;
   unsigned long pline_width = 2;
   int size_width = 0;
   char* string;
   char* string_start;

   assert (this != NULL);
   assert (ARRAY_NOT_NULL (this->hairpin_loop));

   /* find largest no. */
   for (i = 0; i < ARRAY_CURRENT (this->hairpin_loop); i++)
   {
      rprec = 0;
      rprecs = 0;
     
      tmp = ARRAY_ACCESS (this->hairpin_loop, i).i;
      if (tmp < ARRAY_ACCESS (this->hairpin_loop, i).j)
      {
         tmp = ARRAY_ACCESS (this->hairpin_loop, i).j;
      }

      /* get no. of digits */
      if (tmp > 0)
      {
         rprec += floor (log10 (tmp) + 1);
      }
      else
      {
         rprec += 1;
      }
 
      tmps = ARRAY_ACCESS (this->hairpin_loop, i).size;     
      if (tmps > 0)
      {
         rprecs += floor (log10 (tmps) + 1);
      }
      else
      {
         rprecs += 1;
      }
      if (rprecs > size_width)
      {
         size_width = rprecs;
      }
      
      if ((unsigned) rprec > pline_width) 
      {
         pline_width = rprec;
      }
   }
   rprec = pline_width;
   if (i > 0)
   {
      rpreci = floor (log10 (i) + 1);
   }
   else
   {
      rpreci = 1;
   }

   /* add up components of a line */
   pline_width += rprec;
   pline_width += size_width;
   pline_width += rpreci;
   pline_width += 8;            /*:\s\s-\s\s()*/
   pline_width += 1;            /* + \n */

   /* allocate buffer */
   string = (char*) XMALLOC (sizeof (*string) * ((pline_width *
                                      ARRAY_CURRENT (this->hairpin_loop)) + 1));
   if (string == NULL)
   {
      return;
   }
   string_start = string;  

   /* print */
   for (i = 0; i < ARRAY_CURRENT (this->hairpin_loop); i++)
   {
      msprintf (string, "%*lu: %*lu - %*lu (%*lu)\n", rpreci, i, rprec,
                ARRAY_ACCESS (this->hairpin_loop, i).i, rprec,
                ARRAY_ACCESS (this->hairpin_loop, i).j, size_width,
                ARRAY_ACCESS (this->hairpin_loop, i).size);
 
      string += pline_width;      
   }

   string[0] = '\0';
   mfprintf (stream, "%s", string_start);  

   XFREE (string_start);
}

/** @brief Print a certain bulge loop of a secondary structure to a stream.
 *
 * Format is "i/j -  k/l (size)".\n
 * @params[in] stream Output stream to write to.
 * @params[in] i Bulge loop to print.
 * @params[in] this secondary structure.
 */
void
secstruct_fprintf_i_bulge (FILE* stream,
                           const unsigned long i,
                           const SecStruct* this)
{
   unsigned long pi1, pj1, pi2, pj2, size;

   assert (this != NULL);
   assert (ARRAY_NOT_NULL (this->bulge_loop));
   assert (ARRAY_CURRENT (this->bulge_loop) > i);

   pi1 = ARRAY_ACCESS(this->bulge_loop, i).i1;
   pj1 = ARRAY_ACCESS(this->bulge_loop, i).j1;
   pi2 = ARRAY_ACCESS(this->bulge_loop, i).i2;
   pj2 = ARRAY_ACCESS(this->bulge_loop, i).j2;
   size = ARRAY_ACCESS(this->bulge_loop, i).size;

   mfprintf (stream, "%lu/%lu - %lu/%lu (%lu)", pi1, pj1, pi2, pj2, size);
}

/** @brief Print the list of bulges of a secondary structure to a stream.
 *
 * Format is "index of pair: i_1/j_1 - i_2/j_2 (size)".\n
 * @params[in] stream Output stream to write to.
 * @params[in] this secondary structure.
 */
void
secstruct_fprintf_bulges (FILE* stream, const SecStruct* this)
{
   unsigned long i;
   unsigned long tmp, tmps;
   int rprec, rprecs = 0, rpreci;
   unsigned long pline_width = 2;
   int size_width = 0;
   char* string;
   char* string_start;

   assert (this != NULL);
   assert (ARRAY_NOT_NULL (this->bulge_loop));

   /* find largest no. */
   for (i = 0; i < ARRAY_CURRENT (this->bulge_loop); i++)
   {
      rprec = 0;
      rprecs = 0;

      /* width of 'size' component */
      tmps = ARRAY_ACCESS (this->bulge_loop, i).size;

      if (tmps > 0)
      {
         rprecs += floor (log10 (tmps) + 1);
      }
      else
      {
         rprecs += 1;
      }
      if (rprecs > size_width)
      {
         size_width = rprecs;
      }

      /* base pos. component */
      tmp = ARRAY_ACCESS (this->bulge_loop, i).i1;
      if (ARRAY_ACCESS (this->bulge_loop, i).j1 > tmp)
      {
         tmp = ARRAY_ACCESS (this->bulge_loop, i).j1;
      }
      if (ARRAY_ACCESS (this->bulge_loop, i).i2 > tmp)
      {
         tmp = ARRAY_ACCESS (this->bulge_loop, i).i2;
      }
      if (ARRAY_ACCESS (this->bulge_loop, i).j2 > tmp)
      {
         tmp = ARRAY_ACCESS (this->bulge_loop, i).j2;
      }  

      if (tmp > 0)
      {
         rprec += floor (log10 (tmp) + 1);
      }
      else
      {
         rprec += 1;
      }
      
      if ((unsigned) rprec > pline_width) 
      {
         pline_width = rprec;
      }
   }
   rprec = pline_width;
   if (i > 0)
   {
      rpreci = floor (log10 (i) + 1);
   }
   else
   {
      rpreci = 1;
   }

   /* add up components of a line */
   pline_width = (rprec * 4);
   pline_width += rpreci;
   pline_width += size_width;
   pline_width += 10;            /*:\s/\s-\s/\s()*/
   pline_width += 1;            /* + \n */

   /* allocate buffer */
   string = (char*) XMALLOC (sizeof (*string) * ((pline_width *
                                        ARRAY_CURRENT (this->bulge_loop)) + 1));
   if (string == NULL)
   {
      return;
   }
   string_start = string;

   /* print */
   for (i = 0; i < ARRAY_CURRENT (this->bulge_loop); i++)
   {
      msprintf (string, "%*lu: %*lu/%*lu - %*lu/%*lu (%*lu)\n",
                rpreci, i,
                rprec, ARRAY_ACCESS (this->bulge_loop, i).i1,
                rprec, ARRAY_ACCESS (this->bulge_loop, i).j1,
                rprec, ARRAY_ACCESS (this->bulge_loop, i).i2,
                rprec, ARRAY_ACCESS (this->bulge_loop, i).j2,
                size_width, ARRAY_ACCESS (this->bulge_loop, i).size);
 
      string += pline_width;
   }

   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   XFREE (string_start);
}

/** @brief Print a certain internal loop of a secondary structure to a stream.
 *
 * Format is "i/j -  k/l (size/size)".\n
 * @params[in] stream Output stream to write to.
 * @params[in] i Index of internal loop to print.
 * @params[in] this secondary structure.
 */
void
secstruct_fprintf_i_internal (FILE* stream,
                              const unsigned long i,
                              const SecStruct* this)
{
   unsigned long pi1, pj1, pi2, pj2, size1, size2;

   assert (this != NULL);
   assert (ARRAY_NOT_NULL (this->internal_loop));
   assert (ARRAY_CURRENT (this->internal_loop) > i);

   pi1   = ARRAY_ACCESS(this->internal_loop, i).i1;
   pj1   = ARRAY_ACCESS(this->internal_loop, i).j1;
   pi2   = ARRAY_ACCESS(this->internal_loop, i).i2;
   pj2   = ARRAY_ACCESS(this->internal_loop, i).j2;
   size1 = ARRAY_ACCESS(this->internal_loop, i).size1;
   size2 = ARRAY_ACCESS(this->internal_loop, i).size2;

   mfprintf (stream, "%lu/%lu - %lu/%lu (%lu/%lu)",
             pi1, pj1, pi2, pj2, size1, size2);
}

/** @brief Print a list of internal loops of a secondary structure to a stream.
 *
 * Format is "index of pair: i_1/j_1 - i_2/j_2 (size1/size2)".\n
 * @params[in] stream Output stream to write to.
 * @params[in] this secondary structure.
 */
void
secstruct_fprintf_internals (FILE* stream, const SecStruct* this)
{
   unsigned long i;
   unsigned tmp, tmps;
   int rprec, rprecs = 0, rpreci;
   unsigned long pline_width = 2;
   int size_width = 0;
   char* string;
   char* string_start;

   assert (this);
   assert (ARRAY_NOT_NULL (this->internal_loop));

   /* find largest no. */
   for (i = 0; i < ARRAY_CURRENT (this->internal_loop); i++)
   {
      rprec = 0;
      rprecs = 0;

      /* width of 'size' component */
      tmps = ARRAY_ACCESS (this->internal_loop, i).size1;

      if ((  ARRAY_ACCESS (this->internal_loop, i).size2
             * ARRAY_ACCESS (this->internal_loop, i).size2)
          > (unsigned)(tmps * tmps))
      {
         tmps = ARRAY_ACCESS (this->internal_loop, i).size2;
      }

      if (tmps > 0)
      {
         rprecs += floor (log10 (tmps) + 1);
      }
      else
      {
         rprecs += 1;
      }
      if (rprecs > size_width)
      {
         size_width = rprecs;
      }

      /* base pos. component */
      tmp = ARRAY_ACCESS (this->internal_loop, i).i1;
      if (ARRAY_ACCESS (this->internal_loop, i).j1 > tmp)
      {
         tmp = ARRAY_ACCESS (this->internal_loop, i).j1;
      }
      if (ARRAY_ACCESS (this->internal_loop, i).i2 > tmp)
      {
         tmp = ARRAY_ACCESS (this->internal_loop, i).i2;
      }
      if (ARRAY_ACCESS (this->internal_loop, i).j2 > tmp)
      {
         tmp = ARRAY_ACCESS (this->internal_loop, i).j2;
      }  

      if (tmp > 0)
      {
         rprec += floor (log10 (tmp) + 1);
      }
      else
      {
         rprec += 1;
      }
      
      if ((unsigned) rprec > pline_width) 
      {
         pline_width = rprec;
      }
   }
   rprec = pline_width;
   if (i > 0)
   {
      rpreci = floor (log10 (i) + 1);
   }
   else
   {
      rpreci = 1;
   }

   /* add up components of a line */
   pline_width = (rprec * 4);
   pline_width += rpreci;
   pline_width += (2 * size_width);
   pline_width += 11;            /*:\s/\s-\s/\s(/)*/
   pline_width += 1;            /* + \n */

   /* allocate buffer */
   string = (char*) XMALLOC (sizeof (*string) * ((pline_width *
                                     ARRAY_CURRENT (this->internal_loop)) + 1));
   if (string == NULL)
   {
      return;
   }
   string_start = string;


   /* print */
   for (i = 0; i < ARRAY_CURRENT (this->internal_loop); i++)
   {
      msprintf (string, "%*lu: %*lu/%*lu - %*lu/%*lu (%*lu/%*lu)\n",
                rpreci, i,
                rprec, ARRAY_ACCESS (this->internal_loop, i).i1,
                rprec, ARRAY_ACCESS (this->internal_loop, i).j1,
                rprec, ARRAY_ACCESS (this->internal_loop, i).i2,
                rprec, ARRAY_ACCESS (this->internal_loop, i).j2,
                size_width, ARRAY_ACCESS (this->internal_loop, i).size1,
                size_width, ARRAY_ACCESS (this->internal_loop, i).size2);
 
      string += pline_width;
   }

   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   XFREE (string_start);
}

#define SENT1 "Unpaired bases: "
#define DILM1 " - "
#define SENT2 "Stems: "
#define SENT3 "5' dangling ends: "
#define SENT4 "3' dangling ends: "
#define INDT1 "  "
static size_t
printsize_of_multiloop (const MultiLoop* ml)
{
   unsigned long i, j;
   unsigned long tmp;
   size_t store_size = 0;

   /* storage for "Unpaired bases: \d+\n" */
   store_size += strlen (SENT1);
   store_size += 3;             /* "  " + \n */
   if (ml->unpaired > 0)
   {
      store_size += floor (log10 (ml->unpaired) + 1);
   }
   else
   {
      store_size += 1;
   }

   /* storage for stems: \d - \d */
   store_size += strlen (SENT2);
   store_size += 3;             /* "  " + \n */
   if (ml->nstems == 0)
   {
      store_size += 1;
   }
   else
   {
      store_size += floor (log10 (ml->nstems) + 1);
      tmp = 0;
      for (i = 0; i < ml->nstems; i++)
      {
         for (j = 0; j < No_Of_Strands; j++)
         {
            if (ml->stems[i][j] > tmp)
            {
               tmp =  ml->stems[i][j];
            }
         }
      }
      if (tmp > 0)
      {
         tmp = floor (log10 (tmp) + 1);
      }
      else
      {
         tmp = 1;
      }
      store_size += ((floor (log10 (i) + 1) + 4) * ml->nstems);

      store_size += (strlen (INDT1) * ml->nstems);
      store_size += (tmp * No_Of_Strands * ml->nstems);
      store_size += (strlen (DILM1) * (No_Of_Strands - 1) * ml->nstems);
      store_size += ml->nstems;
   }

   /* storage for 5' dangling end */
   store_size += strlen (SENT3);
   store_size += 3;             /* "  " + \n */
   if (ml->ndangle5 == 0)
   {
      store_size += 1;
   }
   else
   {
      store_size += floor (log10 (ml->ndangle5) + 1);
      tmp = 0;
      for (i = 0; i < ml->ndangle5; i++)
      {
         for (j = 0; j < No_Of_Dangles; j++)
         {
            if (ml->dangle5[i][j] > tmp)
            {
               tmp =  ml->dangle5[i][j];
            }
         }
      }
      if (tmp > 0)
      {
         tmp = floor (log10 (tmp) + 1);
      }
      else
      {
         tmp = 1;
      }
      store_size += ((floor (log10 (i) + 1) + 4) * ml->ndangle5);

      store_size += (strlen (INDT1) * ml->ndangle5);
      store_size += (tmp * No_Of_Dangles * ml->ndangle5);
      store_size += (strlen (DILM1) * (No_Of_Dangles - 1) * ml->ndangle5);
      store_size += ml->ndangle5;
   }

   /* storage for 3' dangling end */
   store_size += strlen (SENT4);
   store_size += 3;             /* "  " + \n */
   if (ml->ndangle3 == 0)
   {
      store_size += 1;
   }
   else
   {
      store_size += floor (log10 (ml->ndangle3) + 1);
      tmp = 0;
      for (i = 0; i < ml->ndangle3; i++)
      {
         for (j = 0; j < No_Of_Dangles; j++)
         {
            if (ml->dangle3[i][j] > tmp)
            {
               tmp =  ml->dangle3[i][j];
            }
         }
      }
      if (tmp > 0)
      {
         tmp = floor (log10 (tmp) + 1);
      }
      else
      {
         tmp = 1;
      }
      store_size += ((floor (log10 (i) + 1) + 4) * ml->ndangle3);

      store_size += (strlen (INDT1) * ml->ndangle3);
      store_size += (tmp * No_Of_Dangles * ml->ndangle3);
      store_size += (strlen (DILM1) * (No_Of_Dangles - 1) * ml->ndangle3);
      store_size += ml->ndangle3;
   }

   return store_size;
}

static size_t
sprintf_multiloop (char* str, const MultiLoop* ml)
{
   int prec, preci;
   unsigned long i, j;
   char* string_start = str;

   /* Unpaired bases */
   msprintf (str, "  %s%lu\n", SENT1, ml->unpaired);
   str += strlen (SENT1) + 3;
   if (ml->unpaired > 0)
   {
      prec = floor (log10 (ml->unpaired) + 1);     
   }
   else
   {
      prec = 1;
   }
   str += prec;

   /* stems */
   msprintf (str, "  %s%lu\n", SENT2, ml->nstems);
   str += strlen (SENT2) + 3;
   if (ml->nstems == 0)
   {
      str += 1;     
   }
   else
   {
      prec = floor (log10 (ml->nstems) + 1);
      str += prec;

      /* calc didgits for stems */
      prec = 0;
      for (i = 0; i < ml->nstems; i++)
      {
         for (j = 0; j < No_Of_Strands; j++)
         {
            if (ml->stems[i][j] > (unsigned) prec)
            {
               prec =  ml->stems[i][j];
            }
         }
      }
      if (prec > 0)
      {
         prec = floor (log10 (prec) + 1);
      }
      else
      {
         prec = 1;
      }
      preci = floor (log10 (i) + 1);

      /* print stems */   
      for (i = 0; i < ml->nstems; i++)
      {
         msprintf (str, "  %s%*lu: ", INDT1, preci, i);
         str += (strlen (INDT1) + preci + 4);
         for (j = 0; j < No_Of_Strands; j++)
         {
            if (j == 0)
            {
               msprintf (str, "%*lu", prec, ml->stems[i][j]);
               str += prec;
            }
            else
            {
               msprintf (str, "%s%*lu", DILM1, prec, ml->stems[i][j]);
               str += (strlen (DILM1) + prec);
            }
         }
         msprintf (str, "\n");
         str += 1;
      }
   }

   /* 5' dangle */
   msprintf (str, "  %s%lu\n", SENT3, ml->ndangle5);
   str += strlen (SENT3) + 3;
   if (ml->ndangle5 == 0)
   {
      str += 1;
   }
   else
   {
      prec = floor (log10 (ml->ndangle5) + 1);
      str += prec;
      
      /* calc didgits for 5' dangles */
      prec = 0;
      for (i = 0; i < ml->ndangle5; i++)
      {
         for (j = 0; j < No_Of_Dangles; j++)
         {
            if (ml->dangle5[i][j] > (unsigned) prec)
            {
               prec =  ml->dangle5[i][j];
            }
         }
      }
      if (prec > 0)
      {
         prec = floor (log10 (prec) + 1);
      }
      else
      {
         prec = 1;
      }
      preci = floor (log10 (i) + 1);

      for (i = 0; i < ml->ndangle5; i++)
      {
         msprintf (str, "  %s%*lu: ", INDT1, preci, i);
         str += (strlen (INDT1) + preci + 4);
         for (j = 0; j < No_Of_Dangles; j++)
         {
            if (j == 0)
            {
               msprintf (str, "%*lu", prec, ml->dangle5[i][j]);
               str += prec;
            }
            else
            {
               msprintf (str, "%s%*lu", DILM1, prec, ml->dangle5[i][j]);
               str += (strlen (DILM1) + prec);
            }
         }
         msprintf (str, "\n");
         str += 1;
      }
   }

  /* 3' dangle */
   msprintf (str, "  %s%lu\n", SENT4, ml->ndangle3);
   str += strlen (SENT4) + 3;
   if (ml->ndangle3 == 0)
   {
      str += 1;
   }
   else
   {
      prec = floor (log10 (ml->ndangle3) + 1);
      str += prec;
      
      /* calc didgits for 3' dangles */
      prec = 0;
      for (i = 0; i < ml->ndangle3; i++)
      {
         for (j = 0; j < No_Of_Dangles; j++)
         {
            if (ml->dangle3[i][j] > (unsigned) prec)
            {
               prec =  ml->dangle3[i][j];
            }
         }
      }
      if (prec > 0)
      {
         prec = floor (log10 (prec) + 1);
      }
      else
      {
         prec = 1;
      }
      preci = floor (log10 (i) + 1);

      for (i = 0; i < ml->ndangle3; i++)
      {
         msprintf (str, "  %s%*lu: ", INDT1, preci, i);
         str += (strlen (INDT1) + preci + 4);
         for (j = 0; j < No_Of_Dangles; j++)
         {
            if (j == 0)
            {
               msprintf (str, "%*lu", prec, ml->dangle3[i][j]);
               str += prec;
            }
            else
            {
               msprintf (str, "%s%*lu", DILM1, prec, ml->dangle3[i][j]);
               str += (strlen (DILM1) + prec);
            }
         }
         msprintf (str, "\n");
         str += 1;
      }
   }

   str[0] = '\0';

   return str - string_start;
}

void
secstruct_fprintf_external (FILE* stream, const SecStruct* this)
{
   char* string;
   size_t loop_storage;

   assert (this);

   loop_storage = printsize_of_multiloop (&this->ext_loop);

   string = (char*) XMALLOC (sizeof (*string) * (loop_storage + 1));
   if (string == NULL)
   {
      return;
   }

   sprintf_multiloop (string, &this->ext_loop);

   mfprintf (stream, "%s", string);

   XFREE (string);
}

/** @brief Print a certain multiloop of a secondary structure to a stream.
 *
 * Print a certain multiloop of a secondary structure to a stream..\n
 * @params[in] stream Output stream to write to.
 * @params[in] i Index of loop to print.
 * @params[in] this secondary structure.
 */
void
secstruct_fprintf_i_multiloop (FILE* stream,
                               const unsigned long i,
                               const SecStruct* this)
{
   char* string;
   size_t loop_storage;

   assert (this != NULL);
   assert (ARRAY_NOT_NULL (this->multi_loop));
   assert (ARRAY_CURRENT (this->multi_loop) > i);

   loop_storage = printsize_of_multiloop (&(ARRAY_ACCESS(this->multi_loop, i)));

   string = (char*) XMALLOC (sizeof (*string) * (loop_storage + 1));
   if (string == NULL)
   {
      return;
   }

   sprintf_multiloop (string,
                      &(ARRAY_ACCESS (this->multi_loop, i)));

   mfprintf (stream, "%s", string);

   XFREE (string);
}

void
secstruct_fprintf_multiloops (FILE* stream, const SecStruct* this)
{
   unsigned long i;
   size_t storage = 0;
   int preci;
   char* string;
   char* string_start;

   /* calc storage of loops */
   for (i = 0; i < ARRAY_CURRENT (this->multi_loop); i++)
   {
      storage += printsize_of_multiloop (&(ARRAY_ACCESS (this->multi_loop, i)));
   }
   if (i > 0)
   {
      preci = floor (log10 (i) + 1);
   }
   else
   {
      preci = 1;
   }

   /* add storage for "\d+:\n" */
   storage += ((preci + 2) * ARRAY_CURRENT (this->multi_loop));

   string = (char*) XMALLOC (sizeof (*string) * (storage + 1));
   if (string == NULL)
   {
      return;
   }
   string_start = string;

   /* print */
   for (i = 0; i < ARRAY_CURRENT (this->multi_loop); i++)
   {
      msprintf (string, "%*lu:\n", preci, i);
      string += preci + 2;
      storage = sprintf_multiloop (string,
                                   &(ARRAY_ACCESS (this->multi_loop, i)));
      string += storage;
   }

   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   XFREE (string_start);
}

void
secstruct_fprintf_seqpos_map (FILE* stream, const SecStruct* this)
{
   unsigned long i;
   unsigned int seq_idx_prec, elem_idx_prec;
   char hairpin[] = "Hairpin loop";
   char stack[] = "Stacked base pair";
   char bulge[] = "Bulge loop";
   char internal[] = "Internal loop";
   char multi[] = "Multi loop";
   char external[] = "External loop";
   char mto[] = "-> ";
   char none[] = "Not assigned";
   char panic[] = "(o.O)";
   unsigned long panic_size = strlen (panic);
   char* type_id[N_FEATURES];
   int   type_sz[N_FEATURES];
   char* string;
   char* string_start;
   unsigned long storage = 0;
   unsigned long no_of_chars;

   assert (stream);
   assert (this);
   assert ((this->seq_struct_map) || (this->size == 0));

   type_id[SCSTRCT_HAIRPIN] = hairpin;
   type_id[SCSTRCT_STACK] = stack;
   type_id[SCSTRCT_BULGE] = bulge;
   type_id[SCSTRCT_INTERNAL] = internal;
   type_id[SCSTRCT_MULTI] = multi;
   type_id[SCSTRCT_EXTERNAL] = external;
   type_id[SCSTRCT_MTO] = mto;
   type_id[SCSTRCT_VOID] = none;

   for (i = 0; i < N_FEATURES; i++)
   {
      type_sz[i] = strlen(type_id[i]);
   }

   /* get largest no.s for field widths */
   elem_idx_prec = 0;
   for (i = 0; i < this->size; i++)
   {
      if (elem_idx_prec < this->seq_struct_map[i].idx)
      {
         elem_idx_prec = this->seq_struct_map[i].idx;
      }

      if (this->seq_struct_map[i].type < N_FEATURES)
      {
         storage += type_sz[this->seq_struct_map[i].type];

         if (this->seq_struct_map[i].type == SCSTRCT_MTO)
         {
            storage += 7;
            if ( ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).idx
                 > 0)
            {
               storage += floor (log10 (
                ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).idx)
                                 + 1);
            }
            else
            {
               storage++;
            }
            storage +=
   type_sz[ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).type];

            if (ARRAY_ACCESS(this->multi_map,this->seq_struct_map[i].idx+1).idx
                 > 0)
            {
               storage += floor (log10 (
             ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx + 1).idx)
                                 + 1);
            }
            else
            {
               storage++;
            }
            storage +=
   type_sz[ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx + 1).type];
         }
      }
      else
      {
         THROW_WARN_MSG ("Unrecognised structure element type found at "
                         "sequence position %lu.", i);
         storage += panic_size;
      }
   }

   /* calc field widths */ 
   if (this->size > 0)
   {
      seq_idx_prec = floor (log10 (this->size - 1) + 1);
   }
   else
   {
      seq_idx_prec = 1;
   }
   if (elem_idx_prec > 0)
   {
      elem_idx_prec = floor (log10 (elem_idx_prec) + 1);
   }
   else
   {
      elem_idx_prec = 1;
   }

   /* 2b printed: list of seq pos, each on its own line
      tupel idx with elemidx + type
      align seqidx, elem idx and type
      multi: write list in line */
   /* alloc: seq_idx_prec + elem_idx_prec + type_prec + ":  " + "   " */
   no_of_chars  = seq_idx_prec;
   no_of_chars += 7;                /* ":  " + "   " + "\n" */
   no_of_chars += elem_idx_prec;
   no_of_chars *= this->size;
   no_of_chars += storage;
   string = (char*) XMALLOC (sizeof (*string) * (no_of_chars + 1));
   if (string == NULL)
   {
      return;
   }
   string_start = string;

   /*mprintf ("Allocated: %lu\n", sizeof (*string) * (no_of_chars + 1));*/

   for (i = 0; i < this->size; i++)
   {
      msprintf (string, "%*lu:  ", seq_idx_prec, i);
      string += seq_idx_prec + 3;

      if (this->seq_struct_map[i].type < N_FEATURES)
      {
         msprintf(string,"%*lu   ",elem_idx_prec,this->seq_struct_map[i].idx);
         string += elem_idx_prec + 3;
         msprintf (string, "%s", type_id[this->seq_struct_map[i].type]);
         string += type_sz[this->seq_struct_map[i].type];

         if (this->seq_struct_map[i].type == SCSTRCT_MTO)
         {
           msprintf (string, "%lu  %s; ", 
                 ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).idx,
      type_id[ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).type]);
           string += 4 +
       type_sz[ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).type];
           if ( ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).idx
                > 0)
           {
              storage = floor (log10 (
                 ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx).idx)
                                + 1);
              string += storage;
           }
           else
           {
              string++;
           }

           msprintf (string, "%lu  %s;",
             ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx + 1).idx,
 type_id[ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx + 1).type]);
                  string += 3 +
   type_sz[ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx + 1).type];
           if (ARRAY_ACCESS(this->multi_map,this->seq_struct_map[i].idx + 1).idx
                > 0)
           {
              storage = floor (log10 (
             ARRAY_ACCESS(this->multi_map, this->seq_struct_map[i].idx + 1).idx)
                                + 1);
              string += storage;
           }
           else
           {
              string++;
           }
         }
      }
      else
      {
         msprintf (string, "%*c   ", elem_idx_prec, '?');
         string += elem_idx_prec + 3;
         msprintf (string, "%s",  panic);
         string += panic_size;
      }

      string[0] = '\n';
      string++;
   }

   string[0] = '\0';
   mfprintf (stream, "%s", string_start);

   /*mprintf ("Written: %d\n", string - string_start);*/

   XFREE (string_start);
}
