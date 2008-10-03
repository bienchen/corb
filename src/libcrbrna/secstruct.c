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
 *
 */


#include <config.h>
#include <stdlib.h> /* only as long as EXIT_FAILURE is in use */
#include <libcrbbasic/crbbasic.h>
#include "secstruct.h"

/* structure to store hairpins */
typedef struct {
      unsigned long i;          /* 5' end of loop */
      unsigned long j;          /* 3' end of loop */
      unsigned long size;       /* loop size */
} HairpinLoop;

ARRAY_CREATE_CLASS(HairpinLoop);

/* structure to store stacks */
typedef struct {
/* TODO: only i1, j1 needed (can then be called i and j) */
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

/* indeces for multiloop arrays */
enum stem_indeces {
   P5_Strand = 0,
   P3_Strand,
   No_Of_Strands
};

enum dangle_indeces {
   P5_Dangle = 0,
   P3_Dangle,
   Ne_Dangle,
   No_Of_Dangles
};

/* store multiloops */
typedef struct {
      unsigned long unpaired;
      unsigned long nstems, ndangle5, ndangle3;
      unsigned long (*stems)[No_Of_Strands];
      unsigned long (*dangle5)[No_Of_Dangles];
      unsigned long (*dangle3)[No_Of_Dangles];
} MultiLoop;

ARRAY_CREATE_CLASS(MultiLoop);

/* data structure to store decomposed secondary structures */
struct SecStruct {
      ArrayUlong tetra_loop;    /* store end position of a tetra loop */
      ArrayHairpinLoop hairpin_loop;
      ArrayStackLoop stack;
      ArrayBulgeLoop bulge_loop;
      ArrayIntLoop internal_loop;
      ArrayMultiLoop multi_loop;
      MultiLoop ext_loop;
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
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
SecStruct*
secstruct_new (const char* file, const int line)
{
   SecStruct* this = XOBJ_MALLOC(sizeof (SecStruct), file, line);
   
   if (this != NULL)
   {
      ARRAY_SET_VOID (this->tetra_loop);
      ARRAY_SET_VOID (this->hairpin_loop);
      ARRAY_SET_VOID (this->stack);
      ARRAY_SET_VOID (this->bulge_loop);
      ARRAY_SET_VOID (this->internal_loop);
      ARRAY_SET_VOID (this->multi_loop);

      ARRAY_ULONG_INIT (this->tetra_loop, 1);
      if (ARRAY_IS_NULL (this->tetra_loop))
      {
         secstruct_delete (this);
         return NULL;
      }

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
      ARRAY_DELETE (this->tetra_loop);
      ARRAY_DELETE (this->hairpin_loop);
      ARRAY_DELETE (this->stack);
      ARRAY_DELETE (this->bulge_loop);
      ARRAY_DELETE (this->internal_loop);

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

static void
multiloop_find (unsigned long i, unsigned long j, const unsigned long* pairs,
                unsigned long size, MultiLoop* ml)
{
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

         if (i > 0)
         {
            /* add dangle5 */
            ml->dangle5[ml->ndangle5][P5_Dangle] = i;
            ml->dangle5[ml->ndangle5][P3_Dangle] = pairs[i];
            ml->dangle5[ml->ndangle5][Ne_Dangle] = i - 1;
            ml->ndangle5++;
         }
        
         if (pairs[i] < size - 1)
         {
            /* add dangle3 */
            ml->dangle3[ml->ndangle3][P5_Dangle] = i;
            ml->dangle3[ml->ndangle3][P3_Dangle] = pairs[i];
            ml->dangle3[ml->ndangle3][Ne_Dangle] = pairs[i] + 1;
            ml->ndangle3++;
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

   /*mfprintf (stderr, "UNPAIRED: %lu\n", ml->unpaired);

   mfprintf (stderr, "NSTEMS: %lu\n", ml->nstems);
   for (i = 0; i < ml->nstems; i++)
   {
      mfprintf (stderr, "  5': %2lu 3': %2lu\n",
                ml->stems[i][P5_Strand],
                ml->stems[i][P3_Strand]);
   }

   mfprintf (stderr, "NDANGLE5: %lu\n", ml->ndangle5);
   for (i = 0; i < ml->ndangle5; i++)
   {
      mfprintf (stderr, "  5': %2lu 3': %2lu D: %2lu\n",
                ml->dangle5[i][P5_Dangle],
                ml->dangle5[i][P3_Dangle],
                ml->dangle5[i][Ne_Dangle]);      
   }

   mfprintf (stderr, "Ndangle3: %lu\n", ml->ndangle3); 
   for (i = 0; i < ml->ndangle3; i++)
   {
      mfprintf (stderr, "  5': %2lu 3': %2lu D: %2lu\n",
                ml->dangle3[i][P5_Dangle],
                ml->dangle3[i][P3_Dangle],
                ml->dangle3[i][Ne_Dangle]);      
                }*/
}

/** @brief Decompose an RNA structure stored.
 *
 * Analyses a RNA structure and stores information about its components in the
 * SecStruct object.\n
 * Returns 0 on success, ... else.
 *
 * @param[in] pairs List of pairing partners, where i pairs pairs[i], an
                    unpaired position k is indicated by pairs[k] = NOT_PAIRED
 * @param[in] size Length of the structure
 * @param[in] this SecStruct data object.
 */
int
secstruct_find_interactions (const unsigned long* pairs,
                             const unsigned long size,
                             SecStruct* this)
{
   unsigned long i, p, q, size1, size2, i1, j1, i2, j2;
   unsigned long remains = size;
   int error = 0;

   assert (this);
   /* exterior loop */
   error = multiloop_allocate (remains, &this->ext_loop);
   if (size > 0)
   {
      multiloop_find (0, (size - 1), pairs, size, &this->ext_loop);
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
            ARRAY_PUSH(this->hairpin_loop, hl, HairpinLoop, { error = 1; });
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
               i1 = i;
               j1 = pairs[i];
               i2 = p;
               j2 = pairs[p];
               
               if ((size1 == 0) && (size2 == 0))
               {
                  /* stacking pair of base pairs */
                  StackLoop st = { .i = i1, .j = j1};
                  ARRAY_PUSH(this->stack, st, StackLoop, { error = 1; });
               }
               else if ((size1 == 0) || (size2 == 0))
               {
                  /* bulge loop */
                  BulgeLoop bg = { .i1 = i1, .j1 = j1, .i2 = i2, .j2 = j2,
                                   .size = (size1 > size2 ? size1 : size2) };
                  ARRAY_PUSH(this->bulge_loop, bg, BulgeLoop, { error = 1; });
                  /*mfprintf (stderr, "Found bulge: Start %lu, %lu, stop: "
                    "%lu, %lu, size: %lu\n", i1, j1, i2, j2, bg.size);*/
               }
               else
               {
                  /* generic internal loop */
                  IntLoop in = { .i1 = i1, .j1 = j1, .i2 = i2, .j2 = j2,
                                 .size1 = size1, .size2 = size2};
                  ARRAY_PUSH(this->internal_loop, in, IntLoop, { error = 1; });
                  /*mfprintf (stderr, "Found internal: Start %lu, %lu, stop: "
                            "%lu, %lu, size: %lu, %lu\n", i1, j1, i2, j2,
                            size1, size2);*/
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
                  multiloop_find (i + 1, pairs[i] - 1, pairs, size, &ml);

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
                  /* bienert: ml.dangle5[ml.ndangle5][P5_Dangle] = i;
                     ml.dangle5[ml.ndangle5][P3_Dangle] = pairs[i];
                     ml.dangle5[ml.ndangle5][Ne_Dangle] = i + 1;*/
                  ml.ndangle5++;
                  ml.dangle3[ml.ndangle3][P5_Dangle] = pairs[i];
                  ml.dangle3[ml.ndangle3][P3_Dangle] = i;
                  ml.dangle3[ml.ndangle3][Ne_Dangle] = i + 1;
                     /* bienertml.dangle3[ml.ndangle3][P5_Dangle] = i;
                        ml.dangle3[ml.ndangle3][P3_Dangle] = pairs[i];
                        ml.dangle3[ml.ndangle3][Ne_Dangle] = pairs[i] - 1;*/
                  ml.ndangle3++;
                  
                  ARRAY_PUSH (this->multi_loop, ml, MultiLoop, { error = 1; });
               }
            }
         }
         
         /* move to next paired base */
         i = p;
      }
   }
   
   return error;
}



/*********************************   Access   *********************************/

/** @brief get the no. of hairpin loops of a 2D structure
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

/** @brief get the start base of the ith internal loop of a 2D structure
 *
 * @param[in] i Index of loop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_start_internal (const unsigned long i, const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->internal_loop));
   assert (ARRAY_CURRENT (this->internal_loop) > i);

   return ARRAY_ACCESS(this->internal_loop, i).i1;
}

/** @brief get the last base of the ith internal loop of a 2D structure
 *
 * @param[in] i Index of loop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_end_internal (const unsigned long i, const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->internal_loop));
   assert (ARRAY_CURRENT (this->internal_loop) > i);

   return ARRAY_ACCESS(this->internal_loop, i).j1;
}

/** @brief get the no. of unpaired bases of the ith internal loop
 *
 * @param[in] i Index of loop.
 * @param[in] this Secondary structure.
 */
unsigned long
secstruct_get_i_size_internal (const unsigned long i, const SecStruct* this)
{
   assert (this);
   assert (ARRAY_NOT_NULL (this->internal_loop));
   assert (ARRAY_CURRENT (this->internal_loop) > i);

   return ARRAY_ACCESS(this->internal_loop, i).size1;
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
secstruct_get_i_noof_stems_extloop (const SecStruct* this)
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
   assert (this->ext_loop.ndangle5 > i);

   return (this->ext_loop.dangle5[i][P3_Dangle]);
}
