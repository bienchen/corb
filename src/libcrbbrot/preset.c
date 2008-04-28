/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbrot/preset.c
 *
 *  @brief Storing fixed sites for a sequence matrix
 *
 *  Module: preset
 *
 *  Library: libcrbbrot
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-03-04
 *
 *
 *  Revision History:
 *         - 2008Mar04 bienert: created
 *
 */


#include <assert.h>
#include <libcrbbasic/crbbasic.h>
#include "preset.h"


struct Preset {
      char base;
      unsigned long pos;
};

struct PresetArray {
      Preset* set;
      unsigned long len;
      unsigned long size;
};


/**********************   Constructors and destructors   **********************/

/** @brief Function to wrap creators around.
 *
 * This function implements the constructor of a preset object together with
 * initialisation.\n
 * Returns a @c Preset object or @c NULL on problems.
 *
 * @param[in] base The base to set as fixed.
 * @param[in] pos Position of the fixed site.
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
static __inline__ Preset*
intl_preset_new (const char base,
                 const unsigned long pos,
                 const char* file, const int line)
{
   /* allocate 1 object */
   Preset* obj = XOBJ_MALLOC(sizeof (Preset), file, line);

   if (obj != NULL)
   {
      obj->base = base;
      obj->pos  = pos;
   }

   return obj;
}

/** @brief Create a new preset object.
 *
 * The constructor for @c Preset objects. If compiled with memory checking
 * enabled, @c file and @c line should point to the position where the function
 * was called. Both parameters are automatically set by using the macro
 * @c PRESET_NEW.\n
 * Returns @c NULL on error.
 *
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
Preset*
preset_new (const char* file, const int line)
{
   return intl_preset_new (CHAR_UNDEF, ULONG_UNDEF, file, line);
}

/** @brief Create a new preset object and initialise it.
 *
 * The constructor for @c Preset objects. If compiled with memory checking
 * enabled, @c file and @c line should point to the position where the function
 * was called. Both parameters are automatically set by using the macro
 * @c PRESET_NEW_PRESET.\n
 * Returns @c NULL on error.
 *
 * @param[in] base The base to set as fixed.
 * @param[in] pos Position of the fixed site.
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
Preset*
preset_new_preset (const char base,
                   const unsigned long pos,
                   const char* file, const int line)
{
   return intl_preset_new (base, pos, file, line);
}

/** @brief Function to wrap creators around.
 *
 * The internal constructor for @c PresetArray objects.\n
 * Returns @c NULL on error, a @c PresetArray object else.
 *
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
static PresetArray*
intl_presetarray_new (const unsigned long size,
                      const char* file, const int line)
{
   /* allocate 1 object */
   PresetArray* obj = XOBJ_MALLOC((sizeof (PresetArray)), file, line);

   if (obj != NULL)
   {
      obj->size = size;
      obj->len  = 0;
      obj->set  = XMALLOC (sizeof (Preset) * obj->size);
   }

   return obj;
}

/** @brief Delete a @c Preset object.
 *
 * The destructor for @c Preset objects.
 *
 * @param[in] obj object to be freed.
 */
void
preset_delete (Preset* obj)
{
   XFREE (obj);
}


/** @brief Create a new preset array object.
 *
 * The constructor for @c PresetArray objects. If compiled with memory checking
 * enabled, @c file and @c line should point to the position where the function
 * was called. Both parameters are automatically set by using the macro
 * @c PRESETARRAY_NEW.\n
 * Returns @c NULL on error.
 *
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
PresetArray*
presetarray_new (const char* file, const int line)
{
   return intl_presetarray_new (0, file, line);
}

/** @brief Create a new preset array object with a certain size.
 *
 * The constructor for @c PresetArray objects. If compiled with memory checking
 * enabled, @c file and @c line should point to the position where the function
 * was called. Both parameters are automatically set by using the macro
 * @c PRESETARRAY_NEW_SIZE.\n
 * Returns @c NULL on error.
 *
 * @param[in] size No. of elements in the array.
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
PresetArray*
presetarray_new_size (const unsigned long size,
                      const char* file, const int line)
{
   return intl_presetarray_new (size, file, line);
}

/** @brief Delete a @c PresetArray object.
 *
 * The destructor for @c PresetArray objects.
 *
 * @param[in] obj object to be freed.
 */
void
presetarray_delete (PresetArray* obj)
{
   if (obj != NULL)
   {
      XFREE (obj->set);
      XFREE (obj);
   }
}


/********************************   Altering   ********************************/

/** @brief Add a fixed site to an array.
 *
 * Add a new pair of base and position to an @c PresetArray.\n
 * Returns 0 on success, @c ERR_PA_ALLOC on error.
 *
 * @params[in] base Set type of fixed position.
 * @params[in] pos Fixed position.
 * @params[in] obj @c Preset array.
 */
int
presetarray_add (const char base, const unsigned long pos, PresetArray* obj)
{
   assert (obj);

   /* check whether obj is large enough */
   if (obj->size == obj->len)
   {
      obj->size = obj->size * 2;
      obj->set = XREALLOC (obj->set, (sizeof (Preset) * obj->size));
      if (obj->set == NULL)
      {
         return ERR_PA_ALLOC;
      }
   }

   obj->set[obj->len].base = base;
   obj->set[obj->len].pos = pos;
   obj->len++;

   return 0;
}


/*********************************   Access   *********************************/

/** @brief Get the fixed base of a certain element in a @c PresetArray object.
 *
 * Returns the base of the ith element of @c obj.
 *
 * @param[in] i Position in the array.
 * @param[in] obj The array.
 */
char
presetarray_get_ith_base (const unsigned long i, const PresetArray* obj)
{
   assert (obj);
   assert (obj->len > i);

   return obj->set[i].base;
}

/** @brief Get the fixed position of a element in a @c PresetArray object.
 *
 * Returns the fixed position in the ith element of @c obj.
 *
 * @param[in] i Position in the array.
 * @param[in] obj The array.
 */
unsigned long
presetarray_get_ith_pos (const unsigned long i, const PresetArray* obj)
{
   assert (obj);
   assert (obj->len > i);

   return obj->set[i].pos;
}


/**********************************   Size   **********************************/

/** @brief Get the number of elements in a @c PresetArray object.
 *
 * Returns the 'length', the number of elements in an array.
 *
 * @param[in] obj Array.
 */
unsigned long
presetarray_get_length (const PresetArray* obj)
{
   assert (obj);

   return obj->len;
}



/********************************   Searching   *******************************/


/********************************   Comparison   ******************************/
