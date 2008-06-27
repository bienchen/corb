/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbasic/genarray.h
 *
 *  @brief Macros to create internal arrays
 *
 *  Module:  none
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-06-26
 *
 *
 *  Revision History:
 *         - 2008Jun26 bienert: created
 *
 */


#ifdef __cplusplus
extern "C" {
#endif

#include <config.h>
#include <assert.h>

#ifndef GENARRAY_H
#define GENARRAY_H


/**********************       Define an array type       **********************/

/** @brief Define an array type.
 *
 * This macro creates a new array type of name 'ArraySUFFIX' able to store
 * values of type 'TYPE'. The array itself is stored in a component 'data', the
 * total number of possible elements as 'size' and the index of the current 
 * element stored in the array is identified via 'current'. 'size' and
 * 'current' are stored as type unsigned long.
 *
 * @param[in] TYPE datatype to be stored in the array.
 * @param[in] SUFFIX name of the array ('ArraySUFFIX').
 */
#define ARRAY_CREATE_CLASS_NAMED(TYPE, SUFFIX) \
   typedef struct {                            \
      TYPE *data;                              \
      unsigned long size;                      \
      unsigned long current;                   \
   } Array##SUFFIX

/** @brief Define an array type.
 *
 * This macro creates a new array type of name 'ArrayTYPE' able to store
 * values of type 'TYPE'. The array itself is stored in a component 'data', the
 * total number of possible elements as 'size' and the index of the last 
 * element stored in the array is identified via 'last'. 'size' and 'last' are
 * stored as type unsigned long.
 *
 * @param[in] TYPE datatype to be stored in the array and name.
 */
#define ARRAY_CREATE_CLASS(TYPE) \
    ARRAY_CREATE_CLASS_NAMED(TYPE,TYPE)


/**********************   Constructors and destructors   **********************/

/** @brief Initialise an array.
 *
 * Allocation of memory for the data component of an array and initial setting
 * of its size and current index. It is up to you to check whether the data
 * component was successfully allocated or not afterwards.
 *
 * @param[in] A array.
 * @param[in] SIZE size of the array.
 * @param[in] TYPE datatype of the array elements.
 */
#define ARRAY_INIT(A, SIZE, TYPE)            \
   A.data = XMALLOC(sizeof (TYPE) * (SIZE)); \
   A.current = 0;                            \
   A.size = SIZE

/** @brief Delete an array.
 *
 * Frees the memory allocated by the data component of an array and resets the
 * size counter. The current index is left undefined.
 *
 * @param[in] A array.
 * @param[in] SIZE size of the array.
 * @param[in] TYPE datatype of the array elements.
 */
#define ARRAY_DELETE(A)   \
   if (ARRAY_NOT_NULL(A)) \
   {                      \
      XFREE(A.data);      \
      A.size = 0;         \
   }


/********************************   Altering   ********************************/

/** @brief Add an element to the end of an array.
 *
 * Stores an element in the data component of an array at the current index
 * position. If the index exceeds the size it will be adjusted. It is up to you
 * to check whether the data component was successfully reallocated or not
 * afterwards.
 *
 * @param[in] A array.
 * @param[in] D element to add.
 * @param[in] TYPE datatype of the array elements.
 */
#define ARRAY_PUSH(A, D, TYPE)                                         \
   assert (A.data != NULL);                                                    \
   if (A.current < A.size)                                             \
   {                                                                   \
      A.data[A.current] = D;                                           \
      A.current++;                                                     \
   }                                                                   \
   else                                                                \
   {                                                                   \
      A.data = XREALLOC(A.data, sizeof (TYPE) * ((A.size * 2) + 1));   \
      if (ARRAY_NOT_NULL(A))                                           \
      {                                                                \
         A.size = (A.size * 2) + 1;                                    \
         A.data[A.current] = D;                                        \
         A.current++;                                                  \
      }                                                                \
   }

/** @brief Pops the last element of an array.
 *
 * Stores the last element of a stack into a variable and deletes the entry.
 * You have to assure that the array has at least one element stored for this
 * operation.
 *
 * @param[in] B container to store the element.
 * @param[in] A array.
 */
#define ARRAY_POP(B, A)     \
   assert (A.current > 0);  \
   assert (A.data != NULL); \
   A.current--;             \
   B = A.data[A.current]

/** @brief Store an element at a certain position in an array.
 *
 * Stores the given element at defined position in array. You have to assure,
 * that the index is valid and the data component was allocated before.
 *
 * @param[in] A array.
 * @param[in] B element.
 * @param[in] IDX position.
 */
#define ARRAY_SET(A, B, IDX)  \
   assert (A.current > IDX); \
   assert (A.data != NULL);     \
   A.data[IDX] = B


/*********************************   Access   *********************************/

/** @brief Store the element at a certain position in an array in a variable.
 *
 * Stores the element at defined position in array in a variable. You have to
 * assure, that the index is valid and the data component was allocated before.
 *
 * @param[in] B container.
 * @param[in] A array.
 * @param[in] IDX position.
 */
#define ARRAY_GET(B, A, IDX)   \
   assert (A.current > (IDX)); \
   assert (A.data != NULL);    \
   B = A.data[(IDX)]

/** @brief Get the current index of an array.
 *
 * Get the current index of an array.
 *
 * @param[in] A array.
 */
#define ARRAY_CURRENT(A) \
   A.current

/** @brief Get the size of an array.
 *
 * Get the size of an array.
 *
 * @param[in] A array.
 */
#define ARRAY_SIZE(A) \
   A.size


/*********************************    Misc    *********************************/

/** @brief Check if the data component of an array is not NULL.
 *
 * Check if the data component of an array is not NULL.
 *
 * @param[in] A array.
 */
#define ARRAY_NOT_NULL(A) \
   A.data != NULL

/** @brief Check if the data component of an array is NULL.
 *
 * Check if the data component of an array is NULL.
 *
 * @param[in] A array.
 */
#define ARRAY_IS_NULL(A) \
   A.data == NULL

#endif /* GENARRAY_H */

#ifdef __cplusplus
}
#endif
