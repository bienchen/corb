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
   A.size = (SIZE)

/** @brief Delete an array.
 *
 * Frees the memory allocated by the data component of an array and resets the
 * size counter. The current index is left undefined.
 *
 * @param[in] A array.
 * @param[in] SIZE size of the array.
 */
#define ARRAY_DELETE(A)   \
   if (ARRAY_NOT_NULL(A)) \
   {                      \
      XFREE(A.data);      \
      A.size = 0;         \
   }


/********************************   Altering   ********************************/

/** @brief Set an array to NULL values 
 *
 * Set all components of an array to null bytes. Intended to be used for
 * multiple arrays integrated in other data structures. Should make the code of
 * the constructor simple since ARRAY_DELETE can be called on "nulled" arrays,
 * safely.
 *
 * @param[in] A array.
 */
#define ARRAY_SET_VOID(A) \
   A.data = NULL; \
   A.current = 0; \
   A.size = 0

/** @brief Decrement the no. of elements in an array.
 *
 * Basically decrement the current index of an array by 1.
 *
 * @param[in] A array.
 */
#define ARRAY_CURRENT_DECREMENT(A) \
   assert (A.current > 0);         \
   A.current--

/** @brief Add an element to the end of an array.
 *
 * Stores an element in the data component of an array at the current index
 * position. If the index exceeds the size it will be adjusted. For a failed
 * reallocation, there is one argument reserved for a block of own code.
 *
 * @param[in] A array.
 * @param[in] D element to add.
 * @param[in] TYPE datatype of the array elements.
 * @param[in] FAILURE block of code to be executed if reallocation fails.
 */
#define ARRAY_PUSH(A, D, TYPE, FAILURE)                  \
   assert (A.data != NULL);                              \
   if (A.current >= A.size)                              \
   {                                                     \
      A.size = (A.size * 2) + 1;                         \
      A.data = XREALLOC(A.data, sizeof (TYPE) * A.size); \
      if (ARRAY_IS_NULL (A))                             \
         FAILURE                                         \
   }                                                     \
   A.data[A.current] = (D);                              \
   A.current++

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
#define ARRAY_SET(A, B, IDX)   \
   assert (A.current > (IDX)); \
   assert (A.data != NULL);    \
   A.data[(IDX)] = (B)

/** @brief Reset an array.
 *
 * Sets the current index of an array to 0.
 *
 * @param[in] A array.
 */
#define ARRAY_RESET(A) \
   assert (A.data != NULL); \
   A.current = 0


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

/** @brief Access a certain element in an array.
 *
 * Provides access to a certain element in an array. Does not check whether the
 * index in in or out of scope or if the data component exists anyway. 
 *
 * @param[in] B container.
 * @param[in] A array.
 * @param[in] IDX position.
 */
#define ARRAY_ACCESS(A, IDX) \
   A.data[(IDX)]


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


/*****************************  Standard Arrays   *****************************/
/**
 * This is a set of globally available arrays for the C standard data types.
 * All macros needing TYPE are mirrored here with type set. The macro name
 * changes to ARRAY_type_... where 'type' is written in uppercase letters.
 * Add missing standard types here.
 * Added so far:
 *  ArrayInt   - int
 *  ArrayUlong - unsigned long
 */

/** @brief ArrayUlong.
 *
 * Array for unsigned long.
 *
 */
ARRAY_CREATE_CLASS_NAMED(int, Int);

/** @brief Initialise an array for int.
 *
 * Allocation of memory for the data component of an array and initial setting
 * of its size and current index. It is up to you to check whether the data
 * component was successfully allocated or not afterwards.
 *
 * @param[in] A array.
 * @param[in] SIZE size of the array.
 */
#define ARRAY_INT_INIT(A, SIZE) \
   ARRAY_INIT(A, SIZE, int)

/** @brief Add an element to the end of an array.
 *
 * Stores an element in the data component of an array at the current index
 * position. If the index exceeds the size it will be adjusted. For a failed
 * reallocation, there is one argument reserved for a block of own code.
 *
 * @param[in] A array.
 * @param[in] D element to add.
 * @param[in] FAILURE block of code to be executed if reallocation fails.
 */
#define ARRAY_INT_PUSH(A, D, FAILURE) \
   ARRAY_PUSH(A, D, int, FAILURE)


/** @brief Array for unsigned long.
 */
ARRAY_CREATE_CLASS_NAMED(unsigned long, Ulong);

/** @brief Initialise an array for unsigned long.
 *
 * Allocation of memory for the data component of an array and initial setting
 * of its size and current index. It is up to you to check whether the data
 * component was successfully allocated or not afterwards.
 *
 * @param[in] A array.
 * @param[in] SIZE size of the array.
 */
#define ARRAY_ULONG_INIT(A, SIZE) \
   ARRAY_INIT(A, SIZE, unsigned long)

/** @brief Add an element to the end of an array.
 *
 * Stores an element in the data component of an array at the current index
 * position. If the index exceeds the size it will be adjusted. For a failed
 * reallocation, there is one argument reserved for a block of own code.
 *
 * @param[in] A array.
 * @param[in] D element to add.
 * @param[in] FAILURE block of code to be executed if reallocation fails.
 */
#define ARRAY_ULONG_PUSH(A, D, FAILURE) \
   ARRAY_PUSH(A, D, unsigned long, FAILURE)


#endif /* GENARRAY_H */

#ifdef __cplusplus
}
#endif
