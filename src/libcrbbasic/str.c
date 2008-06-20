/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbasic/str.c
 *
 *  @brief String class.
 *
 *  Module: str
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-02-05
 *
 *
 *  Revision History:
 *         - 2008Feb05 bienert: created
 *
 */

 
#include <config.h>
#include <stddef.h>
#include <limits.h>
#include <assert.h>
#include "gcckywrd.h"
#include "inc_strg.h"
#include "inc_bool.h"
#include "errormsg.h"
#include "memmgr.h"
#include "mprintf.h"
#include "str.h"


struct Str {
      unsigned long len;
      size_t size;
      char* data;
};


/**********************   Constructors and destructors   **********************/

/** @brief Create a new string object.
 *
 * The constructor for @c Str objects. If compiled with enabled
 * memory checking, @c file and @c line should point to the position where the
 * function was called. Both parameters are automatically set by using the
 * macro @c STR_NEW.\n
 * Returns @c NULL on error.
 *
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
Str*
str_new (const char* file, const int line)
{
   /* allocate 1 object */
   Str* obj = XOBJ_MALLOC(sizeof (Str), file, line);
   
   if (obj != NULL)
   {
      obj->len  = 0;
      obj->size = 0;
      obj->data = NULL;
   }

   return obj;
}

/** @brief Function to wrap creators around.
 *
 * This function implements the constructor of a string together with
 * initialisation. Setting a c string or a string on creation is nearby the
 * same thing. This function deals with the only real difference: How to get
 * the size of the new string.\n
 * Returns a string object or @c NULL on problems.
 *
 * @param[in] cstr The initial string.
 * @param[in] size Size of the initial string (length + 1 for the term. '\0').
 * @param[in] file Calling file.
 * @param[in] line Calling line.
 */
static __inline__ Str*
istr_new_init (const char* cstr, size_t size, const char* file, const int line)
{
   Str* this = str_new (file, line);

   if (this != NULL)
   {
      this->size = size;
      this->data = XMALLOC(sizeof (char) * this->size);
      if (this->data != NULL)
      {
         this->len = this->size - 1;
         strncpy (this->data, cstr, this->len);
         this->data[this->len] = '\0';
      }
      else
      {
         XFREE (this);
         this = NULL;
      }
   }

   return this;
}


/** @brief Create a new string object from a char array.
 *
 * A constructor an initialiser for @c Str objects. If compiled with enabled
 * memory checking, @c file and @c line should point to the position where the
 * function was called. Both parameters are automatically set by using the
 * macro @c STR_NEW_CSTR. @c cstr is stored as a copy of the input parameter.\n
 * Returns @c NULL if failed to allocate memory for the object.
 *
 * @param[in] cstr c string used for initialisation (must not be @c NULL).
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
Str*
str_new_cstr (const char* cstr, const char* file, const int line)
{
   return istr_new_init (cstr, strlen (cstr) + 1, file, line);
}


/** @brief Create a new string object of certain size filled with a single char.
 *
 * A constructor an initialiser for @c Str objects. If compiled with enabled
 * memory checking, @c file and @c line should point to the position where the
 * function was called. Both parameters are automatically set by using the
 * macro @c STR_NEW_CHAR. @c cstr is stored as a copy of the input parameter.\n
 * Returns @c NULL if failed to allocate memory for the object.
 *
 * @param[in] c char used to fill the string.
 * @param[in] l desired length of the string.
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
Str*
str_new_char (const char c,
              const unsigned long l,
              const char* file,
              const int line)
{
   Str* this = str_new (file, line);

   if (this != NULL)
   {
      this->size = l + 1;
      this->data = XMALLOC(this->size * sizeof (char));
      if (this->data != NULL)
      {
         this->len = this->size - 1;
         memset (this->data, c, sizeof (*(this->data)) * this->len);
         this->data[this->len] = '\0';
      }
      else
      {
         XFREE (this);
         this = NULL;
      }
   }

   return this;
}


/** @brief Create a new string object as copy of a string.
 *
 * A constructor an initialiser for @c Str objects. If compiled with enabled
 * memory checking, @c file and @c line should point to the position where the
 * function was called. Both parameters are automatically set by using the
 * macro @c STR_NEW_STR. @c str is stored as a copy of the input parameter.\n
 * Returns @c NULL on error.
 *
 * @param[in] str string used for initialisation (must not be @c NULL).
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
Str*
str_new_str (const Str* str, const char* file, const int line)
{
   return istr_new_init (str->data, str->size, file, line);
}


/** @brief Delete a string.
 *
 * The destructor for @c Str objects.
 *
 * @param[in] this object to be freed.
 */
void
str_delete (Str* this)
{

   if (this != NULL)
   {
     /* this->size = 0; */
     /* this->len = 0; */
     XFREE(this->data);
     XFREE(this);
   }
}


/********************************   Altering   ********************************/

/** @brief Copy a character string to a string.
 * 
 * This function is intended as base for setting/ copying strings. These
 * functions could be wrapped around this one.\n
 * Returns 0 on success, @c ERR_STR_ALLOC on problems with allocating memory
 * for the character string.
 *
 * @param[in/out] dest String object to be altered.
 * @param[in] cstr C string to be copied.
 * @param[in] len Length of the c string.
 */
static __inline__ int
istr_cpy_cstr (Str* dest, const char* cstr, const size_t len)
{
   /* check need for reallocation */
   if (dest->size <= len)
   {
      dest->size = len + 1;
      dest->data = XREALLOC (dest->data, (sizeof (char) * dest->size));

      if (dest->data == NULL)
      {
         return ERR_STR_ALLOC;
      }      
   }

   /* copy */
   strncpy (dest->data, cstr, len);
   dest->len = len;
   dest->data[dest->len] = '\0';

   return 0;   
}


/** @brief Copy a string.
 *
 * Copy one string into another. @c dest could be either an empty (new) string
 * or an already used one. In the latter case the content will be overwritten.\n
 * Returns 0 on success, @c ERR_STR_ALLOC on problems with allocating memory
 * for the character string.
 *
 * @param[out] dest Copy @c src to (must not be @c NULL).
 * @param[in] src   Copy from (must not be @c NULL).
 */
int
str_cpy (Str* dest, const Str* src)
{
   return istr_cpy_cstr (dest, src->data, src->len);
}


/** @brief Assign a char array to a string.
 *
 * Copy a c string into a string. Overwrites the old content of the string.\n
 * Returns 0 on success, @c ERR_STR_ALLOC on memory allocation problems.\n
 *
 * @param[out] dest String to copy @c src to (must not be @c NULL).
 * @param[in]  src C string to be copied (must not be @c NULL).
 */
int
str_set (Str* dest, const char* src)
{
   return istr_cpy_cstr (dest, src, strlen (src));
}


/** @brief Set a character at a certain position.
 *
 * Sets character @c c at position @c i in string @c dest. Does not perform any
 * boundary checking. You should be really careful using this function because
 * you might overwrite the terminating '\0' of a string. While this is no
 * problem for copying/ appending strings, it could be for outputting them.\n
 *
 * @param[out] dest String to be altered (must not be @c NULL).
 * @param[in]  i Position to be altered.
 * @param[in]  c Character to be set.
 */
void
str_set_i (Str* dest, const size_t i, const char c)
{
   dest->data[i] = c;
}


/** @brief Set a certain character at a certain position in a string.
 *
 * Sets character @c c at position @c i in string @c dest.\n
 * Returns 0 on success, @c ERR_STR_OOR (out of range), if
 * @c i >= @c str_length (@c dest).
 *
 * @param[out] dest String to be altered (must not be @c NULL).
 * @param[in]  i    Position to be altered.
 * @param[in]  c    Character to be set.
 */
int
str_at (Str* dest, const unsigned long i, const char c)
{
   if (i >= dest->len)
   {
      THROW_ERROR_MSG ("Attempt to write beyond string capacity (%lu): %lu",
                       dest->len, i);
      return ERR_STR_OOR;
   }

   dest->data[i] = c;

   return 0;
}


/** @brief General appending function.
 *
 * This function appends a c string to a string. It is intended to be used by
 * the public appending functions.\n
 * Returns 0 on success, @c ERR_STR_ALLOC on problems with allocating memory
 * for the character string.
 *
 * @param[in/out] dest String to append to.
 * @param[in] cstr C string to be appended.
 * @param[in] len Length of @c cstr.
 */
static __inline__ int
istr_append (Str* dest, const char* cstr, const unsigned long len)
{
   size_t con_size = dest->len + len + 1;

   /* check need of reallocation */
   if (dest->size < con_size)
   {
      dest->size = con_size;
      dest->data = XREALLOC (dest->data, (sizeof (char) * dest->size));
      
      if (dest->data == NULL)
      {
         return ERR_STR_ALLOC;
      }
   }

   strncat (dest->data, cstr, len);
   dest->len = con_size - 1;

   return 0;
}


/** @brief Append a string to another.
 *
 * Concatenates @c str2 on the end of @c str1.\n
 * Returns 0 on success, @c ERR_STR_ALLOC on memory problems.
 *
 * @param[in/out] str1 First string (must not be @c NULL).
 * @param[in] str2 Second string (must not be @c NULL).
 */
int
str_append_str (Str* str1, const Str* str2)
{
   return istr_append (str1, str2->data, str2->len);
}


/** @brief Append a character string to a string.
 *
 * Concatenates @c cstr to the end of @c str.\n
 * Returns 0 on success, @c ERR_STR_ALLOC on memory problems.
 *
 * @param[in/out] str String object.
 * @param[in] cstr Character string to be appended.
 */
int
str_append_cstr (Str* str, const char* cstr)
{
   return istr_append (str, cstr, strlen (cstr));
}


/** @brief Assign a substring of a character string to a string object.
 *
 * Assigns the substring of length @c len starting at position @c start of the
 * character string @c cstr to @c str.\n
 * Returns @c ERR_STR_ALLOC on memory problems.
 *
 * @param[in/out] str String object.
 * @param[in] cstr Source for the substring.
 * @param[in] start Start position in @c cstr.
 * @param[in] len Length of the substring.
 */
int
str_assign_csubstr (Str* str,
                    const char* cstr,
                    const size_t start,
                    const size_t len)
{
   return istr_cpy_cstr (str, (cstr + start), len);
}


/** @brief Assign a substring of a character string to a string object.
 *
 * Assigns the substring of length @c len starting at position @c start of the
 * character string @c cstr to @c str.\n
 * Returns @c ERR_STR_ALLOC on memory problems.
 *
 * @param[in/out] str1 String object.
 * @param[in] str2 Source for the substring.
 * @param[in] start Start position in @c str2.
 * @param[in] len Length of the substring.
 */
int
str_assign_substr (Str* str1,
                   const Str* str2,
                   const size_t start,
                   const size_t len)
{
   return istr_cpy_cstr (str1, (str2->data + start), len);
}


/** @brief Remove all characters from string.
 *
 * @param[in/out] str String object.
 */
void
str_clear (Str* str)
{
   if (str->len > 0)
   {
      str->data[0] = '\0';
   }
   str->len = 0;
}


/*********************************   Access   *********************************/

/** @brief Get a character string stored in a string object.
 *
 * Returns the character string stored in a string object.
 *
 * @param[in] str String object.
 */
const char*
str_get (const Str* str)
{
   return str->data;
}


/** @brief Return the character at position i.
 *
 * The character is returned without any boundary checking.
 *
 * @param[in] str String object.
 * @param[in] i Position.
 */
char
str_get_i (const Str* str, const unsigned long i)
{
   return str->data[i];
}


/**********************************   Size   **********************************/

/** @brief Returns the length of a string.
 *
 * @param[in] str String object.
 */
unsigned long
str_length (const Str* str)
{
   return str->len;
}


/** @brief Returns the size of a string.
 *
 * The size of a string is not the same as the length. It describes the number
 * of characters which can be stored in a string (without reallocation).
 *
 * @param[in] str String object.
 */
size_t
str_capacity (const Str* str)
{
   return str->size;
}


/** @brief Check whether a string is empty.
 *
 * Returns @c true, if @c str is empty, @c false else. This is NOT a test if
 * @c str is @c NULL, but for a string to hold a character chain or not.
 *
 * @param[in] str String object.
 */
bool
str_empty (const Str* str)
{
   if (str->len == 0)
   {
      return true;
   }

   return false;
}


/** @brief Force a resizing of a string.
 *
 * Changes size (capacity) of @c str to @c newsize, padding with @c padchar if
 * necessary.\n
 * Returns 0 on success, ERR_STR_ALLOC on memory problems.
 *
 * @param[in/out] str String object.
 * @param[in] newsize New size of the string object.
 * @param[in] padchar A character to append for growing strings.
 */
int
str_resize (Str* str, const size_t newsize, const char padchar)
{
   unsigned long i;

   /* resize */
   str->size = newsize;
   str->data = XREALLOC (str->data, (sizeof (char) * str->size));
      
   if (str->data == NULL)
   {
      return ERR_STR_ALLOC;
   }

   /* padding if newsize > length */
   if (str->size > 0)
   {
      for (i = str->size - 1; i >= str->len; i--)
      {
         str->data[i] = padchar;
      }
   }
   str->len = str->size - 1;
   str->data[str->len] = '\0';

   return 0;
}


/********************************   Searching   *******************************/

/** @brief Find the leftmost position of a character in a string.
 *
 * Returns the position of leftmost occurrence of @c c in @c str. For positions
 * the counting starts at 1. Therefore a return value of 0 means that @c c does
 * not occur in @c str.
 *
 * @param[in] str String object.
 * @param[in] c Character to search.
 */
unsigned long
str_find_c (const Str* str, const char c)
{
   unsigned long i;

   for (i = 0; i < str->len; i++)
   {
      if (str->data[i] == c)
      {
         return i;
      }
   }

   return 0;
}


/** @brief Find the rightmost position of a character in a string.
 *
 * Returns the position of rightmost occurrence of @c c in @c str. For positions
 * the counting starts at 1. Therefore a return value of 0 means that @c c does
 * not occur in @c str.
 *
 * @param[in] str String object.
 * @param[in] c Character to search.
 */
unsigned long
str_rfind_c (const Str* str, const char c)
{
   unsigned long i;

   if (str->len > 0)
   {
      for (i = (str->len - 1); i > 0; i--)
      {
         if (str->data[i] == c)
         {
            return i;
         }
      }
   }

   return 0;
}


/** @brief Find the leftmost start position of a character string in a string.
 *
 * Returns the start position of @c cstr in @c str. For positions the counting
 * starts at 1. Therefore a return value of 0 means that @c cstr does
 * not occur in @c str.
 *
 * @param[in] str String object.
 * @param[in] cstr Character string.
 * @param[in] len Length of the string.
 */
static __inline__ unsigned long
istr_find_cstr (const Str* str, const char* cstr, size_t len)
{
   unsigned long i;
   unsigned long j = 0;
   unsigned long jshift = 0;

   if (str->len >= len)
   {
      for (i = 0; i < str->len; i++)
      {
         if (cstr[j] == str->data[i])
         {
            if ((j > 0) && (cstr[jshift] == str->data[i]))
            {
               jshift++;
            }
            else
            {
               jshift = 0;
            }
            j++;            
         }
         else
         {
            j = jshift;
            i += jshift;
            jshift = 0;
         }

         if (j == len)
         {
            return i + 2 - len;
         }
      }
   }
   
   return 0;
}

/*while ((j < len) && (cstr[j] == str->data[i]))
  {
  if ((j > 0) && (cstr[jshift] == str->data[i]))
  {
  jshift++;
  }
  else
  {
  jshift = 0;
  }
  j++;
  i++;
  }*/

/** @brief Find the leftmost occurrence of a character string in a string.
 *
 * Returns the first start position of @c cstr in @c str. For positions the
 * counting starts at 1. Therefore a return value of 0 means that @c cstr does
 * not occur in @c str.
 *
 * @param[in] str String object.
 * @param[in] cstr Character string.
 */
unsigned long
str_find_cstr (const Str* str, const char* cstr)
{
   return istr_find_cstr (str, cstr, strlen (cstr));
}


/** @brief Find the leftmost occurrence of a string in a string.
 *
 * Returns the first start position of @c str2 in @c str1. For positions the
 * counting starts at 1. Therefore a return value of 0 means that @c cstr does
 * not occur in @c str.
 *
 * @param[in] str1 String object.
 * @param[in] str2 Search string.
 */
unsigned long
str_find_str (const Str* str1, const Str* str2)
{
   return istr_find_cstr (str1, str2->data, (size_t) str2->len);
}


/** @brief Find the rightmost start position of a character string in a string.
 *
 * Returns the start position of @c cstr in @c str. For positions the counting
 * starts at 1. Therefore a return value of 0 means that @c cstr does
 * not occur in @c str.
 *
 * @param[in] str String object.
 * @param[in] cstr Character string.
 * @param[in] len Length of the string.
 */
static __inline__ unsigned long
istr_rfind_cstr (const Str* str, const char* cstr, size_t len)
{
   unsigned long i;
   unsigned long j = len;
   unsigned long jshift = j;

   if ((len > 0) && (str->len >= len))
   {
      for (i = str->len; i > 0; i--)
      {
         if (cstr[j-1] == str->data[i-1])
         {
            if ((j < len) && (cstr[jshift - 1] == str->data[i-1]))
            {
               jshift--;
            }
            else
            {
               jshift = len;
            }
            j--;
         }
         else
         {
            j = jshift;
            i = i-(len-jshift);
            jshift = len;
         }

         if (j == 0)
         {
            return i;
         }
      }
   }

   return 0;
}

/*static __inline__ unsigned long
istr_rfind_cstr (const Str* str, const char* cstr, size_t len)
{
   unsigned long i, j = 1;
   unsigned long jshift = len, ishift;

   if ((len > 0) && (str->len >= len))
   {
      for (i = str->len; i > len - 1; i--)
      {
         j = jshift;
         ishift = 0;
         jshift = len;

         while ((j > 0) && (cstr[j-1] == str->data[i-1-(len-j)]))
         {
            if ((j < len) && (cstr[jshift - 1] == str->data[i-1-(len-j)]))
            {
               jshift--;
               ishift = len - j;
            }
            else
            {
               jshift = len;
               ishift = len - j;
            }

            j--;
         }
         if (j == 0)
         {
            return i + 1 - len;
         }
         i = i - ishift;
      }
   }
   
   return 0;
}*/


/** @brief Find the rightmost occurrence of a character string in a string.
 *
 * Returns the last start position of @c cstr in @c str. For positions the
 * counting starts at 1. Therefore a return value of 0 means that @c cstr does
 * not occur in @c str.
 *
 * @param[in] str String object.
 * @param[in] cstr Character string.
 */
unsigned long
str_rfind_cstr (const Str* str, const char* cstr)
{
   return istr_rfind_cstr (str, cstr, strlen (cstr));
}


/** @brief Find the rightmost occurrence of a string in a string.
 *
 * Returns the last start position of @c str2 in @c str1. For positions the
 * counting starts at 1. Therefore a return value of 0 means that @c cstr does
 * not occur in @c str.
 *
 * @param[in] str1 String object.
 * @param[in] str2 Search string.
 */
unsigned long
str_rfind_str (const Str* str1, const Str* str2)
{
   return istr_rfind_cstr (str1, str2->data, (size_t) str2->len);
}


/** @brief Preprocess a character set
 *
 * For each character in the alphabet (here defined by the size of type
 * @c char), the @c table holds a bit. If the character is available in the @c
 * set, the bit is set to 1. All other bits are set to 0.
 *
 * @param[in/out] table Table to process.
 * @param[in] size Size of the table.
 * @param[in] set The set of characters.
 * @param[in] setlen Size of the set.
 * @param[in] bits No. of bits in the base type of the table.
 */
static __inline__ void
init_alpha_table_true (char* table,
                       const size_t size,
                       const char* set,
                       const size_t setlen,
                       const size_t bits)
{
   unsigned long i;
   int b;

   memset(table, 0, size * sizeof(*table));

   for (i = 0; i < setlen; i++)
   {
      b = (int) set[i] / bits; /* calc bin */

      /* set character in bin */
      table[b] |= (char) 1<<(set[i] - (b * bits));
   }
}


/** @brief Preprocess a character set
 *
 * For each character in the alphabet (here defined by the size of type
 * @c char), the @c table holds a bit. If the character is available in the @c
 * set, the bit is set to 0. All other bits are set to 1.
 *
 * @param[in/out] table Table to process.
 * @param[in] size Size of the table.
 * @param[in] set The set of characters.
 * @param[in] setlen Size of the set.
 * @param[in] bits No. of bits in the base type of the table.
 */
static __inline__ void
init_alpha_table_false (char* table,
                        const size_t size,
                        const char* set,
                        const size_t setlen,
                        const size_t bits)
{
   unsigned long i;
   int b;

   memset(table, (~0), size * sizeof(*table));

   for (i = 0; i < setlen; i++)
   {
      b = (int) set[i] / bits; /* calc bin */

      /* set character in bin */
      table[b] &= (char) ~(1<<(set[i] - (b * bits)));
   }
}

/** @brief Evaluate a set of characters against a string starting at the front.
 *
 * Returns the first position of a character from @c set matching a certain
 * criteria with @c str. The search starts at @c start. Positions in the string
 * are counted from 1. Therefore a return value of 0 means that no character
 * from the set occurs in the string. As an exception, 0 is a valid start
 * position and is treated in the same way as @c start = 1. According to the
 * initialisation of the table with @c init_table, the first position where a
 * character from a set occurs or not can be found.
 *
 * @param[in] str The string.
 * @param[in] set The character set.
 * @param[in] setlen Size of the character set.
 * @param[in] start Start position for the search.
 * @param[in] init_table Initialisation function for the alphabet table.
 */
static unsigned long
istr_eval_first_of (const Str* str,
                    const char* set,
                    const size_t setlen,
                    const unsigned long start,
                    void (*init_table) (char*,
                                        const size_t,
                                        const char*,
                                        const size_t,
                                        const size_t))
{
   /* create hash table of alphabet's size */
   size_t bits = CHAR_BIT;
   size_t bins = (CHAR_MAX + 1) / bits;
   char alpha_dist[(CHAR_MAX + 1) / CHAR_BIT];
   unsigned long i;
   int b;

   init_table (alpha_dist, bins, set, setlen, bits);

   /* searching */
   if (start == 0)
   {
      i = 0;
   }
   else
   {
      i = start - 1;
   }
   for (; i < str->len; i++)
   {
      b = (int) str->data[i] / bits; /* calc bin */

      if ((char) 1<<(str->data[i]-(b*bits)) & alpha_dist[b])
      {
         return i + 1;
      }
   }

   return 0;
}


/** @brief Find first position of a character from a set in a string
 *
 * Searches the first position of a character from @c set in @c str. The search
 * starts at @c start. Positions in the string are counted from 1. Therefore a
 * return value of 0 means that no character from the set occurs in the string.
 * As an exception, 0 is a valid start position and is treated in the same way
 * as @c start = 1.\n
 * Returns the position of the first occurrence of a member of the set or 0.
 *
 * @param[in] str The string.
 * @param[in] set The character set.
 * @param[in] start Start position for the search.
 */
unsigned long
str_find_first_of_cstr (const Str* str,
                        const char* set,
                        const unsigned long start)
{
   return istr_eval_first_of (str,
                              set,
                              strlen (set),
                              start,
                              init_alpha_table_true);
}


/** @brief Find first position of a character from a set in a string
 *
 * Searches the first position of a character from @c set in @c str. The search
 * starts at @c start. Positions in the string are counted from 1. Therefore a
 * return value of 0 means that no character from the set occurs in the string.
 * As an exception, 0 is a valid start position and is treated in the same way
 * as @c start = 1.\n
 * Returns the position of the first occurrence of a member of the set or 0.
 *
 * @param[in] str The string.
 * @param[in] set The character set.
 * @param[in] start Start position for the search.
 */
unsigned long
str_find_first_of_str (const Str* str,
                       const Str* set,
                       const unsigned long start)
{
   return istr_eval_first_of (str,
                              set->data,
                              set->len,
                              start,
                              init_alpha_table_true);
}


/** @brief Find first position in a string where no character from a set matches
 *
 * Searches the first position in @c str which does not contain any character
 * out of @c set. The search starts at @c start. Positions in the string are
 * counted from 1. Therefore a return value of 0 means that no character from
 * the set occurs in the string. As an exception, 0 is a valid start position
 * and is treated in the same way as @c start = 1.\n
 * Returns the position of the first "non-occurrence" of a member of the set or
 * 0.
 *
 * @param[in] str The string.
 * @param[in] set The character set.
 * @param[in] start Start position for the search.
 */
unsigned long
str_find_first_not_of_cstr (const Str* str,
                            const char* set,
                            const unsigned long start)
{
   return istr_eval_first_of (str,
                              set,
                              strlen (set),
                              start,
                              init_alpha_table_false);
}


/** @brief Find first position in a string where no character from a set matches
 *
 * Searches the first position in @c str which does not contain any character
 * out of @c set. The search starts at @c start. Positions in the string are
 * counted from 1. Therefore a return value of 0 means that all characters from
 * @c str occur in @c set. As an exception, 0 is a valid start position
 * and is treated in the same way as @c start = 1.\n
 * Returns the position of the first "non-occurrence" of a member of the set or
 * 0.
 *
 * @param[in] str The string.
 * @param[in] set The character set.
 * @param[in] start Start position for the search.
 */
unsigned long
str_find_first_not_of_str (const Str* str,
                           const Str* set,
                           const unsigned long start)
{
   return istr_eval_first_of (str,
                              set->data,
                              set->len,
                              start,
                              init_alpha_table_false);
}


/** @brief Evaluate a set of characters against a string starting from the end.
 *
 * Returns the last position of a character from @c set matching a certain
 * criteria with @c str. The search runs up to position @c last. Positions in
 * the string are counted from 1. Therefore a return value of 0 means that all
 * characters from @c str occur in @c set. As an exception, 0 is valid for
 * @c last and is treated like @c last = 1.
 *
 * @param[in] str The string.
 * @param[in] set The character set.
 * @param[in] setlen Size of the character set.
 * @param[in] last Last position to be considered in the search.
 * @param[in] init_table Initialisation function for the alphabet table.
 */
static unsigned long
istr_eval_last_of (const Str* str,
                   const char* set,
                   const size_t setlen,
                   const unsigned long last,
                   void (*init_table) (char*,
                                       const size_t,
                                       const char*,
                                       const size_t,
                                       const size_t))
{
   /* create hash table of alphabet's size */
   size_t bits = CHAR_BIT;
   size_t bins = (CHAR_MAX + 1) / bits;
   char alpha_dist[(CHAR_MAX + 1) / CHAR_BIT];
   unsigned long i;
   int b;

   init_table (alpha_dist, bins, set, setlen, bits);

   /* searching */
   if (last <= str->len)
   {
      i = last;
      if (i == 0)
      {
         i++;
      }
   }
   else
   {
      return 0;
   }

   for (; i > 0; i--)
   {
      b = (int) str->data[i-1] / bits; /* calc bin */

      if ((char) 1<<(str->data[i-1]-(b*bits)) & alpha_dist[b])
      {
         return i;
      }
   }

   return 0;
}


/** @brief Find last position of a character from a set in a string
 *
 * Searches the last position of a character from @c set in @c str. The search
 * runs up to position @c last. Positions in the string are counted from 1.
 * Therefore a return value of 0 means that no character from the set occurs in
 * the string. As an exception, 0 is valid for @c last and is treated in the
 * same way as @c last = 1.\n
 * Returns the position of the first occurrence of a member of the set or 0.
 *
 * @param[in] str The string.
 * @param[in] set The character set.
 * @param[in] last Last position to be considered in the search.
 */
unsigned long
str_find_last_of_cstr (const Str* str,
                       const char* set,
                       const unsigned long last)
{
   return istr_eval_last_of (str,
                             set,
                             strlen (set),
                             last,
                             init_alpha_table_true);
}


/** @brief Find last position of a character from a set in a string
 *
 * Searches the last position of a character from @c set in @c str. The search
 * runs up to position @c last. Positions in the string are counted from 1.
 * Therefore a return value of 0 means that no character from the set occurs in
 * the string. As an exception, 0 is valid for @c last and is treated in
 * the same way as @c last = 1.\n
 * Returns the position of the first occurrence of a member of the set or 0.
 *
 * @param[in] str The string.
 * @param[in] set The character set.
 * @param[in] last Last position to be considered in the search.
 */
unsigned long
str_find_last_of_str (const Str* str,
                      const Str* set,
                      const unsigned long last)
{
   return istr_eval_last_of (str,
                             set->data,
                             set->len,
                             last,
                             init_alpha_table_true);
}


/** @brief Find last position in a string where no character from a set matches
 *
 * Searches the last position in @c str which does not contain any character
 * out of @c set. The search runs up to position @c last. Positions in the
 * string are counted from 1. Therefore a return value of 0 means that all
 * character from @c str are found in @c set. As an exception, 0 is valid for
 * @c last and is treated in the same way as @c last = 1.\n
 * Returns the position of the first "non-occurrence" of a member of the set or
 * 0.
 *
 * @param[in] str The string.
 * @param[in] set The character set.
 * @param[in] last Last position to be considered in the search.
 */
unsigned long
str_find_last_not_of_cstr (const Str* str,
                           const char* set,
                           const unsigned long last)
{
   return istr_eval_last_of (str,
                             set,
                             strlen (set),
                             last,
                             init_alpha_table_false);
}


/** @brief Find first position in a string where no character from a set matches
 *
 * Searches the first position of in @c str which does not contain any
 * character out of @c set. The search starts at @c start. Positions in the
 * string are counted from 1. Therefore a return value of 0 means that all
 * characters from @c str occur in @c set. As an exception, 0 is a valid for
 * @c last and is treated in the same way as @c last = 1.\n
 * Returns the position of the first "non-occurrence" of a member of the set or
 * 0.
 *
 * @param[in] str The string.
 * @param[in] set The character set.
 * @param[in] start Start position for the search.
 */
unsigned long
str_find_last_not_of_str (const Str* str,
                          const Str* set,
                          const unsigned long last)
{
   return istr_eval_last_of (str,
                             set->data,
                             set->len,
                             last,
                             init_alpha_table_false);
}


/********************************   Comparison   ******************************/

/** @brief Compare two substrings.
 *
 * Compares two substrings and returns an integer indicating the differences.
 * The substrings are fetched from a string object and a character string. 
 * @c str_start, @c cstr_start define the starting points and @c str_len,
 * @c cstr_len the length of the substrings. If one of the parameters is out of
 * bounds, it is automatically set to the maximal possible value. E.g. if
 * @c str_start is larger than the length of the string, it will be set to the
 * string length. Positions start at 1 but 0 is acceptable as an exception for
 * the first position.\n
 * Returns a value less than zero when @c str is less than @c cstr. A return
 * value of zero signals equal strings and a value greater than zero signals,
 * that @c str is greater than @c cstr.
 *
 * @param[in] str String object.
 * @param[in] str_start Start position in the string.
 * @param[in] str_len Length of substring to compare.
 * @param[in] cstr Character string.
 * @param[in] strlen_cstr Length of the character string
 * @param[in] cstr_start Start position in the character string.
 * @param[in] cstr_len Length of substring to compare.
 */
static __inline__ int
istr_compare (const Str* str,
              size_t str_start,
              size_t str_len,
              const char* cstr,
              const size_t strlen_cstr, 
              size_t cstr_start,
              size_t cstr_len)
{
   size_t i;
   size_t mlen;

   /* boundary checking for the string */
   if (str_start != 0)
   {
      str_start--;
   }

   if (str_start > str->len)
   {
      str_start = str->len;
   }

   if (str->len - str_start < str_len)
   {
      str_len = str->len - str_start;
   }

   /* boundary checking for the cstring */
   if (cstr_start != 0)
   {
      cstr_start--;
   }

   if (cstr_start > strlen_cstr)
   {
      cstr_start = strlen_cstr;
   }

   if (strlen_cstr - cstr_start < cstr_len)
   {
      cstr_len = strlen_cstr - cstr_start;
   }

   /* define loop size */
   mlen = cstr_len;
   if (mlen > str_len)
   {
      mlen = str_len;
   }

   /* compare substrings */
   for (i = 0; i < mlen; i++)
   {
      if (str->data[i+str_start] < cstr[i+cstr_start])
      {
         return -1;
      }
      if (str->data[i+str_start] > cstr[i+cstr_start])
      {
         return 1;
      }
   }

   /* for equal substrings the larger one wins */
   if (str_len < cstr_len)
   {
      return -1;
   }
   else if (str_len > cstr_len)
   {
      return 1;
   }

   return 0;
}


/** @brief Compare a string with a character string.
 *
 * Compares two strings and returns an integer indicating the differences. A
 * return value less than zero signals, that @c str is less than @c cstr, a
 * return value of zero signals equal strings and a greater than zero signals,
 * that @c str is greater than @c cstr.
 *
 * @param[in] str String object.
 * @param[in] cstr Character string.
 *
 */
int
str_compare_cstr (const Str* str, const char* cstr)
{
   size_t len = strlen (cstr);
   return istr_compare (str, 0, str->len, cstr, len, 0, len);
}


/** @brief Compare two strings.
 *
 * Compares two strings and returns an integer indicating the differences. A
 * return value less than zero signals, that @c str1 is less than @c str2, a
 * return value of zero signals equal strings and a greater than zero signals,
 * that @c str1 is greater than @c str2.
 *
 * @param[in] str1 String object.
 * @param[in] str2 String object.
 *
 */
int
str_compare_str (const Str* str1, const Str* str2)
{
   return istr_compare (str1,
                        0,
                        str1->len,
                        str2->data,
                        str2->len,
                        0,
                        str2->len);
}


/** @brief Compare two substrings.
 *
 * Compares two substrings and returns an integer indicating the differences.
 * The substrings are fetched from a string object and a character string. 
 * @c str_start, @c cstr_start define the starting points and @c str_len,
 * @c cstr_len the length of the substrings. If one of the parameters is out of
 * bounds, it is automatically set to the maximal possible value. E.g. if
 * @c str_start is larger than the length of the string, it will be set to the
 * string length. Positions start at 1 but 0 is acceptable as an exception for
 * the first position.\n
 * Returns a value less than zero when @c str is less than @c cstr. A return
 * value of zero signals equal strings and a value greater than zero signals,
 * that @c str is greater than @c cstr.
 *
 * @param[in] str String object.
 * @param[in] str_start Start position in the string.
 * @param[in] str_len Length of substring to compare.
 * @param[in] cstr Character string.
 * @param[in] cstr_start Start position in the character string.
 * @param[in] cstr_len Length of substring to compare.
 */
int
str_compare_csubstr (const Str* str, size_t str_start, size_t str_len,
                     const char* cstr, size_t cstr_start, size_t cstr_len)
{
   return istr_compare (str,
                        str_start,
                        str_len,
                        cstr,
                        strlen (cstr),
                        cstr_start,
                        cstr_len);
}


/** @brief Compare two substrings.
 *
 * Compares two substrings and returns an integer indicating the differences.
 * The substrings are fetched from two string objects. @c str1_start,
 * @c str2_start define the starting points and @c str1_len, @c str2_len the
 * length of the substrings. If one of the parameters is out of bounds, it is
 * automatically set to the maximal possible value. E.g. if @c str1_start is
 * larger than the length of the string, it will be set to the string's
 * length. Positions start at 1 but 0 is acceptable as an exception for
 * the first position.\n
 * Returns a value less than zero when @c str1 is less than @c str2. A return
 * value of zero signals equal strings and a value greater than zero signals,
 * that @c str1 is greater than @c str2.
 *
 * @param[in] str1 String object.
 * @param[in] str1_start Start position in the string.
 * @param[in] str1_len Length of substring to compare.
 * @param[in] str2 String object.
 * @param[in] str2_start Start position in the character string.
 * @param[in] str2_len Length of substring to compare.
 */
int
str_compare_substr (const Str* str1, size_t str1_start, size_t str1_len,
                    const Str* str2, size_t str2_start, size_t str2_len)
{
   return istr_compare (str1,
                        str1_start,
                        str1_len,
                        str2->data,
                        str2->len,
                        str2_start,
                        str2_len);
}
