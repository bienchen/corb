/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbasic/str.h
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


#ifdef __cplusplus
extern "C" {
#endif

#include "inc_bool.h"

#ifndef STR_H
#define STR_H


/* error numbers */
enum str_retvals{
   ERR_STR_ALLOC = 1,      /* (re)allocation problems */
   ERR_STR_OOR             /* out of range of string capacity */
};


typedef struct Str Str;


/**********************   Constructors and destructors   **********************/

extern Str*
str_new (const char*, const int);

#define STR_NEW str_new (__FILE__, __LINE__)

extern Str*
str_new_cstr (const char*, const char*, const int);

#define STR_NEW_CSTR(CSTR) str_new_cstr (CSTR, __FILE__, __LINE__)

extern Str*
str_new_str (const Str*, const char*, const int);

#define STR_NEW_STR(STR) str_new_str (STR, __FILE__, __LINE__)

extern void
str_delete (Str*);


/********************************   Altering   ********************************/

extern int
str_cpy (Str*, const Str*);

extern int
str_set (Str*, const char*);

extern int
str_at (Str*, const unsigned long, const char);

extern void
str_set_i (Str*, const size_t, const char);

extern int
str_append_str (Str*, const Str*);

extern int
str_append_cstr (Str*, const char*);

extern int
str_assign_substr (Str*, const Str*, const size_t, const size_t);

extern int
str_assign_csubstr (Str*, const char*, const size_t, const size_t);

extern void
str_clear (Str*);


/*********************************   Access   *********************************/

extern const char*
str_get (const Str*);

extern char
str_get_i (const Str*, const unsigned long);


/******************* Size *******************/

extern unsigned long
str_length (const Str*);

extern size_t
str_capacity (const Str*);

extern bool
str_empty (const Str*);

extern int
str_resize (Str*, const size_t, const char);


/******************* Searching *******************/

extern unsigned long
str_find_c (const Str*, const char);

extern unsigned long
str_rfind_c (const Str*, const char);

extern unsigned long
str_find_cstr (const Str*, const char*);

extern unsigned long
str_find_str (const Str*, const Str*);

extern unsigned long
str_rfind_cstr (const Str*, const char*);

extern unsigned long
str_rfind_str (const Str*, const Str*);

extern unsigned long
str_find_first_of_cstr (const Str*, const char*, const unsigned long);

extern unsigned long
str_find_first_of_str (const Str*, const Str*, const unsigned long);

extern unsigned long
str_find_first_not_of_cstr (const Str*, const char*, const unsigned long);

extern unsigned long
str_find_first_not_of_str (const Str*, const Str*, const unsigned long);

extern unsigned long
str_find_last_of_cstr (const Str*, const char*, const unsigned long);

extern unsigned long
str_find_last_of_str (const Str*, const Str*, const unsigned long);

extern unsigned long
str_find_last_not_of_cstr (const Str*, const char*, const unsigned long);

extern unsigned long
str_find_last_not_of_str (const Str*, const Str*, const unsigned long);

/******************* Comparison *******************/
extern int
str_compare_cstr (const Str*, const char*);

extern int
str_compare_str (const Str*, const Str*);

extern int
str_compare_csubstr (const Str*, size_t, size_t, const char*, size_t, size_t);

extern int
str_compare_substr (const Str*, size_t, size_t, const Str*, size_t, size_t);

#endif /* STR_H */

#ifdef __cplusplus
}
#endif
