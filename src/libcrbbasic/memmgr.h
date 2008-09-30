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
 *  @file libcrbbasic/memmgr.h
 *
 *  @brief Replacements for malloc, calloc, realloc and free.
 *
 *  Module: memmgr
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-01-18
 *
 *
 *  Revision History:
 *         - 2008Jan18 bienert: created
 *
 */

#include <config.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MEMMGR_H
#define MEMMGR_H


enum memman_retvals{
   ERR_MM_ALLOC = 1   /* Failed allocation of memory for the memory manager */
};

void*
xmalloc (const size_t, const char*, const int)
#if defined(__GNUC__) && (__GNUC__ > 2)
     __attribute__((malloc))
#endif /* __GNUC__ */
;

void*
checked_xmalloc (const size_t, const char*, const int)
#if defined(__GNUC__) && (__GNUC__ > 2)
     __attribute__((malloc))
#endif /* __GNUC__ */
;

void*
xrealloc (void*, const size_t, const char*, const int)
#if defined(__GNUC__) && (__GNUC__ > 2)
     __attribute__((malloc))
#endif /* __GNUC__ */
;

void*
checked_xrealloc (void*, const size_t, const char*, const int)
#if defined(__GNUC__) && (__GNUC__ > 2)
     __attribute__((malloc))
#endif /* __GNUC__ */
;

void*
xcalloc (const size_t, const size_t, const char*, const int)
#if defined(__GNUC__) && (__GNUC__ > 2)
     __attribute__((malloc))
#endif /* __GNUC__ */
;

void*
checked_xcalloc (const size_t, const size_t, const char*, const int)
#if defined(__GNUC__) && (__GNUC__ > 2)
     __attribute__((malloc))
#endif /* __GNUC__ */
;

void**
xmalloc_2d (const size_t, const size_t, const size_t, const char*, const int)
#if defined(__GNUC__) && (__GNUC__ > 2)
     __attribute__((malloc))
#endif /* __GNUC__ */
;

void*
xmalloc_rnd (const size_t, const size_t, const size_t*, const char*, const int)
#if defined(__GNUC__) && (__GNUC__ > 2)
     __attribute__((malloc))
#endif /* __GNUC__ */
;

void*
xmalloc_nd (const size_t, const char*, const int, const size_t, ...)
#if defined(__GNUC__) && (__GNUC__ > 2)
     __attribute__((malloc))
#endif /* __GNUC__ */
;

void
xfree (void* ptr);

void
checked_xfree (void* ptr);

void
xfree_2d (void** matrix);

void
xfree_nd (const size_t, void**);

void
free_memory_manager (void);

/* Calling macros */
#ifdef MEMCHECK
#define XMALLOC(SIZE)                \
   checked_xmalloc(SIZE, __FILE__, __LINE__)
#define XOBJ_MALLOC(SIZE, F, L)                 \
   checked_xmalloc((SIZE), F, L)
#define XREALLOC(PTR,SIZE)                      \
   checked_xrealloc(PTR,SIZE, __FILE__, __LINE__)
#define XOBJ_REALLOC(PTR,SIZE,F,L)            \
   checked_xrealloc(PTR,SIZE, F, L)   
#define XCALLOC(N,SIZE)                      \
   checked_xcalloc(N,SIZE, __FILE__, __LINE__)
#define XOBJ_CALLOC(N,SIZE,F,L)                  \
   checked_xcalloc(N,SIZE, F, L)
#define XFREE(PTR)                   \
   checked_xfree(PTR)
#define FREE_MEMORY_MANAGER                      \
   free_memory_manager()
#else
#define XMALLOC(SIZE)                          \
   xmalloc(SIZE, __FILE__, __LINE__)
#define XOBJ_MALLOC(SIZE, F, L)                 \
   xmalloc(SIZE, F, L)
#define XREALLOC(PTR,SIZE)                      \
   xrealloc(PTR,SIZE, __FILE__, __LINE__)
#define XOBJ_REALLOC(PTR,SIZE,F,L)            \
   xrealloc(PTR,SIZE, F, L)   
#define XCALLOC(N,SIZE)                      \
   xcalloc(N,SIZE, __FILE__, __LINE__)
#define XOBJ_CALLOC(N,SIZE,F,L)                  \
   xcalloc(N,SIZE, F, L)
#define XFREE(PTR)          \
   xfree(PTR)
#define FREE_MEMORY_MANAGER 
#endif /* MEMCHECK */

#define XMALLOC_2D(ROWS,COLS,SIZE) \
   xmalloc_2d (ROWS, COLS, SIZE, __FILE__, __LINE__)

#define XOBJ_MALLOC_2D(ROWS,COLS,SIZE,F,L)               \
   xmalloc_2d (ROWS, COLS, SIZE, F, L)

#define XMALLOC_RND(SIZE, N, DIM) \
   xmalloc_rnd (SIZE, N, DIM, __FILE__, __LINE__)

#define XOBJ_MALLOC_RND(SIZE, N, DIM, F, L)            \
   xmalloc_rnd (SIZE, N, DIM, F, L)

#define XMALLOC_ND(SIZE,DIM,...) \
   xmalloc_nd (SIZE, __FILE__, __LINE__, DIM, __VA_ARGS__)

#define XFREE_2D(PTRPTR) \
   xfree_2d (PTRPTR)

#define XFREE_ND(N, PTRPTR)                      \
   xfree_nd (N, PTRPTR)

#endif /* MEMMGR_H */

#ifdef __cplusplus
}
#endif
