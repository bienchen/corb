/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbasic/memmgr.c
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
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include "inc_bool.h"
#include "inc_pthr.h"
#include "inc_strg.h"
#include "mprintf.h"
#include "errormsg.h"
#include "gcckywrd.h"
#include "memmgr.h"


/** @brief Information about a memory block.
 *
 */
typedef struct
{
      unsigned long addr; /**< address of the block */
      size_t size;        /**< size of the block */
      char* file;         /**< file where block was allocated */
      int line;           /**< line in @c file where block was allocated */
} MemBlk;

/** @brief List of memory blocks.
 *
 */
typedef struct
{
      unsigned long next_free;  /**< next available position */
      unsigned long size;       /**< size of the list */
      MemBlk* lst;              /**< actual list */
} MemBlkLst;

/** @brief Takes all components of the memory management system.
 *
 */
typedef struct
{
      MemBlkLst* addr_tbl;             /**< hash table for the pointers */
      unsigned long tbl_size;          /**< size of the hash table */
#ifdef HAVE_PTHREAD
      pthread_mutex_t mm_mutex;         /**< thread safety */
#endif
} MemManSys;

/** @brief Memory manager.
 *
 */
static MemManSys __memman = {NULL,     /* address table */
                             1009,     /* size of the table (prime number) */
#ifdef HAVE_PTHREAD
                             PTHREAD_MUTEX_INITIALIZER
#endif
                            }; /* init mutex */


/** @brief Calculates a hash value for an address.
 *
 * Calulates the hash value for a given pointer in @c unsigned long format.
 * Immediately applies the modulo operator with the size of the address table
 * to the result.\n
 * Returns the index where the memory block should be stored.
 *
 * @param[in] ptr The pointer.
 */
static __inline__ unsigned long
hash_address(unsigned long ptr)
{
  ptr = (~ptr) + (ptr << 21); // ptr = (ptr << 21) - ptr - 1;
  ptr = ptr ^ (ptr >> 24);
  ptr = (ptr + (ptr << 3)) + (ptr << 8); // ptr * 265
  ptr = ptr ^ (ptr >> 14);
  ptr = (ptr + (ptr << 2)) + (ptr << 4); // ptr * 21
  ptr = ptr ^ (ptr >> 28);
  ptr = ptr + (ptr << 31);
  return ptr%__memman.tbl_size;
}


/** @brief Search a memory block.
 *
 * Returns a pointer to a block or @c Null.
 *
 * @param[in] hashkey The hash value to be searched.
 * @param[in] ptr Address of the block to be deleted.
 */
static __inline__ MemBlk*
mblist_search_block (const unsigned long hashkey, const unsigned long ptr)
{
   unsigned long i;

   /* check that list exists */
   if (__memman.addr_tbl[hashkey].lst != NULL)
   {
      /* search list */
      for (i = 0; i < __memman.addr_tbl[hashkey].next_free; i++)
      {
         if (__memman.addr_tbl[hashkey].lst[i].addr == ptr)
         {
            return &(__memman.addr_tbl[hashkey].lst[i]);
         }
      }

      return NULL;
   }
   
   return NULL;
}


/** @brief Delete a memory block from the hash table.
 *
 * @param[in] hashkey Hash value of the block to be removed.
 * @param[in] ptr Address of the block to be removed from the table.
 * @param[in] file Calling file.
 * @param[in] line Called from line.
 */
static __inline__ void
mblist_delete_block (const unsigned long hashkey,
                     const unsigned long ptr)
{
   unsigned long i;

   /* check that list exists */
   if (__memman.addr_tbl[hashkey].lst != NULL)
   {
      /* search list */
      for (i = 0; i < __memman.addr_tbl[hashkey].next_free; i++)
      {
         if (__memman.addr_tbl[hashkey].lst[i].addr == ptr)
         {
            if (__memman.addr_tbl[hashkey].next_free > 0)
            {
               xfree (__memman.addr_tbl[hashkey].lst[i].file);
               __memman.addr_tbl[hashkey].lst[i]
    = __memman.addr_tbl[hashkey].lst[__memman.addr_tbl[hashkey].next_free - 1];
               __memman.addr_tbl[hashkey].next_free--;
            }
            return;
         }
      }
   }

   /* block not found! */
   THROW_WARN_MSG ("Request to remove unknown block from memory bookkeeper "
                   "(%p).", (void*) ptr);
}


/** @brief Add a memory block to the memory block list.
 *
 * Add a pointer to the management system.\n
 * Returns 0 on success, ERR_MM_ALLOC, on failure. ERR_MM_ALLOC means that the
 * allocation of memory for the table itself failed, so a serious problem has
 * occured.
 *
 * @param[in] hashkey Index of the list in the hash table.
 * @param[in] ptr The pointer to be stored.
 * @param[in] size Size of the block.
 * @param[in] file File calling for allocation on
 * @param[in] line Line.
 */
static __inline__ int
mblist_add_block (const unsigned long hashkey,
                  const unsigned long ptr,
                  const size_t size,
                  const char* file,
                  const int line)
{
   /* do we need allocation of list first? */
   if (__memman.addr_tbl[hashkey].next_free == __memman.addr_tbl[hashkey].size)
   {
      __memman.addr_tbl[hashkey].size += 100;
      __memman.addr_tbl[hashkey].lst = xrealloc(__memman.addr_tbl[hashkey].lst,
                                                __memman.addr_tbl[hashkey].size
                                                * sizeof (MemBlk),
                                                __FILE__,
                                                __LINE__);
      if (__memman.addr_tbl[hashkey].lst == NULL)
      {
         THROW_ERROR_MSG ("Memory bookkeeping failed for allocating call from "
                          "\"%s\"at line %d", file, line);
         return ERR_MM_ALLOC;
      }
   }

   /* store block */
   __memman.addr_tbl[hashkey].lst[__memman.addr_tbl[hashkey].next_free].addr
      = ptr;
   __memman.addr_tbl[hashkey].lst[__memman.addr_tbl[hashkey].next_free].size
      = size;
   __memman.addr_tbl[hashkey].lst[__memman.addr_tbl[hashkey].next_free].line
      = line;
   __memman.addr_tbl[hashkey].lst[__memman.addr_tbl[hashkey].next_free].file
      = xmalloc ((strlen (file) + 1), __FILE__, __LINE__);
   strcpy (
      __memman.addr_tbl[hashkey].lst[__memman.addr_tbl[hashkey].next_free].file,
      file);

   __memman.addr_tbl[hashkey].next_free++;

   return 0;
}


/** @brief Add a memory block to the memory block list.
 *
 * Add a pointer to the management system. On problems, this function just
 * writes warnings to @c stderr.\n
 * Returns @c ERR_MM_ALLOC on problems allocating memory for the hash table.
 *
 * @param[in] ptr The pointer to be stored.
 * @param[in] size Size of the block.
 * @param[in] file File calling for allocation on.
 * @param[in] line Line.
 */
static int
mm_sys_add_block (const void* ptr,
                  const size_t size,
                  const char* file,
                  const int line)
{
   unsigned long i;
   unsigned long hashkey;
   unsigned long int_ptr = (unsigned long) ptr;
   MemBlk* block;

   assert (file);

   if (ptr == NULL)
   {
      call_warn_msgr (file, line,
                      "Request to add an empty block to the memory "
                      "bookkeeper.");
   }

   /* do we have to allocate the hash table first? */
   if (__memman.addr_tbl == NULL)
   {
      __memman.addr_tbl = xmalloc (__memman.tbl_size * sizeof (MemBlkLst),
                                   __FILE__,
                                   __LINE__);
      if (__memman.addr_tbl == NULL)
      {
         THROW_WARN_MSG ("Memory bookkeeping failed for allocating call from "
                          "\"%s\"at line %d", file, line);
         return ERR_MM_ALLOC;
      }

      /* init lists */
      for (i = 0; i < __memman.tbl_size; i++)
      {
         __memman.addr_tbl[i].lst = NULL;
         __memman.addr_tbl[i].size = __memman.addr_tbl[i].next_free = 0;
      }
   }

   hashkey = hash_address (int_ptr);

   /* check whether ptr already stored */
   block = mblist_search_block (hashkey, int_ptr);
   if (block != NULL)
   {
      call_warn_msgr (file, line,
                      "Request to add block %p to memory bookkeeper, already "
                      "allocated in \"%s\" at line %d.",
                      (void*) block->addr,
                      block->file,
                      block->line);
   }

   return mblist_add_block (hashkey, int_ptr, size, file, line);
}


/** @brief Delete a memory block from the memory address hash table.
 *
 * @param[in] ptr The pointer of the block to be deleted.
 * @param[in] file Calling file.
 * @param[in] line Called from line.
 */
static void
mm_sys_delete_block (const void* ptr)
{
   unsigned long int_ptr = (unsigned long) ptr;
   unsigned long hashkey = hash_address (int_ptr);
   if (__memman.addr_tbl != NULL)
   {
      mblist_delete_block (hashkey, int_ptr);
   }
}


/** @brief Final check of the memory manager.
 *
 * This function checks the memory hash table for non-freed blocks and frees
 * the table itself. To be called with @c FREE_MEMORY_MANAGER.\n
 */
void
free_memory_manager (void)
{
   unsigned long i;
   unsigned long j;

#ifdef HAVE_PTHREAD
   pthread_mutex_lock (&__memman.mm_mutex);
#endif

   if (__memman.addr_tbl != NULL)
   {
      for (i = 0; i < __memman.tbl_size; i++)
      {
         if (__memman.addr_tbl[i].lst != NULL)
         {
            for (j = 0; j < __memman.addr_tbl[i].next_free; j++)
            {
               THROW_WARN_MSG ("Memory block of size %lu not freed, allocated "
                               "at \"%s\" line %d.",
                               (unsigned long) __memman.addr_tbl[i].lst[j].size,
                               __memman.addr_tbl[i].lst[j].file,
                               __memman.addr_tbl[i].lst[j].line);
               xfree (__memman.addr_tbl[i].lst[j].file);
            }
            xfree (__memman.addr_tbl[i].lst);
            __memman.addr_tbl[i].lst = NULL;
         }
      }
   }

   xfree (__memman.addr_tbl);
   __memman.addr_tbl = NULL;

#ifdef HAVE_PTHREAD
   pthread_mutex_unlock (&__memman.mm_mutex);
#endif

#ifdef HAVE_PTHREAD
   pthread_mutex_destroy(&__memman.mm_mutex);
#endif
}


/** @brief Wrapper function for malloc.
 *
 * Do not call this directly. Use the macro @c XMALLOC(size_t size) instead.
 * This will asure the memory manager to work properly and set the calling file
 * and line for you.\n
 * Returns NULL on failure, a pointer to the allocated memory else.
 *
 * @param[in] size Number of bytes to be allocated.
 * @param[in] file Calling file. Must not be @c NULL.
 * @param[in] line Called from line.
 */
void*
xmalloc (const size_t size, const char* file, const int line)
{
   void* ptr;

   errno = 0;  /* following recommendation form the man page */
   ptr = malloc (size);
   
   if (ptr == NULL)
   {
      /* examining if malloc set an error */
      if (errno != 0)
      {
         call_error_msgr (file, line,
                          "Memory allocation of %lu bytes failed:",
                          (unsigned long) size);
      }
      else
      {
         call_error_msgr (file, line,
                          "Memory allocation of %lu bytes failed.",
                          (unsigned long) size);
      }

      return ptr;
   }

   return ptr;
}


/** @brief Memory checking wrapper function for malloc.
 *
 * Do not call this directly. Use the macro @c XMALLOC(size_t size) instead.
 * This will asure the memory manager to work properly and set the calling file
 * and line for you.\n
 * Returns NULL on failure, a pointer to the allocated memory else.
 *
 * @param[in] size Number of bytes to be allocated.
 * @param[in] file Calling file. Must not be @c NULL.
 * @param[in] line Called from line.
 */
void*
checked_xmalloc (const size_t size, const char* file, const int line)
{
   void* ptr;

#ifdef HAVE_PTHREAD
   pthread_mutex_lock (&__memman.mm_mutex);
#endif


   ptr = xmalloc (size, file, line);
   
   if (ptr != NULL)
   {
      if (mm_sys_add_block (ptr, size, file, line) == ERR_MM_ALLOC)
      {
         /* serious problem at allocating memory for the memory manager,
            signal error */
         xfree (ptr);

#ifdef HAVE_PTHREAD
         pthread_mutex_unlock (&__memman.mm_mutex);
#endif
         return NULL;
      }
   }

#ifdef HAVE_PTHREAD
   pthread_mutex_unlock (&__memman.mm_mutex);
#endif
   return ptr;
}


/** @brief Wrapper function for realloc.
 * 
 * Do not call this directly. Use the macro @c XREALLOC(void* ptr, size_t size)
 * instead. This will asure the memory manager to work properly and set the
 * calling file and line for you.\n
 * Returns NULL on failure, a pointer to the allocated memory else.
 *
 * @param[in] ptr  Pointer to already allocated memory.
 * @param[in] size Number of bytes to allocate.
 * @param[in] file Calling file. May be @c NULL.
 * @param[in] line Called from line.
 */
void*
xrealloc (void* ptr, const size_t size, const char* file, const int line)
{
    void* new_ptr;

    if (size == 0)
    {
        xfree (ptr);
        return NULL;
    }

    errno = 0;
    new_ptr = realloc (ptr, size);

    if (new_ptr == NULL) 
    {
      /* examining if realloc set an error */
      if (errno != 0)
      {
         call_error_msgr (file, line,
                          "Memory reallocation of %lu bytes failed:",
                          (unsigned long) size);
      }
      else
      {
         call_error_msgr (file, line,
                          "Memory reallocation of %lu bytes failed.",
                          (unsigned long) size);
      }

      return new_ptr;
    }

   return new_ptr;
}


/** @brief Memory checking wrapper function for realloc.
 * 
 * Do not call this directly. Use the macro @c XREALLOC(void* ptr, size_t size)
 * instead. This will asure the memory manager to work properly and set the
 * calling file and line for you.\n
 * Returns NULL on failure, a pointer to the allocated memory else.
 *
 * @param[in] ptr  Pointer to already allocated memory.
 * @param[in] size Number of bytes to allocate.
 * @param[in] file Calling file. May be @c NULL.
 * @param[in] line Called from line.
 */
void*
checked_xrealloc (void* ptr,
                  const size_t size,
                  const char* file,
                  const int line)
{
    void* new_ptr;

#ifdef HAVE_PTHREAD
   pthread_mutex_lock (&__memman.mm_mutex);
#endif

    if (size == 0)
    {
        XFREE (ptr);

#ifdef HAVE_PTHREAD
        pthread_mutex_unlock (&__memman.mm_mutex);
#endif
        return NULL;
    }

    if (ptr != NULL)
    {
       mm_sys_delete_block (ptr);
    }

    errno = 0;
    new_ptr = xrealloc (ptr, size, file, line);

    if (new_ptr == NULL) 
    {
      /* examining if realloc set an error */
      if (errno != 0)
      {
         call_error_msgr (file, line,
                          "Memory reallocation of %lu bytes failed:",
                          (unsigned long) size);
      }
      else
      {
         call_error_msgr (file, line,
                          "Memory reallocation of %lu bytes failed.",
                          (unsigned long) size);
      }

#ifdef HAVE_PTHREAD
      pthread_mutex_unlock (&__memman.mm_mutex);
#endif
      return new_ptr;
    }

    mm_sys_add_block (new_ptr, size, file, line);

#ifdef HAVE_PTHREAD
    pthread_mutex_unlock (&__memman.mm_mutex);
#endif
    return new_ptr;
}


/** @brief Wrapper for @c calloc.
 *
 * Do not call this directly. Use the macro
 * @c XCALLOC(size_t nmemb, size_t size, char* file, int line)
 * instead. This will asure the memory manager to work properly and set the
 * calling file and line for you.\n
 * Returns NULL on failure, a pointer to the allocated memory else.
 *
 * @param[in] nmemb Number of elements.
 * @param[in] size Size of one element.
 * @param[in] file Calling file. May be @c NULL.
 * @param[in] line Called from line.
 */
void*
xcalloc (const size_t nmemb,
         const size_t size,
         const char* file,
         const int line)
{
    void *ptr;

    errno = 0;
    ptr = calloc (nmemb, size);
    if (ptr == NULL)
    {
       if (errno != 0)
       {
         call_error_msgr (file, line,
                          "Memory allocation of %lu bytes failed:",
                          (unsigned long) size);
       }
       else
       {
         call_error_msgr (file, line,
                          "Memory allocation of %lu bytes failed.",
                          (unsigned long) size);
       }
    }

    return ptr;
}


/** @brief Memory checking wrapper for @c calloc.
 *
 * Do not call this directly. Use the macro
 * @c XCALLOC(size_t nmemb, size_t size, char* file, int line)
 * instead. This will asure the memory manager to work properly and set the
 * calling file and line for you.\n
 * Returns NULL on failure, a pointer to the allocated memory else.
 *
 * @param[in] nmemb Number of elements.
 * @param[in] size Size of one element.
 * @param[in] file Calling file. May be @c NULL.
 * @param[in] line Called from line.
 */
void*
checked_xcalloc (const size_t nmemb,
                 const size_t size,
                 const char* file,
                 const int line)
{
    void *ptr;

#ifdef HAVE_PTHREAD
    pthread_mutex_lock (&__memman.mm_mutex);
#endif

    ptr = xcalloc (nmemb, size, file, line);

    if (ptr != NULL)
    {
       mm_sys_add_block (ptr, size, file, line);
    }

#ifdef HAVE_PTHREAD
    pthread_mutex_unlock (&__memman.mm_mutex);
#endif
    return ptr;
}


/** @brief Allcoate memory for a 2-dimensional matrix
 *
 * Do not call this directly. Use the macro
 * @c XMALLOC_2D(size_t n_rows, size_t n_cols, size_t size) instead. This will
 * asure the memory manager to work properly and set the calling file and line
 * for you.\n
 * Allocates a matrix with @c n_rows and @c n_cols with cells of @c size bytes
 * size.\n
 * Returns a pointer to the allocated memory on success, @c NULL else.
 *
 * @params[in] n_rows No. of rows of the matrix.
 * @params[in] n_cols No. of columns of the matrix.
 * @params[in] size No. of bytes per entry of the matrix.
 * @params[in] file Calling file.
 * @params[in] line Called from line.
 */
void**
xmalloc_2d (const size_t n_rows,
            const size_t n_cols,
            const size_t size,
            const char* file, const int line)
{
   size_t i;
   void** matrix;

   /* allocate array of pointers */
   matrix = XOBJ_MALLOC (n_rows * sizeof (*matrix), file, line);

   /* mfprintf (stderr, "Size: %lu\n", (unsigned long) size); */

   if (matrix != NULL)
   {
      /* allocate the whole space for matrix cells at once */
      matrix[0] = XOBJ_MALLOC (n_rows * n_cols * size, file, line);

      if (matrix[0] != NULL)
      {
         for (i = 1; i < n_rows; i++)
         {
            /* set start points for each row */
            matrix[i] = ((char*) matrix[i-1]) + (n_cols * size);
         }
      }
   }

   return matrix;
}


/** @brief Wrapper function for free.
 *
 * Do not call this directly. Use the macro @c XFREE(void* ptr) instead. This
 * will asure the memory manager to work properly.
 *
 * @param[in] ptr Address to be freed.
 */
void
xfree (void* ptr)
{
   free (ptr);
}


/** @brief Memory checking wrapper function for free.
 *
 * Do not call this directly. Use the macro @c XFREE(void* ptr) instead. This
 * will asure the memory manager to work properly.
 *
 * @param[in] ptr Address to be freed.
 * @param[in] file Calling file.
 * @param[in] line Called from line.
 */
void
checked_xfree (void* ptr)
{
#ifdef HAVE_PTHREAD
   pthread_mutex_lock (&__memman.mm_mutex);
#endif
   if (ptr != NULL)
   {
      mm_sys_delete_block (ptr);
      xfree (ptr);
   }
#ifdef HAVE_PTHREAD
   pthread_mutex_unlock (&__memman.mm_mutex);
#endif
}

/** brief Free 2D memory.
 *
 * Free the memory held by a 2D matrix. May be called directly. For consistent
 * style (all memory allocation & releasing via macros), we provide
 * XFREE_2D(void**).
 *
 * @params[in] matrix Matrix to be freed
*/
void
xfree_2d (void** matrix)
{
   /*if (matrix != NULL)
     {*/
   XFREE (*matrix);
   XFREE (matrix);
      /*}*/
}
