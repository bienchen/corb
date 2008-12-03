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
 *  @file libcrbbasic/gfile.c
 *
 *  @brief Generic file handling
 *
 *  Module: gfile
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-11-01
 *
 *
 *  Revision History:
 *         - 2008Nov01 bienert: created
 *
 */


#include <config.h>
#include <stddef.h>
#include <assert.h>
#include "errormsg.h"
#include "memmgr.h"
#include "str.h"
#include "mprintf.h"
#include "gfile.h"


struct GFile {
   GFileType type;
   union {
      FILE* uc; 
   } fileptr;
   Str* path;
};

static void
gfile_delete (GFile* this)
{
   if (this != NULL)
   {
      str_delete (this->path);
      XFREE(this);
   }
}


/** @brief Get the type of a file.
 *
 * Retrive the type of a file. The file name may contain the whole path\n
 * Returns GFILE_UNCOMPRESSED for uncompressed files or if @c file is @c NULL.
 *
 * @param[in] file File name.
 * @param[in] length Length of file name.
 */
GFileType
gfile_determine_type (const char* file __attribute__((unused)),
                      unsigned long length __attribute__ ((unused)))
{
   return GFILE_UNCOMPRESSED;
}

/** @brief Open a file.
 *
 * Opens an uncompressed, bzip2- or gzip-compressed file. @c type provides the
 * possibility to force a file to be opened as a certain type. Usually @c type
 * should be passed as @c GFILE_VOID. In that case @c gfile_open determines the
 * type by the file extension. @c mode accepts the same modes as @c fopen. If
 * compiled with enabled memory checking, @c file and @c line should point to
 * the position where the function was called. Both parameters are
 * automatically set by using the macro @c GFILE_OPEN.\n
 * Returns a valid file pointer or @c NULL on error.
 *
 * @param[in] filepath File name.
 * @param[in] length Length of file name.
 * @param[in] type Type of file.
 * @param[in] mode Mode to open file with.
 */
GFile*
gfile_open (const char* filepath,
            const unsigned long length,
            const GFileType type,
            const char* mode,
            const char* file, const int line)
{
   GFile* gfile = XOBJ_MALLOC(sizeof (GFile), file, line);
   
   assert (filepath);
   assert (mode);

   if (gfile != NULL)
   {
      gfile->path = str_new_cstr (filepath, file, line);
      if (gfile->path == NULL)
      {
         gfile_delete (gfile);
         gfile = NULL;
      }
   }

   if (gfile != NULL)
   {
      gfile->type = type;

      /* check for file type */
      if (gfile->type == GFILE_VOID)
      {
         gfile->type = gfile_determine_type (filepath, length);
      }
      
      /* open file */
      switch (gfile->type)
      {
         case GFILE_UNCOMPRESSED:
            gfile->fileptr.uc = fopen (filepath, mode);
            if (gfile->fileptr.uc == NULL)
            {
               gfile_delete (gfile);
               gfile = NULL;
            }
            break;            
         default:
            THROW_ERROR_MSG ("Opening file \"%s\" failed: Unknown file type",
                             filepath);
            gfile_delete (gfile);
            return NULL;
      }
   }

   if (gfile == NULL)
   {
      THROW_ERROR_MSG ("Opening file \"%s\" failed:", filepath);
   }

   return gfile;
}

/** @brief Close a file.
 *
 * Close a generic file and free all occupied ressources.\n
 * Returns @c EOF on error, 0 else.
 *
 * @param[in] gfile GFile pointer
 */
int
gfile_close (GFile* gfile)
{
   int ret_val = 0;

   if (gfile == NULL)
   {
      return ret_val;
   }

   switch (gfile->type)
   {
      case GFILE_UNCOMPRESSED:
         ret_val = fclose (gfile->fileptr.uc);
         break;            
      default:
         THROW_ERROR_MSG ("Closing file \"%s\" failed: Unknown file type",
                          str_get (gfile->path));
         gfile_delete (gfile);
         return EOF;
   }

   if (ret_val == EOF)
   {
      THROW_ERROR_MSG ("Closing file \"%s\" failed:", str_get (gfile->path));
   }

   gfile_delete (gfile);

   return ret_val;
}

/** @brief Read from file.
 *
 * Reads @c nobj items of size @c size from stream @c stream into array @c
 * ptr. Unlike the error handling of ANSI C @c fread, we provide an error
 * variable @c error defaulting to 0. Possible values on error are:
 * GFILE_UNKNOWN_TYPE for file of unknown compression state, ...\n
 * Returns the number of items read. Values smaller than @c nobj indicate
 * end-of-file or error.
 *
 * @param[in] gfile GFile pointer
 */
enum gfile_errors {
   GFILE_UNKNOWN_TYPE = 1,
   GFILE_READ_ERROR,
   GFILE_MEM_ERROR
};

size_t
gfile_fread (int* error, void* ptr, size_t size, size_t nobj, GFile *stream)
{
   size_t obj_read = EOF;

   assert (error);
   assert (ptr);
   assert (stream);

   switch (stream->type)
   {
      case GFILE_UNCOMPRESSED:
         obj_read = fread(ptr, size, nobj, stream->fileptr.uc);
         if (obj_read != nobj)
         {
            /* check for error */
            if (ferror(stream->fileptr.uc))
            {
               THROW_ERROR_MSG ("Reading from file \"%s\" failed:",
                                str_get (stream->path));
               *error = GFILE_READ_ERROR;
            }
         }
         break;
      default:
         THROW_ERROR_MSG ("Reading from file \"%s\" failed: Unknown file type",
                          str_get (stream->path));
         *error = GFILE_UNKNOWN_TYPE;
   }

   return obj_read;
}

/** @brief Read from file.
 *
 * Does the same as @c gfile_fread() on a character array BUT translates all
 * tabs into whitespaces.\n
 * For return values/ errors please refer to @c gfile_fread().
 *
 * @param[in] gfile GFile pointer
 */
size_t
gfile_fread_tab (int* error, char* ptr, size_t size, size_t nobj, GFile *stream)
{
   size_t n_read;
   size_t i;

   n_read = gfile_fread (error, ptr, size, nobj, stream);

   for (i = 0; i < n_read; i++)
   {
      if (ptr[i] == '\t')
      {
         ptr[i] = ' ';
      }
   }

   return n_read;
}

/** @brief Read from file.
 *
 * Basically does the same as @c gfile_fread() with paying attention to
 * comments. Only comments introduced by single character @c symbol are
 * recognised. This enables shell style comments, not blocks like in C.
 * Everything after @c symbol will be ignored. This means, if @c nobj is
 * greater than 1, we continue reading on the line following the comment.
 * Thereby the newline is translated into a whitepace.\n
 * Additionally all tabs are translated into whitespaces.\n
 * For return values/ errors please refer to @c gfile_fread().
 *
 */
size_t
gfile_fread_comment (const char comment,
                     int* error,
                     char* ptr,
                     const size_t nobj,
                     GFile* stream)
{
   size_t i, j;                 /* iterator */
   size_t cr;                   /* size of current read */

   cr = gfile_fread (error, ptr, sizeof (*ptr), nobj, stream);

   for (i = 0; i < cr; i++)
   {
      if (ptr[i] != '\n')
         mfprintf (stderr, "%lu: |%c|\n", (unsigned long)i, ptr[i]);
      else
         mfprintf (stderr, "%lu: |\\n|\n", (unsigned long)i);
   }

   for (i = 0; i < cr; i++)
   {
      if (ptr[i] == comment)
      {
         ptr[i] = '\n';
         mfprintf (stderr, "Found comment at %lu\n", (unsigned long)i);

         /* search for \n */
         /* step 1: Search current buffer */
         j = i + 1;
         while ((j < cr) && (ptr[j] != '\n'))
         {
            j++;
         }
         mfprintf (stderr, "Found \\n at %lu\n", (unsigned long)j);

         /* compress buffer */
         /* look out for new comments after the \n */
         j++;
         for (; j < cr; j++)
         {
            i++;

           mfprintf (stderr, "Swap ");
            if (ptr[i] != '\n')
               mfprintf (stderr, "|%c| ", ptr[i]);
            else
               mfprintf (stderr, "|\\n|");
            mfprintf (stderr, " with ");
            if (ptr[j] != '\n')
               mfprintf (stderr, "|%c|", ptr[j]);
            else
               mfprintf (stderr, "|\\n|");
            mfprintf (stderr, "\n");

            ptr[i] = ptr[j];
         }
         cr = i + 1;

         /* step 2: Search file */
         mfprintf (stderr, "cr: %lu j: %lu i: %lu\n", (unsigned long) cr, (unsigned long) j, (unsigned long) i);

         /* if j read to the end of buffer, we have to search file for '\n' */

      }
   }

   /* 4: delete from current buffer to \n */
   /* 5: read from file */
   /* 6: delete tabs */

   return cr;
}
/* size_t */
/* gfile_fread_comment (const char comment __attribute__((unused)), */
/*                      int* error, */
/*                      char* ptr, */
/*                      const size_t nobj, */
/*                      GFile* stream) */
/* { */
/*    size_t i; */
/*    size_t n_read = 0; */
/*    size_t cr; */

/*    while ((n_read < nobj) && (!(*error))) */
/*    { */
/*       cr = gfile_fread (error, */
/*                         (ptr + n_read), */
/*                         sizeof (*ptr), */
/*                         (nobj - n_read), */
/*                         stream); */
   
/*       for (i = 0; i < cr; i++) */
/*       { */
/*          if (ptr[i] == '\t') */
/*          { */
/*             ptr[i] = ' '; */
/*          } */
/*          else if (ptr[i] == comment) */
/*          { */
/*             /\* skip to '\n' *\/ */
/*             /\* step 1: Search in current buffer *\/ */
            

/*             /\* step 2: Search in file *\/ */
/*             cr = 1; */
/*             while ((cr == 1) && (ptr[i] != '\n') && (!(*error))) */
/*             { */
/*                cr = gfile_fread (error, */
/*                                  (ptr + n_read), */
/*                                  sizeof (*ptr), */
/*                                  1, */
/*                                  stream); */
/*             } */
/*             ptr[n_read] = '\n'; */
            
/*             if ((cr == 0) || (*error)) */
/*             { */
/*                return n_read; */
/*             } */
/*          } */

/*          n_read++; */
/*       } */
/*    } */

/*    return n_read; */
/* } */

/** @brief Read from file.
 *
 * Basically does the same as @c gfile_fread() with paying attention to
 * comments. For this function we assume shell style comments, introduced by
 * "#". Everything after this symbol up to '\n' is ignored. If we encounter
 * @c EOF after "#", we just stop there.\n
 * Additionally all tabs are translated into whitespaces.\n
 * For return values/ errors please refer to @c gfile_fread().
 *
 */
size_t
gfile_fread_comment_sh (int* error,
                        char* ptr,
                        size_t size,
                        size_t nobj,
                        GFile *stream)
{
   size_t n_read;
   unsigned long i;

   n_read = gfile_fread (error, ptr, size, nobj, stream);

   if (!(*error))
   {
      for (i = 0; i < n_read; i++)
      {
         if (ptr[i] == '\t')
         {
            ptr[i] = ' ';
         }
      }
   }

   return n_read;
}

/** @brief Read a line form file.
 *
 * Basic function for reading lines from text files. The data read is stored in
 * @c ptr. The size of @c ptr is determined by @c size. If the line to be read
 * is longer than @c size, @c ptr is reallocated to fit the line and @c size is
 * updated. If the line read is shorter than @c size @c ptr is not changed in
 * size. Each line read is terminated the null character while any newline
 * character is chopped of. Since a text file can contain null characters you
 * should not iterate over a line until '\0' occurs in @c ptr.\n
 * @c error is used to signal read failures. For successful reading @c error is
 * 0. Has only to be checked if 0 is returned, otherwise @c error is 0.\n
 * Returns 0 if the end of file is reached or on error. Otherwise the length of
 * the line read. If 0 is returned, @c ptr is unchanged.
 *
 * @param[out] error Container for error values.
 */
static __inline__ int
s_gfile_store_char (const char c,
                    char** buf,
                    size_t* size,
                    const unsigned long pos)
{
   /* care for enough memory */
   /* we assume that this function is used iteratively, therefore before
      pos > size, we ran through pos == size and we do not have to determine
      the right size to allocate */
   if (pos == *size)
   {
      *size += 78;        /* just add a usual line width */
      *buf = XREALLOC (*buf, sizeof (**buf) * (*size));
      if (*buf == NULL)
      {
         return GFILE_MEM_ERROR;
      }
   }
   
   /* store symbol */
   (*buf)[pos] = c;

   return 0;
}

unsigned long
gfile_getline (int* error, char** buf, size_t* size, GFile* stream)
{
   int c;                       /* current character */
   unsigned long length = 0;    /* length of current line */

   assert (error);
   assert (buf);
   assert (size);
   assert (stream);
   assert (stream->fileptr.uc);

   if (stream->type == GFILE_UNCOMPRESSED)
   {
      while ((c = fgetc (stream->fileptr.uc)) != EOF)
      {
         *error = s_gfile_store_char (c, buf, size, length);
         if (*error)
         {
            THROW_ERROR_MSG ("Reading file \"%s\" stopped.",
                             str_get (stream->path));
            return 0;
         }

         /* check line end criteria */
         if ((*buf)[length] == '\n')
         {
            (*buf)[length] = '\0';
            return length + 1;
         }

         length++;
      }
      
      /* handle errors of fgetc */
      if (ferror(stream->fileptr.uc))
      {
         THROW_ERROR_MSG ("Reading from file \"%s\" failed:",
                          str_get (stream->path));
         *error = GFILE_READ_ERROR;
         return 0;
      }

      *error = s_gfile_store_char ('\0', buf, size, length);
      if (*error)
      {
         THROW_ERROR_MSG ("Reading file \"%s\" stopped.",
                          str_get (stream->path));
         return 0;
      }
   }
   else
   {
      THROW_ERROR_MSG ("Reading from file \"%s\" failed: Unknown file type",
                       str_get (stream->path));
      *error = GFILE_UNKNOWN_TYPE;
      return 0;
   }

   return length;
}
