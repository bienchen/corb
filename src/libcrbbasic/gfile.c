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
#include <errno.h>
#include <assert.h>
#include <ctype.h>
#include "crb_unused.h"
#include "inc_strg.h"
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

static const char* cmp_type_str[] = {"bz", "gz", "bz2"};

#define N_TYPE_STRS 3

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
gfile_get_type (const char* file, unsigned long length)
{
   CRB_UNUSED (file);
   CRB_UNUSED (length);

   return GFILE_UNCOMPRESSED;
}

static unsigned long
s_gfile_ext_from_list (const char* path,
                       const unsigned long length,
                       const char** list,
                       const unsigned long n_entries)
{
   unsigned long i, j;
   unsigned long ext_len;
   
   assert ((length && path) || (!length && !path));
   assert ((n_entries && list) || (!n_entries && !list));

   for (i = 0; i < n_entries; i++)
   {
      ext_len = strlen (list[i]);
      j = ext_len;

      /* path len has to be > ext. len to carry file name and extension*/
      if ((length > ext_len) && (path[length - 1 - ext_len] == '.'))
      {
         while ((j > 0)
                && (   tolower (list[i][j - 1])
                    == tolower (path[length - ext_len + (j - 1)])))
         {
            j--;
         }

         if (j == 0)
         {
            return i;
         }         
      }
   }

   return n_entries;
}

const char*
gfile_get_type_str (const char* path, unsigned long length)
{
   unsigned long entry;

   entry = s_gfile_ext_from_list (path, length, cmp_type_str, N_TYPE_STRS);

   if (entry < N_TYPE_STRS)
   {
      return cmp_type_str[entry];
   }

   return NULL;
}

unsigned long
gfile_ext_from_list (const char* path,
                     const unsigned long length,
                     const char** list,
                     const unsigned long n_entries)
{
   unsigned long ext_len = 0;
   unsigned long entry;

   /* try to find compression type, if any */
   entry = s_gfile_ext_from_list (path, length, cmp_type_str, N_TYPE_STRS);
   if (entry < N_TYPE_STRS)
   {
      ext_len = strlen (cmp_type_str[entry]);
      ext_len++;                /* count the '.' */
   }

   /* now find the index of the 'true' extension */
   entry = s_gfile_ext_from_list (path, length - ext_len, list, n_entries);

   return entry;
}

/** @brief Get the path of a file in a GFile object.
 *
 * Returns the path of a file connected to a GFile object. The path includes the
 * name of the file.
 *
 * @param[in] this GFile object.
 */
Str*
gfile_get_path (const GFile* this)
{
   assert (this);
   assert (this->path);

   return this->path;
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
         gfile->type = gfile_get_type (filepath, length);
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

/** @brief Reset the current file position to the beginning.
 *
 * Sets the file position of a file stream to the beginning. Can be used for
 * rereading a file.\n
 * Returns
 *        @c GFILE_UNKNOWN_TYPE if the type of file is unsupported,
 *        @c GFILE_REWIND_ERROR for problems on rewinding,
 *        0 otherwise.
 *
 * @param[in] stream File to be reseted
 */
int
gfile_rewind (GFile* stream)
{
   assert (stream);
   assert (stream->fileptr.uc);
   assert (stream->path);

   if (stream->type == GFILE_UNCOMPRESSED)
   {
      errno = 0;
      rewind (stream->fileptr.uc);
      if (errno != 0)
      {
         THROW_ERROR_MSG ("Rewinding file \"%s\" failed:",
                       str_get (stream->path));
         return GFILE_REWIND_ERROR;
      }
   }
   else
   {
      THROW_ERROR_MSG ("Rewinding file \"%s\" failed: Unknown file type",
                       str_get (stream->path));
      return GFILE_UNKNOWN_TYPE;
   }

   return 0;
}

/** @brief Store a character in a buffer and enlarge buffer if needed.
 *
 * Function assumes that it is used iteratively. Therefore only reallocs
 * @c buf, if @c pos is equal to @c size.\n
 * Returns @c GFILE_MEM_ERROR on reallocation problems.
 *
 * @param[in] c    Character to store.
 * @param[out] buf Storage.
 * @param[in/out]  Storage size.
 * @param[in] pos  Position to store @c c at in @c buf.
 */
static __inline__ int
s_gfile_store_char (const int c,
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
   (*buf)[pos] = (char) c;

   return 0;
}

/** @brief Read a line from file up to a certain delimiter.
 *
 * Basic function for reading lines from text files. A line is terminated by a
 * set of delimiters to be provided. If the delimiter list @c delim is empty,
 * the whole file will be read into @c buf. While reading, characters can be
 * translated according to a given table @c tr. Thereby the character to be
 * translated has to be stored in column @c GFILE_TR_FROM and its complement in
 * @c GFILE_TR_TO. Alltogether the table has @c GFILE_TR_N columns. The data
 * read is stored in @c buf. The size of @c buf is determined by @c size. If
 * the line to be read is longer than @c size, @c buf is reallocated to fit the
 * line and @c size is updated. If the line read is shorter than @c size @c buf
 * is not changed in size. Each line read is terminated by a null character
 * while any delimiter character is chopped of. Since a text file can contain
 * null characters you should not iterate over a line until '\0' occurs in
 * @c buf.\n
 * @c error is used to signal read failures. For successful reading @c error is
 * 0. Has only to be checked if 0 is returned, otherwise @c error is 0.\n
 * Returns 0 if the end of file is reached or on error. Otherwise the length of
 * the line read. If 0 is returned, @c ptr is undefined but still allocated.\n
 * Since the translation and delimiter tables are searched linearily, they
 * should not be to large.
 *
 * @param[out] error     Container for error values.
 * @param[out] buf       Storage for the line read.
 * @param[in/out] size   Size of @c buf.
 * @param[in] tr         Translation table. May be @c NULL if @c tr_size = 0.
 * @param[in] tr_size    No. of entries in @c tr.
 * @param[in] delim      Delimiters list. May be @c NULL if @c delim_size = 0.
 * @param[in] delim_size No. of delimiters.
 * @param[in] stream     File to be read.
 */
static unsigned long
gfile_getdelim_tr (int* error,
                   char** buf,
                   size_t* size,
                   char tr[][GFILE_TR_N],
                   const unsigned long tr_size,
                   const char* delim,
                   const unsigned long delim_size,
                   const GFile* stream)
{
   int c;                       /* current character */
   unsigned long length = 0;    /* length of current line */
   unsigned long i;

   assert (error);
   assert (buf);
   assert ((*size == 0) || (*buf));
   assert (size);
   assert (stream);
   assert (stream->fileptr.uc);
   assert (stream->path);
   assert ((tr_size == 0) || ((tr != NULL) && (*tr != NULL)));
   assert ((delim_size == 0) || (delim != NULL));

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
         for (i = 0; i < delim_size; i++)
         {
            if ((*buf)[length] == delim[i])
            {
               /* read until newline */
               while (((char) c != CRB_LF) && (c != EOF))
               {
                  c = fgetc (stream->fileptr.uc);
               }

               (*buf)[length] = '\0';
               return length + 1;
            }
         }

         /* translate */
         for (i = 0; i < tr_size; i++)
         {
            if ((*buf)[length] == tr[i][GFILE_TR_FROM])
            {
               (*buf)[length] = tr[i][GFILE_TR_TO];
            }
         }

         length++;
      }
      
      /* handle errors of fgetc */
      if (ferror(stream->fileptr.uc))
      {
         THROW_ERROR_MSG ("Reading from file \"%s\" failed:",
                          str_get (stream->path));
         *error = GFILE_READ_ERROR;
         clearerr (stream->fileptr.uc);
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

/** @brief Read a line from file.
 *
 * Basic function for reading lines from text files. The data read is stored in
 * @c buf. The size of @c buf is determined by @c size. If the line to be read
 * is longer than @c size, @c buf is reallocated to fit the line and @c size is
 * updated. If the line read is shorter than @c size @c buf is not changed in
 * size. Each line read is terminated by a null character while any newline
 * character is chopped of. Since a text file can contain null characters you
 * should not iterate over a line until '\0' occurs in @c buf.\n
 * @c error is used to signal read failures. For successful reading @c error is
 * 0. Has only to be checked if 0 is returned, otherwise @c error is 0.\n
 * Returns 0 if the end of file is reached or on error. Otherwise the length of
 * the line read. If 0 is returned, @c ptr is undefined but still allocated.
 *
 * @param[out] error   Container for error values.
 * @param[out] buf     Storage for the line read.
 * @param[in/out] size Size of @c buf.
 * @param[in] stream   File to be read.
 */
unsigned long
gfile_getline_verbatim (int* error, char** buf, size_t* size, GFile* stream)
{
   char delim[] = {CRB_LF};

   return gfile_getdelim_tr (error, buf, size,
                             NULL, 0,
                             delim, sizeof (delim) / sizeof (*delim),
                             stream);
}

/** @brief Read a line from file and translate tabulators into whitespaces.
 *
 * Basically does the same as @c gfile_getline_verbatim() BUT translates
 * tabulators into whitespaces. For more details please refer to the
 * @c gfile_getline_verbatim().
 *
 * @param[out] error   Container for error values.
 * @param[out] buf     Storage for the line read.
 * @param[in/out] size Size of @c buf.
 * @param[in] stream   File to be read.
 */
unsigned long
gfile_getline_tab (int* error, char** buf, size_t* size, GFile* stream)
{
   char delim[] = {CRB_LF};
   char tr[][GFILE_TR_N]  = {{CRB_TAB,' '}};

   return gfile_getdelim_tr (error, buf, size,
                             tr, sizeof (tr) / sizeof (*tr),
                             delim, sizeof (delim) / sizeof (*delim),
                             stream);   
}

/** @brief Read a line from file.
 *
 * Basically does the same as @c gfile_getline_verbatim() but a line ends on
 * newline or the shell command symbol '#'. Additionally translates
 * tabulators into whitespaces. For more details please refer to the
 * @c gfile_getline_verbatim().
 *
 * @param[out] error   Container for error values.
 * @param[out] buf     Storage for the line read.
 * @param[in/out] size Size of @c buf.
 * @param[in] stream   File to be read.
 */
unsigned long
gfile_getline (int* error, char** buf, size_t* size, const GFile* stream)
{
   char delim[] = {CRB_LF, CRB_COM};
   char tr[][GFILE_TR_N]  = {{CRB_TAB,' '}};

   return gfile_getdelim_tr (error, buf, size,
                             tr, sizeof (tr) / sizeof (*tr),
                             delim, sizeof (delim) / sizeof (*delim),
                             stream);   
}

/* returns a value < 0 on error. No specific error val, can't be handled
   anyway */
int
#ifdef __GNUC__
__attribute__ ((format (printf, 2, 3)))
#endif /* __GNUC__ */
gfile_printf (GFile* gfile, const char *format, ...)
{
   int ret_val = -1;
   va_list ap;

   assert (gfile);
   assert (gfile->fileptr.uc);
   assert (gfile->path);

   errno = 0;

   switch (gfile->type)
   {
      case GFILE_UNCOMPRESSED:
         /* fetch errors */
         va_start (ap, format);
         ret_val = mvfprintf (gfile->fileptr.uc, format, ap);
         va_end (ap);

         if (ret_val < 0)
         {
            /* check for EAGAIN, EINTR -> try once again */
            /* exit with error */
            if (ferror (gfile->fileptr.uc))
            {
               THROW_WARN_MSG ("Problem while writing file \"%s\":",
                               str_get (gfile->path));
               clearerr (gfile->fileptr.uc);

               if ((errno == EAGAIN) || (errno == EINTR))
               {
                  va_start (ap, format);
                  ret_val = mvfprintf (gfile->fileptr.uc, format, ap);
                  va_end (ap);
               }
            }

            if (ret_val < 0)
            {
               THROW_ERROR_MSG ("Writing file \"%s\" failed:",
                                str_get (gfile->path));
            }
         }
         break;            
      default:
         THROW_ERROR_MSG ("Writing file \"%s\" failed: Unknown file type",
                          str_get (gfile->path));
         return -1;
   }

   return ret_val;
}
