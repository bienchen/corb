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
            break;            
         default:
            THROW_ERROR_MSG ("Opening file \"%s\" failed: Unknown file type",
                             filepath);
            gfile_delete (gfile);
            gfile = NULL;
      }
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
         ret_val = EOF;
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
 * ...\n
 * Returns the number of items read. Values smaller than @c nobj indicate
 * end-of-file or error.
 *
 * @param[in] gfile GFile pointer
 */
size_t
gfile_read (int *error, void* ptr, size_t size, size_t nobj, GFile *stream)
{
   
   return 0;
}
