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
#include "errormsg.h"
#include "memmgr.h"
#include "gfile.h"


struct GFile {
   GFileType type;
};


/**********************   Constructors and destructors   **********************/
/** @brief Create a new file object.
 *
 * The constructor for @c GFile objects. If compiled with enabled
 * memory checking, @c file and @c line should point to the position where the
 * function was called. Both parameters are automatically set by using the
 * macro @c GFILE_NEW.\n
 * Returns @c NULL on error.
 *
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
GFile*
gfile_new (const char* file, const int line)
{
   /* allocate 1 object */
   GFile* this = XOBJ_MALLOC(sizeof (GFile), file, line);
   
   if (this != NULL)
   {
      this->type = GFILE_VOID;
   }

   return this;
}

/** @brief Delete a file object.
 *
 * The destructor for @c GFile objects.
 *
 * @param[in] this object to be freed.
 */
void
gfile_delete (GFile* this)
{
   if (this != NULL)
   {
     XFREE(this);
   }
}

/********************************   Altering   ********************************/

/*********************************   Access   *********************************/

/*********************************    Size    *********************************/

/********************************* Searching  *********************************/

/********************************* Comparison *********************************/


