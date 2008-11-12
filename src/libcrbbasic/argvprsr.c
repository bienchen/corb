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
 *  @file libcrbbasic/argvprsr.c
 *
 *  @brief Parsing the command line (argv)
 *
 *  Module: argvprsr
 *
 *  Library: crbbasic
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-02-04
 *
 *
 *  Revision History:
 *         - 2008Feb04 bienert: created
 *
 */


#include <config.h>
#include <inc_strg.h>
#include "argvprsr.h"
#include "memmgr.h"

struct ArgvParser {
      char* mt_mail;
      char* mt_name;
      /* array for the options */
      /* dependencies? */
};


/** @brief Create a new @c ArgvParser object.
 *
 * The constructor for @c ArgvParser objects. If compiled with enabled
 * memory checking, @c file and @c line should point to the position where the
 * function was called. Both parameters are automatically set by using the
 * macro @c ARGVPARSER_NEW.\n
 * Returns @c NULL on error.
 *
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
ArgvParser*
argvparser_new (const char* file, const int line)
{
   /* allocate 1 parser object */
   ArgvParser* obj = XOBJ_MALLOC(sizeof (*obj), file, line);

   if (obj != NULL)
   {
      obj->mt_mail = NULL;
      obj->mt_name = NULL;
   }

   return obj;
}


/** @brief Delete a @c ArgvParser object.
 *
 * The destructor for @c ArgvParser objects. If compiled with enabled
 * memory checking, @c file and @c line should point to the position where the
 * function was called. Both parameters are automatically set by using the
 * macro @c ARGVPARSER_DELETE.
 *
 * @param[out] obj object to be freed.
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
void
argvparser_delete (ArgvParser* obj)
{
   XFREE (obj);
}

/** @brief Get the maintainer's name.
 *
 * Retrieve the maintainer name stored in an @c ArgvParser. If none was set up
 * to now, @c NULL is returned.
 *
 * @param[in] argvparser The object to be read.
 */
char* argvparser_get_maintainername (ArgvParser* argvparser)
{
   return argvparser->mt_name;
}

/** @brief Set the maintainer's name.
 *
 * Store the maintainer name in an @c ArgvParser. This will then be the name
 * printed at the end of a usage message summond by @c -h or @c --help in
 * @c argv. If the name was set before, it will be overwritten. You can check
 * for this state by @c argvparser_get_maintainername. The name will be stored
 * as a copy. Providing the name as null pointer will delete any presettings.\n
 * Returns 0 on success, ERR_AP_ALLOC if memory allocation fails.
 *
 * @param[in] name The name to be stored.
 * @param[in] argvparser The object to be read.
 */
int argvparser_set_maintainername (char* name, ArgvParser* argvparser)
{
   if (name == NULL)            /* just delete name */
   {
      XFREE (argvparser->mt_name);
      argvparser->mt_name = NULL;
   }
   else                         /* set a name */
   {
      /* regardless of what is the content of mt_name, we just realloc it to
         the right size. On mt_name == NULL, realloc should act like malloc() */
      argvparser->mt_name = XREALLOC (argvparser->mt_name,
                                        sizeof (*argvparser->mt_name) 
                                      * (strlen (name) + 1));
      if (argvparser->mt_name == NULL)
      {
         return ERR_AP_ALLOC;
      }

      strcpy (argvparser->mt_name, name);
   }

   return 0;
}

/** @brief Set maintainer information shown in the usage information.
 *
 * Set the email address and program author/ maintainer, shown on calling the
 * argv parser with -h, --help. The information is stored as a copy and
 * released by deleting the argv parser. If this information was set before, it
 * will be overwritten. To check whether the information already exists, use
 * @c argvparser_get_maintainername/ argvparser_get_maintainermail.\n
 * Returns 0 on success.
 *
 * @param[in]  name The author's name.
 * @param[in]  mail The author's email address.
 * @param[out] argvparser Object to be extended.
 */
/*int
argvparser_set_maintainerinfo (const char* name,
                               const char* mail,
                               const ArgvParser* argvparser)
{
   if (argvparser->mt_name != NULL) 
   {
   
   }
   else
   {
   }
  
   }*/

/* if not set, everything NULL */
/* get_maintainermail */


/* start parser functions */
/* only args */

/* args + operand */
/* args + operandS */

/* add arguments/ options functions */
