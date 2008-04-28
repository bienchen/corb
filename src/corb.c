/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file corb.c
 *
 *  @brief CoRB - Collection of RNAanalysis Binaries
 *
 *  Module: corb
 *
 *  Library: none
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-02-25
 *
 *
 *  Revision History:
 *         - 2008Feb25 bienert: created
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <libcrbbasic/crbbasic.h>
#include "crb_cmdline.h"
#include <libcrbapps/brot.h>

static int
verify_tool (const Str* tool)
{
   /* verify name */
   if (str_compare_cstr (tool, "brot"))
   {
      THROW_ERROR_MSG ("Unknown application: \"%s\", try `%s --help` for more "
                       "information.",str_get (tool), get_progname());
      return 1;      
   }

   return 0;
}

static Str*
parse_toolname (char* argv_string)
{
   Str* tool;
   size_t i;
   size_t name_len = 0;
   size_t argv_len = strlen (argv_string);

   tool = STR_NEW;

   if (tool == NULL)
   {
      return tool;
   }

   for (i = 0; i < argv_len; i++)
   {
      if (argv_string[i] != ' ')
      {
         name_len++;
      }
      else
      {
         i = argv_len;
      }
   }

   if (str_assign_csubstr (tool, argv_string, 0, name_len) == ERR_STR_ALLOC)
   {
      XFREE (tool);
      return NULL;
   }

   return tool;
}

int main(int argc,char *argv[])
{
   Str* tool = NULL;
   struct gengetopt_args_info crb_args;
   int retval = 0;

   if (set_progname ("corb"))
   {
      return EXIT_FAILURE;
   }

   crb_cmdline_parser_init (&crb_args);
   retval = crb_cmdline_parser (argc, argv, &crb_args);
   
   if (retval == 0)
   {
      retval = crb_cmdline_parser_required();
   }   

   /* postprocess parsed options */
   if (crb_args.inputs_num != 1)
   {
      THROW_ERROR_MSG ("Exactly one applicion string is needed, try "
                       "`%s --help` for more information.",
                       get_progname());
      retval = 1;
   }

   /* parse toolname from command string */
   if (retval == 0)
   {
      tool = parse_toolname (crb_args.inputs[0]);
      if (tool == NULL)
      {
         THROW_ERROR_MSG ("No applicion name found in string provided: \"%s\", "
                          "try %s --help` for more information.",
                          crb_args.inputs[0], get_progname());
         retval = 1;
      }
   }

   /* inspect tool chosen */
   if (retval == 0)
   {
      retval = verify_tool (tool);
   }

   if (retval == 0)
   {
      retval = add_name_2_progname (str_get (tool));
   }

   /* start tool */
   if (retval == 0)
   {
      retval = brot_main(crb_args.inputs[0]/*  + str_length (tool) */);
   }

   /* finalise */
   crb_cmdline_parser_free (&crb_args);
   str_delete (tool);
   free_progname();
   FREE_MEMORY_MANAGER;

   if (retval == 0)
   {
      return EXIT_SUCCESS;
   }
   else
   {
      return EXIT_FAILURE;
   }
}
