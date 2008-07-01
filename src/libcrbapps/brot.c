/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file brot.c
 *
 *  @brief BROT, Basic RNAsequence Optimisation Tool
 *
 *  Module: brot
 *
 *  Library: libcrbapps
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-02-28
 *
 *
 *  Revision History:
 *         - 2008Feb28 bienert: created
 *
 */


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <libcrbbasic/crbbasic.h>
#include <libcrbbrot/crbbrot.h>
#include <libcrbrna/crbrna.h>
#include "brot_cmdline.h"
#include "brot.h"


static unsigned long*
parse_input_structure (const char* structure)
{
   size_t struct_len;
   unsigned long open_pair = 0;
   unsigned long i;
   unsigned long* pairlist = NULL;
   unsigned long* p_stack;
   int error = 0;

   assert (structure);

   struct_len = strlen (structure);

   /* allocate list */
   pairlist = XCALLOC (struct_len, sizeof (*pairlist));
   if (pairlist == NULL)
   {
      return pairlist;
   }

   p_stack = XCALLOC (struct_len, sizeof (*pairlist));

   i = 0;
   while ((!error) && (i < struct_len))
   {
      if (  (structure[i] != '(')
          &&(structure[i] != ')')
          &&(structure[i] != '.'))
      {
         THROW_ERROR_MSG ("Non-valid symbol found in structure string. "
                          "Allowed characters: \'(\', \')\', \'.\'. Found "
                          "\'%c\' at position %lu.", structure[i],
                         i);
         error = 1;
      }

      switch (structure[i])
      {
         case '(' :
            p_stack[open_pair] = i;
            open_pair++;
            break;     

         case ')' :
            if (open_pair == 0)
            {
               THROW_ERROR_MSG ("Mismatched nucleotide (closing base pair "
                                "partner) found at position %lu in structure "
                                "string",
                                i);
               error = 1;               
            }
            else
            {
               open_pair--;
               pairlist[i] = p_stack[open_pair] + 1;
               pairlist[p_stack[open_pair]] = i + 1;
            }
            break;

         case '.' : ;
            break;     

         default  : 
            THROW_ERROR_MSG ("Non-valid symbol found in structure string. "
                             "Allowed characters: \'(\', \')\', \'.\'. Found "
                             "\'%c\' at position %lu.",
                             structure[i], i);
            error = 1;
            break;
      }

      i++;
   }

   if (open_pair != 0)
   {
      THROW_ERROR_MSG ("Mismatched nucleotide (opening base pair "
                       "partner) found in structure string");
      error = 1;
   }
   
   XFREE (p_stack);

   if (error)
   {
      XFREE (pairlist);
      return NULL;
   }

   return pairlist;
}


static PresetArray*
parse_base_presettings (const struct brot_args_info* args_info)
{
   unsigned long i;
   unsigned long position = 0;
   unsigned long* already_set;
   PresetArray* ps;
   char* endptr;
   char base = CHAR_UNDEF;
   size_t struct_len;
   size_t len;
   int error = 0;

   struct_len = strlen (args_info->inputs[1]);

   /* allocate presettings */
   ps = PRESETARRAY_NEW_SIZE (args_info->fixed_nuc_given);
   if (ps == NULL)
   {
      return ps;
   }

   already_set = XMALLOC (sizeof (unsigned long) * struct_len);
   already_set = memset (already_set, 0, sizeof (unsigned long) * struct_len);

   /* check presettings */
   /* for (i = 0; i < args_info->fixed_nuc_given; i++) */
   i = 0;
   while ((i < args_info->fixed_nuc_given) && (!error))
   {
      /* check alphabet */
      len = strlen (args_info->fixed_nuc_arg[i]);
      if (  (len > 2)
          &&(args_info->fixed_nuc_arg[i][1] == ':'))
      {
         base = transform_base_2_number (args_info->fixed_nuc_arg[i][0]);
         if (base == CHAR_UNDEF)
         {
            error = 1;
         }
      }
      else
      {
         THROW_ERROR_MSG ("Found fixed base of wrong format: \'%s\'. Try "
                          "`%s --help` for more information.", 
                          args_info->fixed_nuc_arg[i], get_progname());
         error = 1;
      }

      if (!error)
      {
         position = strtoul (args_info->fixed_nuc_arg[i] + 2, &endptr, 10);
         if (*endptr != '\0')
         {
            THROW_ERROR_MSG ("Fixed bases require a positive integer as "
                             "position, found: %s. Try `%s --help` for more "
                             "information.", 
                             args_info->fixed_nuc_arg[i], get_progname());
            error = 1;
         }
      }

      if ((position >= struct_len) && (!error))
      {
         THROW_ERROR_MSG ("Preset position (%lu) larger than or equal to "
                          "structure length (%zu)", 
                          position, struct_len);
         error = 1;        
      }

      /* check whether position is already set */
      if ((already_set[position] > 0) && (!error))
      {
         THROW_ERROR_MSG ("Presetting conflict for position %lu (%c and %c "
                          "given), only one base allowed",
                          position,
                          presetarray_get_ith_base (already_set[position]-1,ps),
                          transform_number_2_base (base));
         error = 1;
      }

      if (!error)
      {
         /* transform base characters to numbers */
         error = presetarray_add (base, position, ps);
         already_set[position] = presetarray_get_length (ps);
      }

      i++;
   }

   XFREE (already_set);

   if (error)
   {
      presetarray_delete (ps);
      return NULL;
   }

   return ps;
}

static PresetArray* 
brot_cmdline_parser_postprocess (const struct brot_args_info* args_info)
{
   /* check input structure */

   if (args_info->inputs_num == 1)
   {
      THROW_ERROR_MSG ("RNA structure required as argument, try "
                       "`%s --help` for more information.", get_progname());
      return NULL;
   }
   if (args_info->inputs_num != 2)
   {
      THROW_ERROR_MSG ("Only one RNA structure allowed as argument, try "
                       "`%s --help` for more information.", get_progname());
      return NULL;
   }

   /* check steps */
   if (args_info->steps_given)
   {
      if (args_info->steps_arg < 0)
      {
         THROW_ERROR_MSG ("Option \"--steps\" requires positive integer as"
                          "argument, found: %ld", args_info->steps_arg);
         return NULL; 
      }
   }
 
   /* check temperature */
   if (args_info->temp_given)
   {
      if (args_info->temp_arg < 0)
      {
         THROW_ERROR_MSG ("Option \"--temp\" requires positive floating point "
                          "value as argument, found: %2.2f",
                          args_info->temp_arg);
         return NULL; 
      }
   }
  
   return parse_base_presettings (args_info);

}

/*
- Input:
  - parameters -> "grow@implementation"
- Output: Sequence + x
- flow:
  - read input & verify
  - set matrix-functions
  - simulate
    - init sequence matrix
    - simulate
    - read out sequence
    - output
  - output
*/

int
brot_main(const char *cmdline)
{
   struct brot_args_info brot_args;
   PresetArray* presets = NULL;
   SeqMatrix* sm = NULL;
   int retval = 0;
   unsigned long i;
   unsigned long* pairlist = NULL;
   float** score_matrix;

   /* command line parsing */
   brot_cmdline_parser_init (&brot_args);

   retval = brot_cmdline_parser_string (cmdline, &brot_args, get_progname());

   if (retval == 0)
   {
      retval = brot_cmdline_parser_required (&brot_args, get_progname());
   }

   /* postprocess arguments */
   if (retval == 0)
   {
      /* get fixed sites */
      presets = brot_cmdline_parser_postprocess (&brot_args);
      if (presets == NULL)
      {
         retval = 1;
      }
   }

   if (retval == 0)
   {
      /* Analyse & transform structure */
      pairlist = parse_input_structure (brot_args.inputs[1]);
      if (pairlist == NULL)
      {
         retval = 1;
      }
   }

   /* do work */
   if (retval == 0)
   {
      for (i = 0; i < presetarray_get_length (presets); i++)
      {
         mprintf ("Pos: %.2lu : Base: %c (%d)\n",
                                         presetarray_get_ith_pos (i, presets),
                transform_number_2_base (presetarray_get_ith_base (i, presets)),
                                         presetarray_get_ith_base (i, presets));
      }
   }

   /* init matrix */
   if (retval == 0)
   {
      sm = SEQMATRIX_NEW;
      retval = seqmatrix_init (pairlist,
                               strlen (brot_args.inputs[1]),
                               presets,
                               sm);
   }

   score_matrix = create_scoring_matrix ();
   if (score_matrix == 0)
   {
      retval = 1;
   }

   /* simulate */
   if (retval == 0)
   {
      seqmatrix_print_2_stderr (3, sm);
      retval = sequence_matrix_simulate_scmf (brot_args.steps_arg,
                                              brot_args.temp_arg,
                                              sm,
                                              score_matrix);
   }

   /* output */
   if (retval == 0)
   {
      mfprintf (stderr, "COLLATING\n");
      /* retval = seqmatrix_collate_is (0.99,
                                        brot_args.steps_arg / 2,
                                        brot_args.temp_arg,
                                        score_matrix, sm); */
      retval = seqmatrix_collate_mv (sm);
   }

   if (retval == 0)
   {
      seqmatrix_printf_sequence (sm);
      mprintf ("\n");
   }

   /* finalise */
   brot_cmdline_parser_free (&brot_args);
   presetarray_delete (presets);
   seqmatrix_delete (sm);
   XFREE (pairlist);
   XFREE_2D ((void**)score_matrix);

   if (retval == 0)
   {
      return EXIT_SUCCESS;
   }
   else
   {
      return EXIT_FAILURE;
   }
}
