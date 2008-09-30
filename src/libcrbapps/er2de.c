/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbapps/er2de.c
 *
 *  @brief eR2De, evaluate RNA 2D energy
 *
 *  Module: er2de
 *
 *  Library: libcrbapps
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-09-24
 *
 *
 *  Revision History:
 *         - 2008Sep24 bienert: created
 *
 */

#include <config.h>
#include <stdlib.h>
#include <libcrbbasic/crbbasic.h>
#include <libcrbrna/crbrna.h>
#include "er2de_cmdline.h"
#include "er2de.h"


static int
er2de_cmdline_parser_postprocess (const struct er2de_args_info* args_info)
{
   /* check for given input sequence and structure */
   if (args_info->inputs_num == 1)
   {
      THROW_ERROR_MSG ("RNA sequence and structure required as arguments, try "
                       "`%s --help` for more information.", get_progname());
      return 1;
   }

   if (args_info->inputs_num != 3)
   {
      THROW_ERROR_MSG ("Only one RNA sequence and one structure allowed as "
                       "arguments, try %s --help` for more information.",
                       get_progname());
      return 1;
   }

   if (strlen (args_info->inputs[1]) != strlen (args_info->inputs[2]))
   {
      THROW_ERROR_MSG ("Sequence and structure have to be of equal length: "
                       "length(sequence) = %lu, length(structure) = %lu.",
                       (unsigned long) strlen (args_info->inputs[1]),
                       (unsigned long) strlen (args_info->inputs[2]));
      return 1;      
   }

   return 0;
}

int
er2de_main(const char *cmdline)
{
   Rna* rna = NULL;
   Alphabet* sigma = NULL;
   NN_scores* scores = NULL;
   int G = 0;
   struct er2de_args_info erde_args;
   unsigned long tmp = 0;
   int retval = 0;

   /* command line parsing */
   er2de_cmdline_parser_init (&erde_args);

   retval = er2de_cmdline_parser_string (cmdline, &erde_args, get_progname());

   if (retval == 0)
   {
      retval = er2de_cmdline_parser_required (&erde_args, get_progname());
   }

   /* postprocess arguments */
   if (retval == 0)
   {
      retval = er2de_cmdline_parser_postprocess (&erde_args);
   }

   /* store sequence and structure */
   rna = RNA_NEW;
   sigma = ALPHABET_NEW_SINGLE (RNA_ALPHABET, strlen (RNA_ALPHABET) / 2);
   if (   (rna == NULL)
       || (sigma == NULL))
   {
      retval = 1;
   }

   if (retval == 0)
   {      
      retval = RNA_INIT_SEQUENCE_STRUCTURE(erde_args.inputs[1],
                                           erde_args.inputs[2],
                                           strlen (erde_args.inputs[1]),
                                           sigma,
                                           rna);
   }

   if (retval == 0)
   {
      scores = NN_SCORES_NEW_INIT(sigma);
      if (scores == NULL)
      {
         retval = 1;
      }
   }

   if (retval == 0)
   {
      tmp = rna_validate_basepairs_nn_scores (scores, rna);
      if (tmp != rna_get_size (rna))
      {
         THROW_ERROR_MSG ("Base pair \'%c%c\' not covered by the Nearest "
                          "Neighbour model. Formed by positions %lu and %lu.",
                          alphabet_no_2_base (rna_get_sequence_base (tmp, rna),
                                              sigma),
                          alphabet_no_2_base (
                             rna_get_sequence_base (rna_base_pairs_with (tmp,
                                                                         rna),
                                                    rna),
                                              sigma),
                          tmp,
                          rna_base_pairs_with (tmp, rna));
         retval = 1;
      }
   }

   /* explore secondary structure */
   if (retval == 0)
   {
      retval = RNA_SECSTRUCT_INIT (rna);
   }

   /* calculate free energy */
   if (retval == 0)
   {
      G = rna_secstruct_calculate_DG (scores, rna);
      mfprintf (stderr, "G= %5.2f\n", G * 0.01);
   }

   /* finalise */
   er2de_cmdline_parser_free (&erde_args);
   rna_delete (rna);
   alphabet_delete (sigma);
   nn_scores_delete (scores);

   if (retval == 0)
   {
      return EXIT_SUCCESS;
   }
   else
   {
      return EXIT_FAILURE;
   }
}
