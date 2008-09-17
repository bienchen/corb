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

static int
brot_cmdline_parser_postprocess (const struct brot_args_info* args_info)
{
   /* check input structure */
   if (args_info->inputs_num == 1)
   {
      THROW_ERROR_MSG ("RNA structure required as argument, try "
                       "`%s --help` for more information.", get_progname());
      return 1;
   }
   if (args_info->inputs_num != 2)
   {
      THROW_ERROR_MSG ("Only one RNA structure allowed as argument, try "
                       "`%s --help` for more information.", get_progname());
      return 1;
   }

   /* check steps */
   if (args_info->steps_given)
   {
      if (args_info->steps_arg < 0)
      {
         THROW_ERROR_MSG ("Option \"--steps\" requires positive integer as"
                          "argument, found: %ld", args_info->steps_arg);
         return 1; 
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
         return 1; 
      }
   }

   return 0;
}

static int
adopt_site_presettings (const struct brot_args_info* args_info,
                        Alphabet* sigma,
                        SeqMatrix* sm)
{
   unsigned long i;
   unsigned long position = 0;
   unsigned long struct_len;
   char base = CHAR_UNDEF;
   char* endptr;

   assert (sm);

   struct_len = seqmatrix_get_width (sm);

   /* check presettings */
   for (i = 0; i < args_info->fixed_nuc_given; i++)
   {
      /* check format & alphabet */
      if (  (strlen (args_info->fixed_nuc_arg[i]) > 2)
          &&(args_info->fixed_nuc_arg[i][1] == ':'))
      {
         base = alphabet_base_2_no (args_info->fixed_nuc_arg[i][0], sigma);
         if (base == CHAR_UNDEF)
         {
            return 1;
         }
      }
      else
      {
         THROW_ERROR_MSG ("Found fixed base of wrong format: \'%s\'. Try "
                          "`%s --help` for more information.", 
                          args_info->fixed_nuc_arg[i], get_progname());
         return 1;
      }

      position = strtoul (args_info->fixed_nuc_arg[i] + 2, &endptr, 10);
      if (*endptr != '\0')
      {
         THROW_ERROR_MSG ("Fixed bases require a positive integer as "
                          "position, found: %s. Try `%s --help` for more "
                          "information.", 
                          args_info->fixed_nuc_arg[i], get_progname());
         return 1;
      }

      if (position >= struct_len)
      {
         THROW_ERROR_MSG ("Preset position (%lu) larger than or equal to "
                          "structure length (%lu) for presetting \"%s\"", 
                          position, struct_len, args_info->fixed_nuc_arg[i]);
         return 1;        
      }

      if (seqmatrix_is_col_fixed (position, sm))
      {
         THROW_ERROR_MSG ("Presetting conflict for position %lu (\"%c\"): "
                          "Already set", position,
                          alphabet_no_2_base (base, sigma));
         return 1;
      }

      seqmatrix_fix_col (base, position, sm);      
   }

   return 0;
}

static int
simulate_using_nn_scoring (struct brot_args_info* brot_args,
                           SeqMatrix* sm, Scmf_Rna_Opt_data* data)
{
   int error = 0;
   float c_rate = 0;
   char** bp_allowed = NULL;
   char bi, bj;
   unsigned long i, j, k;
   unsigned long alpha_size
      = alphabet_size (scmf_rna_opt_data_get_alphabet(data));
   unsigned long allowed_bp = 0;
   NN_scores* scores =
      NN_SCORES_NEW_INIT(scmf_rna_opt_data_get_alphabet (data));

   if (scores == NULL)
   {
      error = 1;
   }

   /* prepare index of allowed base pairs */
   if (!error)
   {
      bp_allowed = XMALLOC(alpha_size * sizeof (*bp_allowed));
      if (bp_allowed == NULL)
      {
         error = 1;
      }
   }

   if (!error)
   {

      allowed_bp = nn_scores_no_allowed_basepairs (scores);

      /* we need 1 byte for each possible pair + 1byte for the NULL byte for
         each letter in the alphabet */
      bp_allowed[0] = XCALLOC (allowed_bp + alpha_size,
                               sizeof (*(bp_allowed[0])));
      if (bp_allowed[0] == NULL)
      {
         error = 1;
      }
   }

   i = 0;
   k = 0;
   while ((!error) && (i < alpha_size))
   {
      bp_allowed[i] = bp_allowed[0] + (k * sizeof (char));

      for (j = 0; j < allowed_bp; j++)
      {
         nn_scores_get_allowed_basepair (j, &bi, &bj, scores);
         if (i == (unsigned) bi)
         {
            bp_allowed[0][k] = bj + 1;
            k++;
         }
      }

      k++;
      i++;
   }

   /* simulate */
   if (!error)
   {
      scmf_rna_opt_data_set_scores (scores, data);
      scmf_rna_opt_data_set_bp_allowed (bp_allowed, data);

      seqmatrix_set_gas_constant (8.314472, sm);

      if (brot_args->steps_arg > 0)
      {
         c_rate = logf (brot_args->temp_arg / 1.0f);
         c_rate /= (brot_args->steps_arg - 1);
         c_rate = expf ((-1) * c_rate);
      }

      seqmatrix_set_func_calc_cell_energy (scmf_rna_opt_calc_nn, sm);
      /*scmf_rna_opt_data_init_negative_design_energies (data, sm);*/
      seqmatrix_set_pre_col_iter_hook (
        scmf_rna_opt_data_init_negative_design_energies_alt, sm);
      error = seqmatrix_simulate_scmf (brot_args->steps_arg,
                                       brot_args->temp_arg,
                                       c_rate,
                                       0,
                                       0.6,
                                       0.05,
                                       sm,
                                       data);
   }

   /* collate */
   if (!error)
   {
      seqmatrix_set_transform_row (scmf_rna_opt_data_transform_row_2_base, sm);
      seqmatrix_set_fixed_site_hook (scmf_rna_opt_data_update_neg_design_energy,
                                     sm);

      /*seqmatrix_print_2_stdout (2, sm); */
      if (brot_args->steps_arg > 0)
      {
         c_rate = logf (brot_args->temp_arg / 1.0f);
         c_rate /= ((brot_args->steps_arg - 1) / 2);
         c_rate = expf ((-1) * c_rate);
      }

      error = seqmatrix_collate_is (0.99,
                                    brot_args->steps_arg / 2,
                                    brot_args->temp_arg,
                                    c_rate,
                                    0,
                                    0.6,
                                    0.05,
                                    sm,
                                    data);
      /*error = seqmatrix_collate_mv (sm, data);*/
   }

   nn_scores_delete (scores);
   scmf_rna_opt_data_set_scores (NULL, data);
   scmf_rna_opt_data_set_bp_allowed (NULL, data);
   XFREE (bp_allowed[0]);
   XFREE (bp_allowed);

   return error;
}

static int
simulate_using_nussinov_scoring (const struct brot_args_info* brot_args,
                                 SeqMatrix* sm, Scmf_Rna_Opt_data* data)
{
   int error = 0;
   float c_rate = 0;
   float** scores
      = create_scoring_matrix (scmf_rna_opt_data_get_alphabet (data)); 

   if (scores == NULL)
   {
      error = 1;
   }

   /* simulate */
   if (!error)
   {
      scmf_rna_opt_data_set_scores (scores, data);

      if (brot_args->steps_arg > 0)
      {
         c_rate = logf (brot_args->temp_arg / 0.1f);
         c_rate /= (brot_args->steps_arg - 1);
         c_rate = expf ((-1) * c_rate);
      }

      /* seqmatrix_set_func_calc_eeff_row (seqmatrix_calc_eeff_row_scmf, sm);*/
      seqmatrix_set_func_calc_cell_energy (scmf_rna_opt_calc_nussinov, sm);
      error = seqmatrix_simulate_scmf (brot_args->steps_arg,
                                       brot_args->temp_arg,
                                       c_rate,
                                       0,
                                       0.2,
                                       0.05,
                                       sm,
                                       data);
   }

   if (!error)
   {
      seqmatrix_set_transform_row (scmf_rna_opt_data_transform_row_2_base, sm);

      if (brot_args->steps_arg > 0)
      {
         c_rate = logf (brot_args->temp_arg / 0.1f);
         c_rate /= ((brot_args->steps_arg - 1) / 2);
         c_rate = expf ((-1) * c_rate);
      }

      error = seqmatrix_collate_is (0.99,
                                    brot_args->steps_arg / 2,
                                    brot_args->temp_arg,
                                    c_rate,
                                    0,
                                    0.2,
                                    0.05,
                                    sm,
                                    data);
      /* error = seqmatrix_collate_mv (sm, sigma); */
   }

   scmf_rna_opt_data_set_scores (NULL, data);
   XFREE_2D ((void**)scores);

   return error;
}

int
brot_main(const char *cmdline)
{
   struct brot_args_info brot_args;
   SeqMatrix* sm               = NULL;
   int retval                  = 0;
   Scmf_Rna_Opt_data* sim_data = NULL;

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
      retval = brot_cmdline_parser_postprocess (&brot_args);
   }

   /* init simulation data */
   if (retval == 0)
   {
      sim_data = SCMF_RNA_OPT_DATA_NEW_INIT(brot_args.inputs[1], 
                                            strlen (brot_args.inputs[1]),
                                            RNA_ALPHABET,
                                            strlen(RNA_ALPHABET)/2,
            ((-1) * ((logf (1 / 0.000001f)) / (strlen (brot_args.inputs[1])))));
      if (sim_data == NULL)
      {
         retval = 1;
      }
   }

   /* init matrix */
   if (retval == 0)
   {
      sm = SEQMATRIX_NEW;
      if (sm != NULL)
      {
         retval = SEQMATRIX_INIT (
            alphabet_size (scmf_rna_opt_data_get_alphabet (sim_data)),
                                  strlen (brot_args.inputs[1]),
                                  sm);
         /*seqmatrix_print_2_stdout (2, sm);*/
      }
      else
      {
         retval = 1;
      }
   }

   /* fix certain sites in the matrix */
   if (retval == 0)
   {
      retval = adopt_site_presettings (&brot_args,
                                      scmf_rna_opt_data_get_alphabet (sim_data),
                                       sm);
   }

   if (retval == 0)
   {
      if (brot_args.scoring_arg == scoring_arg_NN)
      {
         /* special to NN usage: structure has to be of size >= 2 */
         if (strlen (brot_args.inputs[1]) > 1)
         {
            retval = simulate_using_nn_scoring (&brot_args, sm, sim_data);
         }
         else
         {
            THROW_ERROR_MSG ("Nearest Neighbour model can only be used with "
                             "structures of size greater than 1, size of "
                             "given structure (\"%s\"): %lu",
                             brot_args.inputs[1], 
                             (unsigned long) strlen (brot_args.inputs[1]));
            retval = 1;
         }
      }
      else if (brot_args.scoring_arg == scoring_arg_nussinov)
      {
         retval = simulate_using_nussinov_scoring (&brot_args, sm, sim_data);
      }
   }
  
   if (retval == 0)
   {
      mprintf ("%s\n", scmf_rna_opt_data_get_seq(sim_data));
   }   

   /* finalise */
   brot_cmdline_parser_free (&brot_args);
   seqmatrix_delete (sm);
   scmf_rna_opt_data_delete (sim_data);

   if (retval == 0)
   {
      return EXIT_SUCCESS;
   }
   else
   {
      return EXIT_FAILURE;
   }
}
