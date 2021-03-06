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
 *  ToDo:
 *         - tryout a little equilibration phase for each temperature step
 *           (suggested by Thomas Huber)
 */


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <libcrbbasic/crbbasic.h>
#include <libcrbbrot/crbbrot.h>
#include <libcrbrna/crbrna.h>
#include "brot_cmdline.h"
#include "brot.h"

#define GAS_CONST 8.314472
#define COLLATE_THRESH 0.99

static const char NN_2_SMALL_WARNING[] = "Nearest Neighbour model can only be used with "
                             "structures of size greater than 1, size of "
   "given structure (\"%s\"): %lu";

static int
brot_cmdline_parser_postprocess (const struct brot_args_info* args_info,
                                 const char* cmdline)
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

   if (args_info->verbose_given)
   {
      enable_verbose_messaging();
      print_verbose ("# This is %s %s out of the %s\n# %s\n"
                     "# Input structure: %s\n",
                     BROT_CMDLINE_PARSER_PACKAGE,
                     BROT_CMDLINE_PARSER_VERSION,
                     PACKAGE_STRING,
                     cmdline,
                     args_info->inputs[1]);
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
   print_verbose ("# Max. steps              (-s): %li\n",
                  args_info->steps_arg);
 
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
   print_verbose ("# Start-temp.             (-t): %.2f\n",
                  args_info->temp_arg);

   if (args_info->window_size_given)
   {
      if (args_info->window_size_arg < 0)
      {
         THROW_ERROR_MSG ("Option \"--window_size\" requires positive integer "
                          "as argument, found: %ld",
                          args_info->window_size_arg);
         return 1;         
      }

      if ((unsigned) args_info->window_size_arg > (strlen (args_info->inputs[1])/2 - 1))
      {
         THROW_ERROR_MSG ("Option \"--window_size\" must be less than or equal to half of the size of the input structure. Is: \"%ld\", allowed: \"%lu\"",
                          args_info->window_size_arg, (unsigned long) strlen (args_info->inputs[1])/2 - 1);
         return 1;
      }
   }

   print_verbose ("# Neg. design term        (-d): %.2f\n",
                  args_info->negative_design_scaling_arg);
   print_verbose ("# Het. term               (-h): %.2f\n",
                  args_info->heterogenity_term_scaling_arg);
   print_verbose ("# Het. term window size   (-w): %ld\n",
                  args_info->window_size_arg);
   print_verbose ("# Entropy dropoff thresh. (-e): %.2f\n",
                  args_info->sm_entropy_arg);
   print_verbose ("# Lambda                  (-l): %.2f\n",
                  args_info->lambda_arg);
   print_verbose ("# B_long                  (-o): %.2f\n",
                  args_info->beta_long_arg);
   print_verbose ("# B_short                 (-i): %.2f\n",
                  args_info->beta_short_arg);
   print_verbose ("# Speedup thresh.         (-u): %.2f\n",
                  args_info->speedup_threshold_arg);
   print_verbose ("# Min. cool. factor       (-j): %.2f\n",
                  args_info->min_cool_arg);

   return 0;
}

static int
adopt_site_presettings (const struct brot_args_info* args_info,
                        Alphabet* sigma,
                        Scmf_Rna_Opt_data* data,
                        SeqMatrix* sm)
{
   unsigned long i;
   unsigned long position = 0;
   unsigned long struct_len;
   char base = CHAR_UNDEF;
   char* endptr;

   assert (sm);

   struct_len = seqmatrix_get_width (sm);

   print_verbose ("# Fixed sites             (-n): ");

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

      seqmatrix_fix_col (base, position, data, sm);

      print_verbose ("%s ", args_info->fixed_nuc_arg[i]);
   }
   print_verbose ("\n");

   return 0;
}

static int
simulate_using_simplenn_scoring (struct brot_args_info* brot_args,
                                 SeqMatrix* sm,
                                 Scmf_Rna_Opt_data* data,
                                 GFile* entropy_file,
                                 GFile* simulation_file)
{
   int error = 0;
   char** bp_allowed = NULL;
   char bi, bj;
   unsigned long i, j, k;
   unsigned long alpha_size
      = alphabet_size (scmf_rna_opt_data_get_alphabet(data));
   unsigned long allowed_bp = 0;
   NN_scores* scores =
      NN_SCORES_NEW_INIT(0, scmf_rna_opt_data_get_alphabet (data));

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
                               sizeof (**(bp_allowed)));
      if (bp_allowed[0] == NULL)
      {
         error = 1;
      }
   }

   if (!error)
   {
      i = 0;
      k = 0;
      while (i < alpha_size)
      {
         bp_allowed[i] = bp_allowed[0] + (k * sizeof (**(bp_allowed)));
         
         for (j = 0; j < allowed_bp; j++)
         {
            nn_scores_get_allowed_basepair (j, &bi, &bj, scores);
            if (i == (unsigned) bi)
            {
               bp_allowed[0][k] = (char)(bj + 1);
               k++;
            }
         }
         
         k++;
         i++;
      }
   }

   /* simulate */
   if (!error)
   {
      scmf_rna_opt_data_set_scores (scores, data);
      scmf_rna_opt_data_set_bp_allowed (bp_allowed, data);

      seqmatrix_set_gas_constant (GAS_CONST, sm);

      seqmatrix_set_func_calc_cell_energy (scmf_rna_opt_calc_simplenn, sm);
      /*scmf_rna_opt_data_init_negative_design_energies (data, sm);*/
      seqmatrix_set_pre_col_iter_hook (
        scmf_rna_opt_data_init_negative_design_energies_alt, sm);

      seqmatrix_set_transform_row (scmf_rna_opt_data_transform_row_2_base, sm); /* SB: 27.11.09 moved here */
      seqmatrix_set_get_seq_string (scmf_rna_opt_data_get_seq_sm, sm);

      error = seqmatrix_simulate_scmf (brot_args->steps_arg,
                                       brot_args->temp_arg,
                                       brot_args->beta_long_arg,
                                       brot_args->beta_short_arg,
                                       brot_args->speedup_threshold_arg,
                                       brot_args->min_cool_arg,
                                       /*brot_args->scale_cool_arg,*/
                                       brot_args->lambda_arg,
                                       brot_args->sm_entropy_arg,
                                       entropy_file,
                                       simulation_file,
                                       sm,
                                       data);
   }

   /* collate */
   if (!error)
   {
 /* seqmatrix_set_transform_row (scmf_rna_opt_data_transform_row_2_base, sm); SB 27.11.09 - moved before simulation starts */
      seqmatrix_set_fixed_site_hook (scmf_rna_opt_data_update_neg_design_energy,
                                     sm);

      error = seqmatrix_collate_is (COLLATE_THRESH,
                                    brot_args->steps_arg / 2,
                                    brot_args->temp_arg, 
                                    brot_args->beta_long_arg,
                                    brot_args->beta_short_arg,
                                    brot_args->speedup_threshold_arg,
                                    brot_args->min_cool_arg,
                                    brot_args->lambda_arg,
                                    /* brot_args->sm_entropy_arg, */
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
simulate_using_nn_scoring (struct brot_args_info* brot_args,
                           SeqMatrix* sm,
                           Scmf_Rna_Opt_data* data,
                           GFile* entropy_file,
                           GFile* simulation_file)
{
   int error = 0;
   char** bp_allowed = NULL;
   char bi, bj;
   long int seed;
   unsigned long i, j, k;
   unsigned long alpha_size
      = alphabet_size (scmf_rna_opt_data_get_alphabet(data));
   unsigned long allowed_bp = 0;
   NN_scores* scores =
      /* SB: 17.11.09 50 */
      NN_SCORES_NEW_INIT(50.0f, scmf_rna_opt_data_get_alphabet (data));

   if (scores == NULL)
   {
      error = 1;
   }

   /* prepare index of allowed base pairs */
   if (!error)
   {
      /* "randomise" scoring function */
      print_verbose ("# Random seed             (-r): ");
      if (brot_args->seed_given)
      {
         if (brot_args->seed_arg != 0)
         {
            print_verbose ("%ld\n", brot_args->seed_arg);
            nn_scores_add_thermal_noise (alpha_size,
                                         brot_args->seed_arg,
                                         scores);
         }
         else
         {
            print_verbose ("disabled\n");
         }
      }
      else
      {
         /* if no seed is given, use time */
         seed = (long int) time(NULL);
         print_verbose ("%ld\n", seed);
         nn_scores_add_thermal_noise (alpha_size,
                                      seed,
                                      scores);
      }

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
                               sizeof (**(bp_allowed)));
      if (bp_allowed[0] == NULL)
      {
         error = 1;
      }
   }

   if (!error)
   {
      i = 0;
      k = 0;
      while (i < alpha_size)
      {
         bp_allowed[i] = bp_allowed[0] + (k * sizeof (**(bp_allowed)));
         
         for (j = 0; j < allowed_bp; j++)
         {
            nn_scores_get_allowed_basepair (j, &bi, &bj, scores);
            if (i == (unsigned) bi)
            {
               bp_allowed[0][k] = (char)(bj + 1);
               k++;
            }
         }
         
         k++;
         i++;
      }
   }

   /* decompose secondary structure */
   if (!error)
   {
      error = scmf_rna_opt_data_secstruct_init (data);
   }

   /*mfprintf (stdout, "TREATMENT of fixed sites: If both sites of a pair are "
     "fixed, delete from list? Verbose info!!!\n");*/

   /* set our special function for calc. cols.: Iterate over sec.struct., not
      sequence matrix! */
   if (!error)
   {
      scmf_rna_opt_data_set_scores (scores, data);
      scmf_rna_opt_data_set_bp_allowed (bp_allowed, data);
      scmf_rna_opt_data_set_scales (brot_args->negative_design_scaling_arg,
                                    brot_args->heterogenity_term_scaling_arg,
                                    data);
      scmf_rna_opt_data_set_het_window (brot_args->window_size_arg, data);

      seqmatrix_set_func_calc_eeff_col (scmf_rna_opt_calc_col_nn, sm);
      seqmatrix_set_gas_constant (8.314472, sm);

      seqmatrix_set_transform_row (scmf_rna_opt_data_transform_row_2_base, sm); /* SB 27.11.09 moved here */
      seqmatrix_set_get_seq_string (scmf_rna_opt_data_get_seq_sm, sm);

      error = seqmatrix_simulate_scmf (brot_args->steps_arg,
                                       brot_args->temp_arg,
                                       brot_args->beta_long_arg,
                                       brot_args->beta_short_arg,
                                       brot_args->speedup_threshold_arg,
                                       brot_args->min_cool_arg,
                                       /* brot_args->scale_cool_arg,*/
                                       brot_args->lambda_arg,
                                       brot_args->sm_entropy_arg,
                                       entropy_file,
                                       simulation_file,
                                       sm,
                                       data);
   }

   /* collate */
   if (!error)
   {
      /*seqmatrix_print_2_stdout (2, sm);*/
/* seqmatrix_set_transform_row (scmf_rna_opt_data_transform_row_2_base, sm); SB 27.11.09 moved before simulation*/

      error = seqmatrix_collate_is (COLLATE_THRESH,
                                    brot_args->steps_arg / 2,
                                    brot_args->temp_arg,
                                    brot_args->beta_long_arg,
                                    brot_args->beta_short_arg,
                                    brot_args->speedup_threshold_arg,
                                    brot_args->min_cool_arg,
                                    brot_args->lambda_arg,
                                    /* brot_args->sm_entropy_arg, */
                                    sm,
                                    data);
      /*error = seqmatrix_collate_mv (sm, data);*/
   }

   /* first: iterate scmf on secstruct, not sm! */

   nn_scores_delete (scores);
   scmf_rna_opt_data_set_scores (NULL, data);
   scmf_rna_opt_data_set_bp_allowed (NULL, data);

   XFREE (bp_allowed[0]);
   XFREE (bp_allowed);

   return error;
}

static int
simulate_using_nussinov_scoring (const struct brot_args_info* brot_args,
                                 SeqMatrix* sm, Scmf_Rna_Opt_data* data,
                                 GFile* entropy_file,
                                 GFile* simulation_file)
{
   int error = 0;
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

      /* seqmatrix_set_func_calc_eeff_row (seqmatrix_calc_eeff_row_scmf, sm);*/
      seqmatrix_set_func_calc_cell_energy (scmf_rna_opt_calc_nussinov, sm);

      seqmatrix_set_transform_row (scmf_rna_opt_data_transform_row_2_base, sm); /* SB 27.11.09 moved here */
      seqmatrix_set_get_seq_string (scmf_rna_opt_data_get_seq_sm, sm);

      error = seqmatrix_simulate_scmf (brot_args->steps_arg,
                                       brot_args->temp_arg,
                                       brot_args->beta_long_arg,
                                       brot_args->beta_short_arg,
                                       brot_args->speedup_threshold_arg,
                                       brot_args->min_cool_arg,
                                       /*brot_args->scale_cool_arg,*/
                                       brot_args->lambda_arg,
                                       brot_args->sm_entropy_arg,
                                       entropy_file,
                                       simulation_file,
                                       sm,
                                       data);
   }

   if (!error)
   {
/* seqmatrix_set_transform_row (scmf_rna_opt_data_transform_row_2_base, sm); SB 27.11.09 moved before simulation starts */

      error = seqmatrix_collate_is (COLLATE_THRESH,
                                    brot_args->steps_arg / 2,
                                    brot_args->temp_arg,
                                    brot_args->beta_long_arg,
                                    brot_args->beta_short_arg,
                                    brot_args->speedup_threshold_arg,
                                    brot_args->min_cool_arg,
                                    brot_args->lambda_arg,
                                    /* brot_args->sm_entropy_arg, */
                                    sm,
                                    data);
      /* error = seqmatrix_collate_mv (sm, sigma); */
   }

   scmf_rna_opt_data_set_scores (NULL, data);
   XFREE_2D ((void**)scores);

   return error;
}

static int
brot_settings_2_file (GFile* file, const char* cmdline)
{
   if (gfile_printf (file,
                     "# This is %s %s out of the %s\n"
                     "# %s\n",
                     BROT_CMDLINE_PARSER_PACKAGE,
                     BROT_CMDLINE_PARSER_VERSION,
                     PACKAGE_STRING,
                     cmdline) < 0)
   {
      return 1;
   }

   return 0;
}

int
brot_main(const char *cmdline)
{
   struct brot_args_info brot_args;
   SeqMatrix* sm               = NULL;
   int retval                  = 0;
   Scmf_Rna_Opt_data* sim_data = NULL;
   GFile* entropy_file = NULL;
   GFile* simulation_file = NULL;

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
      retval = brot_cmdline_parser_postprocess (&brot_args, cmdline);
   }

   /* init simulation data */
   if (retval == 0)
   {
      sim_data = SCMF_RNA_OPT_DATA_NEW_INIT(brot_args.inputs[1], 
                                            strlen (brot_args.inputs[1]),
                                            RNA_ALPHABET,
                                            strlen(RNA_ALPHABET)/2,
             ((-1) * ((logf (1 / 0.000001f)) / (strlen (brot_args.inputs[1])))),
                                            brot_args.file_given);
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
                                      scmf_rna_opt_data_get_rna_size(sim_data),
                                      /* strlen (brot_args.inputs[1]), */
                                  sm);
         /*seqmatrix_print_2_stdout (6, sm);*/
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
                                       sim_data,
                                       sm);
   }

   /* open entropy file if name given */
   if (retval == 0)
   {
      print_verbose ("# Entropy file            (-p): ");
      if (brot_args.entropy_output_given)
      {
         print_verbose ("%s", brot_args.entropy_output_arg);
         entropy_file = GFILE_OPEN (brot_args.entropy_output_arg,
                                    strlen (brot_args.entropy_output_arg),
                                    GFILE_VOID, "a");
         if (entropy_file == NULL)
         {
            retval = 1;
         }
         else
         {
            retval = brot_settings_2_file (entropy_file, cmdline);
         }
      }
      /* open simulation/ matrix file if name given */
      print_verbose ("\n# Simulation/ Matrix file (-m): ");
      if (brot_args.simulation_output_given)
      {
         print_verbose ("%s", brot_args.simulation_output_arg);
         simulation_file = GFILE_OPEN (brot_args.simulation_output_arg,
                                       strlen (brot_args.simulation_output_arg),
                                       GFILE_VOID, "a");
         
         if (simulation_file == NULL)
         {
            retval = 1;
         }
         else
         {
            retval = brot_settings_2_file (simulation_file, cmdline);
            if (retval == 0)
            {
               if (gfile_printf (simulation_file, "START\n") < 0)
               {
                  retval = 1;
               }
            }
         }
      }
   }
   
   if (retval == 0)
   {
      print_verbose ("\n# Scoring scheme          (-c): ");

      if (brot_args.scoring_arg == scoring_arg_simpleNN)
      {
         print_verbose ("simpleNN\n");
         /* special to NN usage: structure has to be of size >= 2 */
         if (strlen (brot_args.inputs[1]) > 1)
         {
            retval = simulate_using_simplenn_scoring (&brot_args,
                                                      sm,
                                                      sim_data,
                                                      entropy_file,
                                                      simulation_file);
         }
         else
         {
            THROW_ERROR_MSG (NN_2_SMALL_WARNING,
                             brot_args.inputs[1], 
                             (unsigned long) strlen (brot_args.inputs[1]));
            retval = 1;
         }
      }
      else if (brot_args.scoring_arg == scoring_arg_nussinov)
      {
         print_verbose ("nussinov\n");
         retval = simulate_using_nussinov_scoring (&brot_args,
                                                   sm,
                                                   sim_data,
                                                   entropy_file,
                                                   simulation_file);
      }
      else if (brot_args.scoring_arg == scoring_arg_NN)
      {
         print_verbose ("NN\n");
         /* special to NN usage: structure has to be of size >= 2 */
         if (strlen (brot_args.inputs[1]) > 1)
         {
            retval = simulate_using_nn_scoring (&brot_args,
                                                sm,
                                                sim_data,
                                                entropy_file,
                                                simulation_file);
         }
         else
         {
            THROW_ERROR_MSG (NN_2_SMALL_WARNING,
                             brot_args.inputs[1], 
                             (unsigned long) strlen (brot_args.inputs[1]));
            retval = 1;
         }
      }
   }

   /* close files */
   if (retval == 0)
   {
      retval = gfile_close (entropy_file);
   }
   else
   {
      gfile_close (entropy_file);
   }

   if (simulation_file != NULL)
   {
      if (gfile_printf (simulation_file, "END\n") < 0)
      {
         retval = 1;
      }
   }

   if (retval == 0)
   {
      retval = gfile_close (simulation_file);
   }
   else
   {
      gfile_close (simulation_file);
   }

   if (retval == 0)
   {
      /*seqmatrix_print_2_stdout (2, sm);*/
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
