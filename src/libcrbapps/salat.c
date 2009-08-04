/*
 * Copyright (C) 2009 Stefan Bienert
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbapps/salat.c
 *
 *  @brief saLaT, simply annotate Loop-assembly Topologies
 *
 *  Module: salat
 *
 *  Library: libcrbapps
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2009-07-21
 *
 *
 *  Revision History:
 *         - 2009Jul21 bienert: created
 *
 */

#include <config.h>
#include <stdlib.h>
#include <libcrbbasic/crbbasic.h>
#include <libcrbrna/crbrna.h>
#include "salat_cmdline.h"
#include "salat.h"


static int
salat_cmdline_parser_postprocess (const struct salat_args_info* args_info)
{
   /* check for given input structure */
   if (args_info->inputs_num < 2)
   {
      THROW_ERROR_MSG ("Exactly one RNA secondary structure required as "
                       "argument, try %s --help` for more information.",
                       get_progname());
      return 1;      
   }

   if (args_info->inputs_num > 2)
   {
      THROW_ERROR_MSG ("Only one RNA secondary structure allowed as "
                       "argument, try %s --help` for more information.",
                       get_progname());
      return 1;
   }

   return 0;
}

int
salat_main(const char *cmdline)
{
   int retval = 0;
   struct salat_args_info salat_args;
   Rna* rna = NULL;
   SecStruct* structure;
   unsigned int i;
   unsigned long j, k, p5, p3, pq5, pq3, loop_size, loop_size2;
   unsigned long size;
   unsigned long n, m;
   bool found, printed;

   /* command line parsing */
   salat_cmdline_parser_init (&salat_args);

   retval = salat_cmdline_parser_string (cmdline, &salat_args, get_progname());

   if (retval == 0)
   {
      retval = salat_cmdline_parser_required (&salat_args, get_progname());
   }

   /* postprocess arguments */
   if (retval == 0)
   {
      retval = salat_cmdline_parser_postprocess (&salat_args);
   }

   /* store structure */
   rna = RNA_NEW;
   if (rna == NULL)
   {
      retval = 1;
   }

   if (retval == 0)
   {
      retval = RNA_INIT_PAIRLIST_VIENNA (salat_args.inputs[1],
                                         strlen (salat_args.inputs[1]),
                                         rna);
   }

   /* explore secondary structure */
   if (retval == 0)
   {
      retval = RNA_SECSTRUCT_INIT (rna);
   }

   if (retval == 0)
   {
      structure = rna_get_secstruct (rna);
      size = rna_get_size (rna);
      /* without positions, just output what we have */
      if (salat_args.position_given == 0)
      {
         mfprintf (stdout, "Stacked base pairs:\n");
         secstruct_fprintf_stacks (stdout, structure);
         mfprintf (stdout, "\nHairpin loops:\n");
         secstruct_fprintf_hairpins (stdout, structure);
         mfprintf (stdout, "\nBulge loops:\n");
         secstruct_fprintf_bulges (stdout, structure);
         mfprintf (stdout, "\nInternal loops:\n");
         secstruct_fprintf_internals (stdout, structure);
         mfprintf (stdout, "\nExternal loop:\n");
         secstruct_fprintf_external (stdout, structure);
         mfprintf (stdout, "\nMulti loops:\n");
         secstruct_fprintf_multiloops (stdout, structure);
         mfprintf (stdout, "\n");
      }
      else
      {
         /* search for a certian position */
         for (i = 0; i < salat_args.position_given; i++)
         {
            printed = false;

            if ((unsigned) salat_args.position_arg[i] > size)
            {
               THROW_WARN_MSG ("Position \"%li\" not in scope of structure "
                               "\"%s\", skipping.",
                               salat_args.position_arg[i],
                               salat_args.inputs[1]);
               printed = true;
            }
            else
            {
               mfprintf (stdout, "Position %li:\n", salat_args.position_arg[i]);

               /* stacks */
               n = secstruct_get_noof_stacks (structure);
               for (j = 0; j < n; j++)
               {
                  secstruct_get_i_geometry_stack (&p5, &p3, j, structure);
                  if (   (p5 == (unsigned) salat_args.position_arg[i])
                      || (p3 == (unsigned) salat_args.position_arg[i]))
                  {
                     mfprintf (stdout, "  Stacked base pair:\n  %lu: ", j);
                     secstruct_fprintf_i_stack (stdout, j, structure);
                     mfprintf (stdout, "\n");
                     printed = true;
                  }
               }

               /* hairpins */
               n = secstruct_get_noof_hairpins (structure);
               for (j = 0; j < n; j++)
               {
                  secstruct_get_geometry_hairpin(&p5, &p3, &loop_size, j,
                                                 structure);
                  if (   (p5 <= (unsigned) salat_args.position_arg[i])
                      && (p3 >= (unsigned) salat_args.position_arg[i]))
                  {
                     mfprintf (stdout, "  Hairpin loop:\n  %lu: ", j);
                     secstruct_fprintf_i_hairpin (stdout, j, structure);
                     mfprintf (stdout, "\n");
                     printed = true;
                  }
               }

               /* bulges */
               n = secstruct_get_noof_bulges (structure);
               for (j = 0; j < n; j++)
               {
                  secstruct_get_geometry_bulge (&p5, &p3, &pq5, &pq3,
                                                &loop_size,
                                                j,
                                                structure);

                  if (   ((p5  <= (unsigned) salat_args.position_arg[i])
                        &&(pq5 >= (unsigned) salat_args.position_arg[i])) 
                      || ((p3  >= (unsigned) salat_args.position_arg[i])
                        &&(pq3 <= (unsigned) salat_args.position_arg[i])))
                  {
                     mfprintf (stdout, "  Bulge loop:\n  %lu: ", j);
                     secstruct_fprintf_i_bulge (stdout, j, structure);
                     mfprintf (stdout, "\n");
                     printed = true;
                  }
               }

               /* internals */
               n = secstruct_get_noof_internals (structure);
               for (j = 0; j < n; j++)
               {
                  secstruct_get_geometry_internal (&p5, &p3, &pq5, &pq3,
                                                   &loop_size, &loop_size2,
                                                   j, structure);

                  if (    ((p5  <= (unsigned) salat_args.position_arg[i])
                        && (pq5 >= (unsigned) salat_args.position_arg[i])) 
                       || ((p3  >= (unsigned) salat_args.position_arg[i])
                        && (pq3 <= (unsigned) salat_args.position_arg[i])))
                  {
                     mfprintf (stdout, "  Internal loop:\n  %lu: ", j);
                     secstruct_fprintf_i_internal (stdout, j, structure);
                     mfprintf (stdout, "\n");
                     printed = true;
                  }
               }

               /* externals */
               found = false;
               /* first, search in stems assembling the ext.loop */
               n = secstruct_get_noof_stems_extloop (structure);
               for (j = 0; (j < n) && (found == false); j++)
               {
                  secstruct_get_i_stem_extloop (&p5, &p3, j, structure);

                  if (   (p5 == (unsigned) salat_args.position_arg[i])
                      || (p3 == (unsigned) salat_args.position_arg[i]))
                  {
                     found = true;
                  }
               }

               /* search at 5' end */
               n = secstruct_get_noof_5pdangles_extloop (structure);
               for (j = 0; (j < n) && (found == false); j++)
               {
                  secstruct_get_i_5pdangle_extloop (&p5, &p3, &pq5, j,
                                                    structure);

                  if (   (p5  == (unsigned) salat_args.position_arg[i])
                      || (p3  == (unsigned) salat_args.position_arg[i])
                      || (pq5 == (unsigned) salat_args.position_arg[i]))
                  {
                     found = true;
                  }
               }

               /* search at 3' end */
               n = secstruct_get_noof_3pdangles_extloop (structure);
               for (j = 0; (j < n) && (found == false); j++)
               {
                  secstruct_get_i_3pdangle_extloop (&p5, &p3, &pq3, j,
                                                    structure);

                  if (   (p5  == (unsigned) salat_args.position_arg[i])
                      || (p3  == (unsigned) salat_args.position_arg[i])
                      || (pq3 == (unsigned) salat_args.position_arg[i]))
                  {
                     found = true;
                  }
               }

               if (found == true)
               {
                  mfprintf (stdout, "  External loop:\n");
                  secstruct_fprintf_external (stdout, structure);
                  printed = true;
               }

               /* multis */
               n = secstruct_get_noof_multiloops (structure);
               for (j = 0; j < n; j++)
               {
                  found = false;

                  /* stems */
                  m = secstruct_get_i_noof_stems_multiloop (j, structure);
                  for (k = 0; (k < m) && (found == false); k++)
                  {
                     secstruct_get_i_stem_multiloop (&p5, &p3, k, j, structure);

                     if (   (p5  == (unsigned) salat_args.position_arg[i])
                         || (p3  == (unsigned) salat_args.position_arg[i]))
                     {
                        found = true;
                     }
                  }

                  /* 5' dangle */
                  m = secstruct_get_i_noof_5pdangles_multiloop (j, structure);
                  for (k = 0; (k < m) && (found == false); k++)
                  {
                     secstruct_get_i_5pdangle_multiloop (&p5, &p3, &pq5, k, j,
                                                         structure);

                     if (   (p5  == (unsigned) salat_args.position_arg[i])
                         || (p3  == (unsigned) salat_args.position_arg[i])
                         || (pq5 == (unsigned) salat_args.position_arg[i]))
                     {
                        found = true;
                     }
                  }

                  /* 3' dangle */
                  m = secstruct_get_i_noof_3pdangles_multiloop (j, structure);
                  for (k = 0; (k < m) && (found == false); k++)
                  {
                     secstruct_get_i_3pdangle_multiloop (&p5, &p3, &pq3, k, j,
                                                         structure);

                     if (   (p5  == (unsigned) salat_args.position_arg[i])
                         || (p3  == (unsigned) salat_args.position_arg[i])
                         || (pq3 == (unsigned) salat_args.position_arg[i]))
                     {
                        found = true;
                     }
                  }

                  if (found == true)
                  {
                     mfprintf (stdout, "  Mulit loop:\n");
                     secstruct_fprintf_i_multiloop (stdout, j, structure);
                     printed = true;
                  }
               }
            }

            if (printed == false)
            {
               mfprintf (stdout, "  unpaired\n");
            }
         }
      }
   }

   /* finalise */
   salat_cmdline_parser_free (&salat_args);
   rna_delete (rna);

   if (retval == 0)
   {
      return EXIT_SUCCESS;
   }
   else
   {
      return EXIT_FAILURE;
   }
}
