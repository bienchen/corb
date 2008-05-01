/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbrot/seqmatrix.c
 *
 *  @brief Sequence matrix for SCMF
 *
 *  Module: seqmatrix
 *
 *  Library: libcrbbrot
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-03-05
 *
 *
 *  Revision History:
 *         - 2008Mar05 bienert: created
 *
 */


#include <config.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <libcrbbasic/crbbasic.h>
#include "preset.h"
#include "seqmatrix.h"

/* use until something better is found */
#define SM_ROWS 4               /* HP: 2 letters, RNA: 4 */

struct SeqMatrix {
      char* fixed_sites;
      float** matrix;
      unsigned long* pairlist;
      size_t rows;
      size_t cols;   
};


/**********************   Constructors and destructors   **********************/

/** @brief Create a new sequence matrix.
 *
 * The constructor for @c SeqMatrix objects. If compiled with memory checking
 * enabled, @c file and @c line should point to the position where the function
 * was called. Both parameters are automatically set by using the macro
 * @c SEQMATRIX_NEW.\n
 * Returns @c NULL on error.
 *
 * @param[in] file fill with name of calling file.
 * @param[in] line fill with calling line.
 */
SeqMatrix*
seqmatrix_new (const char* file, const int line)
{
   /* allocate 1 object */
   SeqMatrix* sm = XOBJ_MALLOC(sizeof (SeqMatrix), file, line);

   if (sm != NULL)
   {
      sm->fixed_sites = NULL;
      sm->pairlist    = NULL;
      sm->matrix      = NULL;
   }

   return sm;
}

/** @brief Delete a sequence matrix.
 *
 * The destructor for @c SeqMatrix objects.
 *
 * @param[in] sm object to be freed.
 */
void
seqmatrix_delete (SeqMatrix* sm)
{
   if (sm != NULL)
   {
      XFREE    (sm->fixed_sites);
      XFREE_2D ((void**)sm->matrix);
      XFREE    (sm->pairlist);
      XFREE    (sm);
   }
}


/*********************************   Access   *********************************/
/** @brief Check if a row is set as 'fixed'.
 *
 * Check if a certain row of a sequence matrix is preempted with a fixed site.\n
 * Returns 'true', if the site is set, 'false' otherwise.
 *
 * @param[in] row The row to check.
 * @param[in] sm The sequence matrix
 */
__inline__ bool
seqmatrix_is_col_fixed (unsigned long col, SeqMatrix* sm)
{
   assert (sm);
   assert (sm->matrix);
   assert (sm->fixed_sites);
   /* assert (); */

   if (sm->fixed_sites[(col / CHAR_BIT)] & (1 << (col % CHAR_BIT)))
   {
      return true;
   }

   return false;
}


/********************************   Altering   ********************************/

/** @brief Initialise a sequence matrix.
 *
 * On initialisation, all probabilities of the matrix are set. Additionally, a
 * set of sites can be fixed to exclude them from the simulation.\n
 * Returns 0 on success, @c ERR_SM_ALLOC indicating memory problems.
 *
 * @param[in] size Size of the structure.
 * @param[in] presettings List of fixed sites.
 * @param[in] sm Sequence matrix to initialise.
 */
int
seqmatrix_init (const unsigned long* pairs,
                const unsigned long size,
                const PresetArray* presettings,
                SeqMatrix* sm)
{
   unsigned long i, j;
   unsigned long pslen;
   unsigned long row, col;
   float init_row[size];

   assert (sm);
   assert (sm->fixed_sites == NULL);
   assert (sm->pairlist    == NULL);
   assert (sm->matrix      == NULL);

   /* allocate matrix */
   sm->rows = SM_ROWS;
   sm->cols = size;
   sm->matrix = (float**) XMALLOC_2D (sm->rows,sm->cols,sizeof (*(sm->matrix)));
   if (sm->matrix == NULL)
   {
      return ERR_SM_ALLOC;
   }
   
   /* init sm */
   init_row[0] = 1.0f / SM_ROWS;
   for (i = 1; i < size; i++)
   {
      init_row[i] = init_row[0];
   }
   
   for (i = 0; i < SM_ROWS; i++)
   {
      memcpy (sm->matrix[i], init_row, sizeof (*(sm->matrix)) * size);
   }

   /* copy pairlist */
   sm->pairlist = XCALLOC (size, sizeof (*(sm->pairlist)));

   if (pairs != NULL)
   {
      for (i = 0; i < size; i++)
      {
         sm->pairlist[i] = pairs[i];
      }
   }

   /* run over list of presettings and set sites */
   /* allocate memory & init fixed sites */
   sm->fixed_sites = XCALLOC ((size / (sizeof(*(sm->fixed_sites))*CHAR_BIT))+1,
                              sizeof(*(sm->fixed_sites)));

   if (sm->fixed_sites == NULL)
   {
      return ERR_SM_ALLOC;
   }

   if (presettings != NULL)
   {
      /* use, if rows are adjacent in memory
        init_row[0] = 0.0f;
        for (i = 1; i < SM_ROWS; i++)
        {
        init_row[i] = init_row[0];
        }
      */

      pslen = presetarray_get_length (presettings);
      for (i = 0; i < pslen; i++)
      {
         row = presetarray_get_ith_base (i, presettings);
         col = presetarray_get_ith_pos (i, presettings);

         /* set site as fixed */
         sm->fixed_sites[(col / CHAR_BIT)] |= 
            1 << (col % CHAR_BIT);

         /* set probability in matrix to 1 */
         /* use, if rows are the connected memory areas
           memcpy (sm->matrix[row], init_row, sizeof (*(sm->matrix)) * SM_ROWS);
         */
         for (j = 0; j < SM_ROWS; j++)
            sm->matrix[j][col] = 0.0f; 
         sm->matrix[row][col] = 1.0f;
      }
   }

   return 0;
}


/** @brief Update a column of a sequence matrix
 *
 * Update a column of a sequence matrix in a SCMF simulation.
 * Returns...
 *
 * @params[in] sm Sequence matrix.
 */
static __inline__ int
sequence_matrix_update_col_scmf (unsigned long row, SeqMatrix* sm)
{
   unsigned long j;

   mfprintf (stderr, "R: %lu ", row);
   for (j = 0; j < sm->cols; j++)
   {
      /* check if fixed */
      if (!seqmatrix_is_col_fixed (j, sm))
      {
         mfprintf (stderr, "C:%lu ", j);
      }

      /* update cell */
      
   }

   return 0;
}


/** @brief Update the rows of a sequence matrix
 *
 * Update rows of a sequence matrix in a SCMF simulation.
 * Returns...
 *
 * @params[in] sm Sequence matrix.
 */
static __inline__ int
sequence_matrix_update_row_scmf (SeqMatrix* sm)
{
   unsigned long i;
   int error = 0;

   /* for each row */
   i = 0;
   while ((!error) && (i < sm->rows))
   {
      error = sequence_matrix_update_col_scmf (i, sm);        

      /* update columns */
      i++;
   }

   return error;
}


/** @brief Perform a SCMF simulation on a sequence matrix.
 *
 * Calculate the mean force field for a sequence matrix and update cells. This
 * is done either until the system converges or for a specified number of
 * steps. For results of the simulation, the matrix itself has to be
 * interpreted.\n
 * Returns ...
 *
 * @params[in] steps No. of max. simulation steps.
 * @params[in] sm The sequence matrix.
 */
int
sequence_matrix_simulate_scmf (const unsigned long steps,
                               SeqMatrix* sm)
{
   unsigned long t;
   int error= 0;

   assert (sm);

   /* perform for a certain number of steps */
   t = 0;
   while ((!error) && (t < steps))
   {
      mfprintf (stderr, "Step %lu: ", t);

      /* update rows */
      error = sequence_matrix_update_row_scmf (sm);
      t++;
      mfprintf (stderr, "\n");
   }

   return error;
}

/*********************************   Output   *********************************/

/** @brief Print a sequence matrix to a stream.
 *
 * @params[in] stream Output stream to write to. FILE *stream
 * @params[in] p  Precision of the numbers.
 * @params[in] sm The sequence matrix.
 */
void
seqmatrix_fprintf (FILE* stream, const int p, const SeqMatrix* sm)
{
   unsigned long line_width = 0;
   unsigned long i, j;
   int rprec = 0;               /* precision for row number */
   int cprec = 0;               /* precision for cells */
   char* string;
   char* string_start;

   assert (sm);
   assert (p <= FLT_DIG); /* check that p stays in the precision of float */
   
   /* to circumvent problems with multiple threads, interfering prints...
      we first create the whole matrix in the memory */
   
   /* find the widest line */
   for (i = 0; i < sm->rows; i++)
   {
      for (j = 0; j < sm->cols; j++)
      {
         if (sm->matrix[i][j] > 1.0f)
         {
            rprec = floor (log10 (sm->matrix[i][j]) + 1.0f);
         }
         else
         {
            rprec = 1;
         }

         if (rprec > cprec) 
         {
            cprec = rprec;
         }      
      }
   }
   cprec += p + 1;

   if (sm->rows > 1.0f)
   {
      rprec = floor (log10 (sm->rows) + 1.0f);
   }
   else
   {
      rprec = 1;
   }

   /* add up components of a line */
   line_width += rprec;
   line_width += 4;              /* + ':  |' row start */
   line_width += ((3 + cprec) * sm->cols); /* + '.' + '  |' */
   line_width += 1;              /* + '\n' */

   /* alloc memory for the string */
   string = XMALLOC (sizeof (char) * (((line_width) * sm->rows) + 1));
   string_start = string;

   /* write matrix */
   for (i = 0; i < sm->rows; i++)
   {
      msprintf (string, "%*lu:  | ", rprec, i);
      string += rprec + 4;

      for (j = 0; j < sm->cols; j++)
      {
         msprintf (string, " %*.*f |", cprec, p, sm->matrix[i][j]);
         string += 3 + cprec;
      }

      string[0] = '\n';
      string++;
   }   
   string[0] = '\0';

   string = string_start;
   mfprintf (stream, "%s\n", string);

   XFREE (string);
}


/** @brief Print a sequence matrix to stdout.
 *
 * @params[in] p  Precision of the numbers.
 * @params[in] sm The sequence matrix.
 */
void
seqmatrix_print_2_stdout (const int p, const SeqMatrix* sm)
{
   seqmatrix_fprintf (stdout, p, sm);
}


/** @brief Print a sequence matrix to stderr.
 *
 * @params[in] p  Precision of the numbers.
 * @params[in] sm The sequence matrix.
 */
void
seqmatrix_print_2_stderr (const int p, const SeqMatrix* sm)
{
   seqmatrix_fprintf (stderr, p, sm);
}
