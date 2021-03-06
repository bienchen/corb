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
 *         - 2009Nov bienert: Added new hook for seqmatrix_fix_col() 
 *
 *  ToDo:
 *
 */


#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <libcrbbasic/crbbasic.h>
#include "seqmatrix.h"

/* matrix orders */
/* enum seqmatrix_matrix_no{ */
/*    F_Mtrx = 0,      /\* first matrix *\/ */
/*    S_Mtrx,          /\* second matrix *\/ */
/*    No_Of_Mtrx       /\* overall number of matrices *\/ */
/* }; */

struct SeqMatrix {
   char* fixed_sites;          /* list of fixed sites in the matrix */
   float** prob_m;             /* probability matrix */
   float** calc_m;             /* matrix for calculation of new prob. */
   size_t rows;
   size_t cols;
   float gas_constant;
   int (*calc_eeff_col) (SeqMatrix*,
                         const float,
                         void*);
   int (*calc_eeff_row) (const unsigned long,
                         SeqMatrix*,
                         const float,
                            void*);
   float (*calc_cell_energy) (const unsigned long, const unsigned long,
                              void*,
                              SeqMatrix*);
   int (*transform_row) (const unsigned long, const unsigned long, void*);
   int (*pre_col_iter_hook) (void*, SeqMatrix*);
   /* hook for fixed sites during seqmatrix_calc_eeff_col_scmf(), e.g. to
      jump over fixed sites if your energy function is evaluated linearly */
   int (*fixed_site_hook) (void*, unsigned long, SeqMatrix*);
   /* hook to act upon sites during seqmatrix_fix_col() */
   int (*fixing_site_hook) (void*, unsigned long, SeqMatrix*);
   char* (*get_seq_string) (void*);
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
   SeqMatrix* sm = XOBJ_MALLOC(sizeof (*sm), file, line);

   if (sm != NULL)
   {
      sm->fixed_sites       = NULL;
      sm->calc_eeff_col     = NULL;
      sm->calc_eeff_row     = NULL;
      sm->calc_cell_energy  = NULL;
      sm->transform_row     = NULL;
      sm->pre_col_iter_hook = NULL;
      sm->fixed_site_hook   = NULL;
      sm->fixing_site_hook  = NULL;
      sm->get_seq_string    = NULL;
      sm->prob_m            = NULL;
      sm->calc_m            = NULL;
      sm->gas_constant  = 1;
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
   /* unsigned long i; */

   if (sm != NULL)
   {
      XFREE    (sm->fixed_sites);
      XFREE_2D ((void**)sm->prob_m);
      XFREE_2D ((void**)sm->calc_m);

/*       for (i = 0; i < No_Of_Mtrx; i++) */
/*       { */
/*          XFREE_2D ((void**)sm->matrix[i]); */
/*       } */

      XFREE    (sm);
   }
}


/*********************************   Access   *********************************/

/** @brief Check if a column is set as 'fixed'.
 *
 * Check if a certain column of a sequence matrix is preempted with a fixed
 * site.\n
 * Returns 'true', if the site is set, 'false' otherwise.
 *
 * @param[in] col The column to check.
 * @param[in] sm The sequence matrix
 */
__inline__ bool
seqmatrix_is_col_fixed (const unsigned long col, const SeqMatrix* sm)
{
   assert (sm);
   assert (sm->fixed_sites);

   if (sm->fixed_sites[(col / CHAR_BIT)] & (1 << (col % CHAR_BIT)))
   {
      return true;
   }

   return false;
}

/** @brief get the width of a sequence matrix.
 *
 * Returns the number of columns of a sequence matrix.
 *
 * @param[in] sm sequence matrix.
 */
unsigned long
seqmatrix_get_width (const SeqMatrix* sm)
{
   assert (sm);

   return sm->cols;
}

/** @brief get the no. of rows of a sequence matrix.
 *
 * Returns the number of rows of a sequence matrix.
 *
 * @param[in] sm sequence matrix.
 */
unsigned long
seqmatrix_get_rows (const SeqMatrix* sm)
{
   assert (sm);

   return sm->rows;
}

/** @brief Get the probability of a certain site and state of a sequenc ematrix.
 *
 * Retruns the value of a certain cell of the sequenc ematrix.
 *
 * @params[in] row Row.
 * @params[in] col Column.
 * @params[in] sm Sequence matrix.
 */
float
seqmatrix_get_probability (const unsigned long row, const unsigned long col,
                           const SeqMatrix* sm)
{
   assert (sm);
   assert (sm->prob_m);
   assert (row < sm->rows);
   assert (col < sm->cols);

   return sm->prob_m[row][col];
}

/** @brief Get the effective energye stored in a certain site and state.
 *
 * Retruns the value of a cell of the effective energy matrix.
 *
 * @params[in] row Row.
 * @params[in] col Column.
 * @params[in] sm Sequence matrix.
 */
float
seqmatrix_get_eeff (const unsigned long row, const unsigned long col,
                    const SeqMatrix* sm)
{
   assert (sm);
   assert (sm->calc_m);
   assert (row < sm->rows);
   assert (col < sm->cols);

   return sm->calc_m[row][col];
}

/** @brief Get the gas constant.
 *
 * @params[in] sm Sequence matrix
 */
float
seqmatrix_get_gas_constant (const SeqMatrix* sm)
{
   assert (sm);

   return sm->gas_constant;
}

/********************************   Altering   ********************************/

/** @brief Sets the cells of the effective energy matrix to 0.
 *
 * This function can be used to set all cells of the effective energy matrix
 * to 0. Should usually onlly be used when writing coloumn/ row iteration of
 * the scmf simulation by yourself.
 *
 * @param[in] sm sequence matrix.
 */
void
seqmatrix_set_eeff_matrix_zero (SeqMatrix* sm)
{
   sm->calc_m = 
      (float**) matrix2d_set_zero ((void**)sm->calc_m,
                                   sm->rows, sm->cols,
                                   sizeof (**(sm->calc_m)));
}

/** @brief Switch the currently enabled matrix.
 *
 * For use when writing large parts of the simulation on your own. Switches the matrix whose values are accessible/ editable
 */
/*void
seqmatrix_switch_current_matrix (SeqMatrix* sm)
{
}*/

/** @brief Fix a certain column in a sequence matrix.
 *
 * Sets a certain cell of a sequenc ematrix to 1 and marks it to be preserved.
 * 
 * @param[in] row row to be fixed.
 * @param[in] col column to be fixed.
 * @param[in] sm sequence matrix.
 */
int
seqmatrix_fix_col (const unsigned long row, const unsigned long col,
                   void* data, SeqMatrix* sm)
{
   unsigned long i;

   assert (sm);
   assert (sm->prob_m);
   assert (sm->calc_m);
   assert (sm->fixed_sites);
   assert (col < sm->cols);
   assert (row < sm->rows);

   /* fix site */

   sm->fixed_sites[(col / CHAR_BIT)] = 
      (char) (sm->fixed_sites[(col / CHAR_BIT)] | (1 << (col % CHAR_BIT)));

   /* set everything to 0 */
   for (i = 0; i < sm->rows; i++)
   {
      sm->prob_m[i][col] = 0.0f;
      sm->calc_m[i][col] = 0.0f;
   }

   /* set demand to 1 */
   sm->prob_m[row][col] = 1.0f;
   sm->calc_m[row][col] = 1.0f;

   return sm->fixing_site_hook (data, i, sm);
}

static __inline__ int
seqmatrix_pre_col_iter_hook (void* data, SeqMatrix* sm)
{
   CRB_UNUSED (data);
   CRB_UNUSED (sm);

   return 0;
}

static __inline__ int
seqmatrix_fixed_site_hook (void* data, unsigned long i,SeqMatrix* sm)
{
   CRB_UNUSED (data);
   CRB_UNUSED (sm);
   CRB_UNUSED (i);

   return 0;
}

static __inline__ int 
seqmatrix_fixing_site_hook (void* data, unsigned long i, SeqMatrix* sm)
{
   CRB_UNUSED (data);
   CRB_UNUSED (sm);
   CRB_UNUSED (i);

   return 0;
}

static __inline__ char* 
seqmatrix_get_seq_string (void* data)
{
   CRB_UNUSED (data);

   return NULL;
}

/** @brief Update the columns of a sequence matrix
 *
 * Update cols of a sequence matrix in a SCMF simulation.
 * Returns...
 *
 * @params[in] sm Sequence matrix.
 */
static __inline__ int
seqmatrix_calc_eeff_col_scmf (SeqMatrix* sm,
                              const float t,
                              void* sco)
{
   unsigned long i;
   int error = 0;

   assert (sm);
   assert (sm->calc_eeff_row);

   /* for each col */
   i = 0;
   while ((!error) && (i < sm->cols))
   {
      if (!seqmatrix_is_col_fixed (i, sm))
      {
         error = sm->calc_eeff_row (i, sm, t, sco);
      }
      else
      {
         /*mfprintf (stdout, "\n\n\nIN %lu\n\n\n", i);*/
         error = sm->fixed_site_hook (sco, i, sm);
      }

      i++;
   }

/*    if (sm->curr_matrix == F_Mtrx) */
/*    { */
/*       sm->curr_matrix = S_Mtrx; */
/*    } */
/*    else */
/*    { */
/*       sm->curr_matrix = F_Mtrx; */
/*    } */

   return error;
}

/** @brief Update a row of a sequence matrix
 *
 * Update a row of a sequence matrix in a SCMF simulation.
 * Returns...
 *
 * @params[in] sm Sequence matrix.
 */
static __inline__ int
seqmatrix_calc_eeff_row_scmf (const unsigned long col,
                              SeqMatrix* sm,
                              const float t,
                              void* sco)
{
/*    short int new_matrix = F_Mtrx; */
   unsigned long j;

   assert (sm);
   assert (sm->calc_cell_energy);

   for (j = 0; j < sm->rows; j++)
   {
      sm->calc_m[j][col] = sm->calc_cell_energy (j, col,
                                                 sco,
                                                 sm);

      sm->calc_m[j][col] =
         expf ((-1.0f) * (sm->calc_m[j][col]/(sm->gas_constant*t)));
   }

   return 0;
}

/** @brief Initialise a sequence matrix.
 *
 * On initialisation, all probabilities of the matrix are set. Additionally, a
 * set of sites can be fixed to exclude them from the simulation.\n
 * Returns 0 on success, @c ERR_SM_ALLOC indicating memory problems.
 *
 * @param[in] pairs base pairs.
 * @param[in] size Size of the structure.
 * @param[in] sm Sequence matrix to initialise.
 */
int
seqmatrix_init (const unsigned long rows,
                const unsigned long width,
                SeqMatrix* sm,
                const char* file, const int line)
{
   unsigned long i/* , j */;
   /* float whobble = 0.001f; */

   assert (sm);
   assert (sm->fixed_sites == NULL);
   assert (sm->prob_m      == NULL);
   assert (sm->calc_m      == NULL);

   /* set standard functions */
   sm->calc_eeff_col     = seqmatrix_calc_eeff_col_scmf;
   sm->calc_eeff_row     = seqmatrix_calc_eeff_row_scmf;
   sm->pre_col_iter_hook = seqmatrix_pre_col_iter_hook;
   sm->fixed_site_hook   = seqmatrix_fixed_site_hook;
   sm->fixing_site_hook  = seqmatrix_fixing_site_hook;
   sm->get_seq_string    = seqmatrix_get_seq_string;

   /* allocate matrix */
   sm->rows = rows;
   sm->cols = width;

   sm->prob_m = (float**) XOBJ_MALLOC_2D (sm->rows,sm->cols,
                                          sizeof (**sm->prob_m),
                                          file, line);
   if (sm->prob_m == NULL)
   {
      return ERR_SM_ALLOC;
   }

   sm->calc_m = (float**) XOBJ_MALLOC_2D (sm->rows,sm->cols,
                                          sizeof (**sm->calc_m),
                                          file, line);
   if (sm->calc_m == NULL)
   {
      return ERR_SM_ALLOC;
   }

   /* init sm */
   if (width > 0)
   {
      sm->prob_m[0][0] = 1.0f / sm->rows; /* even distributed init */
      sm->calc_m[0][0] = 0.0f;
      for (i = 1; i < width; i++)
      {
         sm->prob_m[0][i] = sm->prob_m[0][0]; /* even distri */
         sm->calc_m[0][i] = sm->calc_m[0][0];
      }

      for (i = 1; i < sm->rows; i++)
      {
         memcpy (sm->prob_m[i],
                 sm->prob_m[i - 1],
                 sizeof (**sm->prob_m) * width);/* even distri */
         memcpy (sm->calc_m[i],
                 sm->calc_m[i - 1],
                 sizeof (**sm->calc_m) * width);
      }

      /*srand(997654329);
      for (j = 0; j < size; j++)
      {
         for (i = 0; i < SM_ROWS; i++)
         {
            sm->matrix[F_Mtrx][i][j] = rand();
            sm->matrix[S_Mtrx][0][j] += sm->matrix[F_Mtrx][i][j];
         }
         for (i = 0; i < SM_ROWS; i++)
         {
            sm->matrix[F_Mtrx][i][j] = sm->matrix[F_Mtrx][i][j] / sm->matrix[S_Mtrx][0][j];
         }
         sm->matrix[S_Mtrx][0][j] = 0.0f;
         }*/

   }

   /* run over list of presettings and set sites */
   /* allocate memory & init fixed sites */
   sm->fixed_sites = XCALLOC ((width / (sizeof(*(sm->fixed_sites))*CHAR_BIT))+1,
                              sizeof(*(sm->fixed_sites)));

   if (sm->fixed_sites == NULL)
   {
      return ERR_SM_ALLOC;
   }

   return 0;
}

/** @brief Set a certain cell of the effective energy matrix.
 *
 * @param[in] value Value to be set.
 * @param[in] row Row index.
 * @param[in] cell Cell index.
 * @param[in] sm Sequence matrix.
 */
void
seqmatrix_set_eeff (const float value,
                    const unsigned long row,
                    const unsigned long col,
                    SeqMatrix* sm)
{
   assert (sm);
   assert (row < sm->rows);
   assert (col < sm->cols);

   sm->calc_m[row][col] = value;
}

/** @brief Add a number to a certain cell of the effective energy matrix.
 *
 * @param[in] value Value to be set.
 * @param[in] row Row index.
 * @param[in] cell Cell index.
 * @param[in] sm Sequence matrix.
 */
void
seqmatrix_add_2_eeff (const float value,
                      const unsigned long row,
                      const unsigned long col,
                      SeqMatrix* sm)
{
   assert (sm);
   assert (row < sm->rows);
   assert (col < sm->cols);
   /*assert (sm->calc_m[row][col] == 0); just for checking that cell was
     correctly set to 0*/

   sm->calc_m[row][col] += value;
}

/** @brief Set gas constant.
 *
 * @params[in] r Gas constant.
 * @params[in] sm Sequence matrix.
 */
void
seqmatrix_set_gas_constant (const float r, SeqMatrix* sm)
{
   assert (sm);

   sm->gas_constant = r;
}

/** @brief Set the function for calculating columns for the effective energy
 *
 * @params[in] calc_eeff_col Function to use.
 * @params[in] sm Sequence matrix to store the function in.  
 */
void
seqmatrix_set_func_calc_eeff_col (int (*calc_eeff_col) (SeqMatrix*,
                                                        const float,
                                                        void*),
                                  SeqMatrix* sm)
{
   assert (sm);
   sm->calc_eeff_col = calc_eeff_col;
}

/** @brief Set the function for calculating rows for the effective energy
 *
 * @params[in] calc_eeff_row Function to use.
 * @params[in] sm Sequence matrix to store the function in.  
 */
void
seqmatrix_set_func_calc_eeff_row (int (*calc_eeff_row) (const unsigned long,
                                                        SeqMatrix*,
                                                        const float,
                                                        void*),
                                  SeqMatrix* sm)
{
   assert (sm);
   sm->calc_eeff_row = calc_eeff_row;
}

void
seqmatrix_set_func_calc_cell_energy (float (*calc_cell_energy)
                                     (const unsigned long, const unsigned long,
                                      void*,
                                      SeqMatrix*),
                                     SeqMatrix* sm)
{
   assert (sm);
   sm->calc_cell_energy = calc_cell_energy;
}

void
seqmatrix_set_transform_row (int (*transform_row) (const unsigned long,
                                                   const unsigned long,
                                                   void*),
                             SeqMatrix* sm)
{
   assert (sm);
   sm->transform_row = transform_row;
}

void
seqmatrix_set_pre_col_iter_hook (int (*pre_col_iter_hook) (void*, SeqMatrix*),
                                     SeqMatrix* sm)
{
   assert (sm);
   sm->pre_col_iter_hook = pre_col_iter_hook;
}

void
seqmatrix_set_fixed_site_hook (int (*fixed_site_hook) (void*, unsigned long,
                                                       SeqMatrix*),
                                SeqMatrix* sm)
{
   assert (sm);
   sm->fixed_site_hook = fixed_site_hook;
}

void
seqmatrix_set_fixing_site_hook (int (*fixing_site_hook) (void*, unsigned long,
                                                         SeqMatrix*),
                                SeqMatrix* sm)
{
   assert (sm);
   sm->fixing_site_hook = fixing_site_hook;
}

void
seqmatrix_set_get_seq_string (char* (*get_seq_string) (void*), SeqMatrix* sm)
{
   assert (sm);
   sm->get_seq_string = get_seq_string;
}

static __inline__ void
s_seqmatrix_find_lamb_site (SeqMatrix* sm,
                            unsigned long* row,
                            unsigned long* col)
{
   unsigned long i, j;
   float max = 0.0f;

   *col = sm->cols + 1;

   for (j = 0; j < sm->cols; j++)
   {
      /* check if site is fixed */
      if (!seqmatrix_is_col_fixed (j, sm))
      {
         for (i = 0; i < sm->rows; i++)
         {
            if (sm->prob_m[i][j] > max)
            {
               max = sm->prob_m[i][j];
               *col = j;
               *row = i;
            }
         }
      }   
   }
}

/** @brief Transform a sequence matrix into an unambigouos sequence.
 *
 * Creates a sequence out of a sequence matrix in an iterative simulation
 * process. As a start, all unambigouos columns are fixed. Then the matrix is
 * used in a new SCMF process with the "most" unambigouos of the undecided
 * columns fixed by majority voting. This procedure is repeated until... In
 * the end the matrix is compressed into a sequence.
 *
 * @params[in] sm The sequence matrix.
 */
int
seqmatrix_collate_is (const float fthresh,
                      const unsigned long steps,
                      const float temp,
                      const float b_long,
                      const float b_short,
                      const float sc_thresh,
                      const float c_min,
                      const float lambda,
                      /* const float s_thresh, */
                      SeqMatrix* sm,
                      void* data)
{
   unsigned long /*i, j,*/ largest_amb_col = 0, largest_amb_row = 0;
   int retval = 0;

   assert (sm);
   CRB_UNUSED (fthresh);//assert (fthresh <= 1.0f);

   /* Approach: find unambigouos sites and fixate 'em */
   /*           find the largest of the ambigouos sites */
   /*           until all sites are fixed */

   while ((largest_amb_col != sm->cols + 1) && (! retval))
   {
      /* for all columns */
      /* SB: 16-10-09, fixing now happens during the simulation
         for (j = 0; j < sm->cols; j++)
      {
         if (!seqmatrix_is_col_fixed (j, sm))
         {
      *//* for all rows *//*
            for (i = 0; i < sm->rows; i++)
            {
               if (sm->prob_m[i][j] >= fthresh)
               {
                          *//* unambigouos site found, fixate it *//*
                  seqmatrix_fix_col (i, j, sm);

                  i = sm->rows + 1;
               }
            }
         }
      }*/
      
      /* find largest ambigouos site */
      s_seqmatrix_find_lamb_site (sm, &largest_amb_row, &largest_amb_col);

      /* set site to 1/0 and fixate it */
      if (largest_amb_col < sm->cols + 1)
      {
         seqmatrix_fix_col (largest_amb_row, largest_amb_col, data, sm);

         /* simulate */
         /*seqmatrix_print_2_stdout (2, sm);*/
         /* the simulation for collating uses an entropy dropoff of 0, because
            if only a few sites are not fixed, the entropy is to low to trigger
            a simulation step. */
         retval = seqmatrix_simulate_scmf (steps,
                                           temp,
                                           b_long,
                                           b_short,
                                           sc_thresh,
                                           c_min,
                                           /*c_scale,*/
                                           lambda,
                                           0.0f, /* fixed dropoff */
                                           NULL,
                                           NULL,
                                           sm,
                                           data);
      }
   }

   if (!retval)
   {
      retval = seqmatrix_collate_mv (sm, data);
   }

   return retval;
}

/** @brief Transform a sequence matrix into an unambigouos sequence.
 *
 * Collates all rows of a column to a single representative. Thereby the row
 * with the highest share is chosen. The sequence is stored in the SeqMatrix
 * object an can be accessed via...\n
 * Returns @c ERR_SM_ALLOC on problems allocating memory for the sequence. 0
 * otherwise.
 *
 * @params[in] sm The sequence matrix.
 */
int
seqmatrix_collate_mv (const SeqMatrix* sm, void* data)
{
   unsigned long i,j;
   float curr_max_prob;
   unsigned long max_row = 0;

   assert (sm);
   assert (sm->transform_row);

   /* for all columns */
   for (j = 0; j < sm->cols; j++)
   {
      curr_max_prob = -1.0f;
      for (i = 0; i < sm->rows; i++)
      {
         /* find highest number */
         if (curr_max_prob < sm->prob_m[i][j])
         {
            /* write position to seq */
            curr_max_prob = sm->prob_m[i][j];
            max_row = i;
         }
      }
      sm->transform_row (max_row, j, data);
   }

   return 0;
}

static __inline__ float
s_seqmatrix_calc_init_entropy (const SeqMatrix* sm)
{
   unsigned long i, j;
   float s = 0.0f;

   for (j = 0; j < sm->cols; j++)
   {
      if (!seqmatrix_is_col_fixed (j, sm))
      {
         for (i = 0; i < sm->rows; i++)
         {
            if (sm->prob_m[i][j] > FLT_EPSILON)
            {
               s += (sm->prob_m[i][j] * logf (sm->prob_m[i][j]));
            }
         }
      }
   }

   return (s / sm->cols) * (-1.0f);
}

static __inline__
int write_entropy (GFile* file, unsigned long step,
                   float t,
                   float s,
                   float s_short,
                   float s_long,
                   float k)
{
   if (gfile_printf (file, "%lu %f %f %f %f %f %f\n", step,
                     t,
                     s,
                     s_short,
                     s_long,
                     (s_short / s_long),
                     k) < 0)
   {
      return ERR_SM_WRITE;
   }

   return 0;
}

static __inline__
int dummy_write_entropy (GFile* file, unsigned long step,
                         float t,
                         float s,
                         float s_short,
                         float s_long,
                         float k)
{
   CRB_UNUSED (file);
   CRB_UNUSED (step);
   CRB_UNUSED (t);
   CRB_UNUSED (s);
   CRB_UNUSED (s_short);
   CRB_UNUSED (s_long);
   CRB_UNUSED (k);

   return 0;
}

static __inline__
int write_matrix (GFile* file, void* data, const SeqMatrix* sm)
{
   unsigned long i, j; /* row and column indices */
   float max_prob;
   unsigned long max_row = 0;

   assert (file);
   assert (data);
   assert (sm);
   assert (sm->transform_row);

   /* collate seq */
   for (j = 0; j < sm->cols; j++)
   {
      max_prob = -1.0f;

      for (i = 0; i < sm->rows; i++)
      {
         if (max_prob < sm->prob_m[i][j])
         {
            max_prob = sm->prob_m[i][j];
            max_row = i;
         }
      }
      sm->transform_row (max_row, j, data);
   }

   /* write seq as letters */
   if (gfile_printf (file, "%s\n", sm->get_seq_string (data)) < 0)
   {
      return ERR_SM_WRITE;
   }

   /* write matrix */
   for (j = 0; j < sm->cols; j++)
   {
      for (i = 0; i < sm->rows; i++)
      {
         if (gfile_printf (file, " % .6f", sm->prob_m[i][j]) < 0)
         {
            return ERR_SM_WRITE;
         }
      }
      if (gfile_printf (file, "\n") < 0)
      {
         return ERR_SM_WRITE;
      }
   }

   return 0;
}

static __inline__
int dummy_write_matrix (GFile* file,
                        void* data,
                        const SeqMatrix* sm)
{
   CRB_UNUSED (file);
   CRB_UNUSED (data);
   CRB_UNUSED (sm);

   return 0;
}

/** @brief Perform a SCMF simulation on a sequence matrix using the NN.
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
seqmatrix_simulate_scmf (const unsigned long steps,
                         const float t_init,
                         const float b_long,
                         const float b_short,
                         const float sc_thresh,
                         const float c_min,
/*                         const float c_scale __attribute__((unused)),*/
                         const float lambda,
                         const float s_thresh,
                         GFile* entropy_file,
                         GFile* matrix_file,
                         SeqMatrix* sm,
                         void* sco)
{
   unsigned long t = 0;         /* time */
   int error= 0;
   unsigned long i, j;          /* iterator */
   float col_sum;
   float T = t_init;            /* current temperature */
   float c_rate = 0.999999f;/* SB 090715 for testing 1.0f; *//* cooling rate */
   float s_cur/* , s_last */;         /* matrix entropy */
   /* unsigned long s_count; */       /* count times s did not change */
   /* unsigned long s_dropout;  */    /* convergence criterion */
   /* unsigned long largest_amb_row = 0, largest_amb_col; */
   float s_long, s_short;
   int (*output_entropy) (GFile*,
                          unsigned long,
                          float,
                          float,
                          float,
                          float,
                          float) = dummy_write_entropy;
   int (*output_matrix) (GFile*,
                         void*,
                         const SeqMatrix*) = dummy_write_matrix;

   assert (sm);
   assert (sm->calc_eeff_col);
   assert (sco);

   /* init. long and short term avg. entropies */
   s_cur = s_short = s_seqmatrix_calc_init_entropy (sm);
   s_long = s_short * 2;     /* SB 25-11-09, was s_long = s_short */

   /* calculate initial cooling rate */
/*  SB for testing, 2009-03-30  if (steps > 0) */
/*    { */
/*       c_rate = logf (t_init / 1.0f); */
/*       c_rate /= (steps - 1); */
/*       c_rate = expf ((-1) * c_rate); */
/*    } */

   /* if we have an output file, we set a function to write stats */
   if (entropy_file != NULL)
   {
      output_entropy = write_entropy;

      /* write info on simulation */
      if (gfile_printf (entropy_file, "# step | T | S | S_short | S_long | "
                        "(S_short / S_long) | cooling rate\n") < 0)
      {
         error = 1;
      }
      else if (output_entropy (entropy_file,
                               t,
                               T,
                               s_cur,
                               s_short,
                               s_long,
                               c_rate) < 0)
      {
         error = 1;
      }
   }

   if (matrix_file != NULL)
   {
      output_matrix = write_matrix;
   }

/*    s_last = FLT_MAX * (-1.0f); */
/*    s_count = 0; */
/*    s_dropout = 100/\* steps * 0.005 *\/; */

   /* perform for a certain number of steps */
   /* SB 16-09-09 T > 1.0f */
   while ((!error) && (t < steps) && (T > 0.45f) && (s_cur >= s_thresh))
   {
      /* mfprintf (stderr, "Step: %lu\n", t); */
      error = sm->pre_col_iter_hook (sco, sm);

      s_cur = 0.0f;

      /* calculate Eeff */
      if (!error)
      {
         error = sm->calc_eeff_col (sm, T, sco);
      }

      /* update matrix */
      /* for all columns */ 
      if (!error)
      {
         for (j = 0; j < sm->cols; j++)
         {
            /* which are not fixed */
            if (!seqmatrix_is_col_fixed (j, sm))
            {
               /* calc sum of col */
               col_sum = 0.0f;
               for (i = 0; i < sm->rows; i++)
               {
                  col_sum += sm->calc_m[i][j];
               }

               /* for each row */
               for (i = 0; i < sm->rows; i++)
               {
                  sm->calc_m[i][j] = sm->calc_m[i][j] / col_sum;  

                  /* avoid oscilation by Pnew = uPcomp + (1 - u)Pold) */
                  sm->prob_m[i][j] = 
                     (lambda * sm->calc_m[i][j])
                     + ((1 - lambda) * sm->prob_m[i][j]);

                  if (sm->prob_m[i][j] > 0.99f)
                  {
                     seqmatrix_fix_col (i, j, sco, sm);
                     i = sm->rows;
                     /* if this works, check if we can omit fixing in
                        collate_is (we still need to search for the next one to
                        fix but do not need to search for > 0.99!) */
                  }
                  else if (sm->prob_m[i][j] > FLT_EPSILON)
                  {
                     /* calculate "entropy", ignore fixed sites since ln(1)=0 */
                     s_cur += (sm->prob_m[i][j] * logf (sm->prob_m[i][j]));
                  }
               }
            }
         }

         s_cur = (s_cur / sm->cols) * (-1.0f);

         /* cooling: We follow the entropy with a long- and short term avg. If
            s_long and s_short diverge to much, we slow down or speed up
            cooling. */
         s_long  = (b_long * s_long) + ((1 - b_long) * s_cur);
         s_short = (b_short * s_short) + ((1 - b_short) * s_cur);

         if ((s_short / s_long) < sc_thresh)
         {
            /* entropy changes to fast, slow down */
            c_rate = sqrtf (c_rate);
            if (c_rate >= 1.000000f)
            {
               c_rate = 0.999999f;
            }
         }
         else
         {
            /* small changes, speed up */
            if (c_rate > c_min)
            {
               /* SB 090714 for testing c_rate = c_rate * (c_rate * c_scale);*/
               c_rate = c_rate * c_rate; /* SB 090714 for testing */
               /* c_rate = c_rate * 0.95f; working for 192 */
            }
         }

         T = T * c_rate;
         t++;

         error = output_entropy (entropy_file, t, T, s_cur, s_short, s_long,
                                 c_rate);
         if (!error)
         {
            error = output_matrix (matrix_file, sco, sm);
         }
      }
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
   float tmp = 0.0f;
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
         if (sm->prob_m[i][j] < 0.0f)
         {
            tmp = (-1) * sm->prob_m[i][j];
            rprec = 1;
         }
         else
         {
            tmp = sm->prob_m[i][j];
            rprec = 0;
         }

         if (tmp > 1.0f)
         {
            rprec += floorf (log10f (tmp) + 1.0f);
         }
         else
         {
            rprec += 1;
         }

         if (rprec > cprec) 
         {
            cprec = rprec;
         }      
      }
   }
   cprec += p + 1;

   /* SB 16-09-09 - Start */
   if ((sm->cols > 1) && (floorf (log10f ((float) sm->cols) + 1.0f) >= cprec))
   {
      cprec = floorf (log10f ((float) sm->cols) + 1.0f) + 1;
   }
   /* SB 16-09-09 - End */

   if (sm->rows > 1) /* SB 16-09-09 1.0f*/
   {
      rprec = floorf (log10f ((float) sm->rows) + 1.0f);
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
   /* SB 16-09-09 string = XMALLOC (sizeof (*string) * (((line_width) * sm->rows) + 1));*/
   string = XMALLOC (sizeof (*string) * (((line_width) * (sm->rows + 1)) + 1)); /* SB 16-09-09*/
   string_start = string;

   /* SB 16-09-09 - Start */
   /* printing col. no.s: 
      - allocate 1 line extra string-space
      - skip to first '|'?
      - width/2 as a spacer
      - what about no. of digits? */
   /* write col. no.s */
   if (sm->cols > 0)
   {
      msprintf (string, "%*c", (rprec + 4), ' ');
      string += rprec + 4;

      for (j = 0; j < sm->cols; j++)
      {
         msprintf (string, " %*lu   ", (cprec - 1), j);
         string += cprec + 3;
      }
      msprintf (string, "\n");
      string++;
   }
   /* SB 16-09-09 - End */

   /* write matrix */
   for (i = 0; i < sm->rows; i++)
   {
      msprintf (string, "%*lu:  | ", rprec, i);
      string += rprec + 4;

      for (j = 0; j < sm->cols; j++)
      {
         msprintf (string, " %*.*f |",
                   cprec,
                   p,
                   sm->prob_m[i][j]);
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
