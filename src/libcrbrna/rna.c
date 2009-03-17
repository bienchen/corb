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
 *  @file libcrbrna/rna.c
 *
 *  @brief RNA datastructure
 *
 *  Module: rna
 *
 *  Library: crbrna
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-08-19
 *
 *
 *  Revision History:
 *         - 2008Aug19 bienert: created
 *
 */


#include <config.h>
#include <stdlib.h>
#include <limits.h>
#include <stddef.h>
#include <libcrbbasic/crbbasic.h>
/*#include "alphabet.h"*/
#include "secstruct.h"
/*#include "nn_scores.h"*/
#include "rna.h"


/*typedef struct {
      ArrayUlong tetra_loop;
      } struct_comp;*/


struct Rna {
   char*          seq;       /* the nucleotide sequence */
   unsigned long  size;      /* size of the RNA sequence */
   /* char* vienna; */       /* vienna string */
   unsigned long* pairs;     /* base pairs */
   unsigned long pairs_size; /* size of the pairlist */
   SecStruct*    structure; /* decomposed structure */
   Str* info;
};

#define CT_FILE_ERR "ct-file \'%s\', line %lu: "

/**********************   Constructors and destructors   **********************/

/** @brief Create a new rna object.
 *
 * The constructor for @c Rna objects. If compiled with enabled
 * memory checking, @c file and @c line should point to the position where the
 * function was called. Both parameters are automatically set by using the
 * macro @c RNA_NEW.\n
 * Returns @c NULL on error.
 *
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
Rna*
rna_new (const char* file, const int line)
{
   Rna* this = XOBJ_MALLOC(sizeof (*this), file, line);
   
   if (this != NULL)
   {
      this->size      = 0;
      this->seq       = NULL;
      /* this->vienna    = NULL; */
      this->pairs     = NULL;
      this->pairs_size = 0;
      this->structure = NULL;
      this->info = NULL;
   }

   return this;
}

/** @brief Delete a Rna object.
 *
 * The destructor for @c Rna objects.
 *
 * @param[in] this Object to be freed.
 */
void
rna_delete (Rna* this)
{

   if (this != NULL)
   {
     XFREE(this->seq);
     this->size = 0;
     /* XFREE(this->vienna); */
     XFREE(this->pairs);
     this->pairs_size = 0;
     secstruct_delete (this->structure);
     str_delete (this->info);
     XFREE(this);
   }
}

/** @brief Allocate memory for the sequence component of an @c Rna object.
 *
 * Allocate the memory for a sequence stored within an @c Rna data object. If
 * compiled with enabled memory checking, @c file and @c line should point to
 * the position where the function was called. Both parameters are 
 * automatically set by using the macro @c RNA_ALLOCATE_SEQUENCE.\n
 * Returns 0 on success, @c ERR_RNA_ALLOC on memory problems.
 *
 * @param[in] size Size of the sequence/ structure of the RNA.
 * @param[in] this Rna data object.
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
int
rna_alloc_sequence (const unsigned long size, Rna* this,
                    const char* file, const int line)
{
   assert (this);
   assert (this->seq == NULL);
   assert ((this->size == 0) || (this->size == size));

   this->size = size;

   this->seq = XOBJ_CALLOC ((this->size + 1), sizeof (this->seq[0]),
                            file, line);
   if (this->seq == NULL)
   {
      return ERR_RNA_ALLOC;
   }

   return 0;
}

/** @brief Reallocate memory for the sequence component of an @c Rna object.
 *
 * Reallocate the memory for a sequence stored within an @c Rna data object. If
 * compiled with enabled memory checking, @c file and @c line should point to
 * the position where the function was called. Both parameters are 
 * automatically set by using the macro @c RNA_REALLOC_SEQUENCE.\n
 * Returns 0 on success, @c ERR_RNA_ALLOC on memory problems.
 *
 * @param[in] size Size of the sequence/ structure of the RNA.
 * @param[in] this Rna data object.
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
int
rna_realloc_sequence (const unsigned long size, Rna* this,
                        const char* file, const int line)
{
   assert (this);

   this->seq = XOBJ_REALLOC(this->seq, (size + 1) * sizeof (this->seq[0]),
                            file, line);

   if (this->seq == NULL)
   {
      return ERR_RNA_ALLOC;
   }

   if (size > this->size)
   {
      memset ((this->seq + this->size), 0,
              sizeof (this->seq[0]) * ((size + 1) - this->size));   
   }

   this->size = size;

   return 0;
}


/** @brief Store a copy of an RNA sequence in an Rna data object.
 *
 * Initialise the sequence component of an Rna object. The data source is a
 * character string over a certain alphabet. The sequence is transforme into
 * the internal RNA representation and therefore stored as a copy in the
 * object. If compiled with enabled memory checking, @c file and @c line
 * should point to the position where the function was called. Both parameters
 * are automatically set by using the macro @c RNA_INIT_SEQUENCE.\n
 * Returns 0 on success, @c ERR_RNA_ALLOC on memory problems. On problems
 * parsing the sequence string, the following error codes are returned:\n
 * @list
 * @c ERR_RNA_NO_BASE A base in the sequence is not in the given alphabet.
 *
 * @params[in] seq Sequence string.
 * @params[in] length Length of the sequence.
 * @params[in] sigma Alphabet.
 * @params[in] this Rna object to store structure.
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
int
rna_init_sequence (const char* seq,
                   const unsigned long length,
                   Alphabet* sigma,
                   Rna* this,
                   const char* file, const int line)
{
   int error = 0;
   unsigned long i;

   assert (this);
   assert (this->seq == NULL);
   assert ((this->size == 0) || (this->size == length));

   /* allocate sequence */
   error = rna_alloc_sequence (length, this, file, line);

   /* store & validate sequence */
   i = 0;
   while ((!error) && (i < length))
   {
      this->seq[i] = alphabet_base_2_no (seq[i], sigma);

      if (this->seq[i] == CHAR_UNDEF)
      {
         error =  ERR_RNA_NO_BASE;
      }
      i++;
   }

   return error;
}

/** @brief Allocate memory for the pair list component of an @c Rna object.
 *
 * Allocate the memory for a list of pairs stored within an @c Rna data object.
 * If compiled with enabled memory checking, @c file and @c line should point to
 * the position where the function was called. Both parameters are 
 * automatically set by using the macro @c RNA_ALLOCATE_PAIRLIST.\n
 * Returns 0 on success, @c ERR_RNA_ALLOC on memory problems.
 *
 * @param[in] size Size of the sequence/ structure of the RNA.
 * @param[in] this Rna data object.
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
int
rna_allocate_pairlist (const unsigned long size, Rna* this,
                       const char* file, const int line)
{
   assert (this);
   assert (this->pairs == NULL);
   assert (this->pairs_size == 0);  
   
   this->pairs_size = size;
   this->pairs = XOBJ_MALLOC(sizeof (this->pairs[0]) * this->pairs_size,
                             file, line);

   if (this->pairs == NULL)
   {
      return ERR_RNA_ALLOC;
   }

   memset (this->pairs, INT_MAX, sizeof (this->pairs[0]) * this->pairs_size);   

   return 0;
}

/** @brief Reallocate memory for the pair list of an @c Rna object.
 *
 * Reallocate the memory for the list of base pairs stored within an @c Rna
 * data object. If compiled with enabled memory checking, @c file and @c line
 * should point to the position where the function was called. Both parameters
 * are automatically set by using the macro @c RNA_REALLOC_PAIRLIST.\n
 * Returns 0 on success, @c ERR_RNA_ALLOC on memory problems.
 *
 * @param[in] size Size of the structure.
 * @param[in] this Rna data object.
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
int
rna_realloc_pairlist (const unsigned long size, Rna* this,
                      const char* file, const int line)
{
   assert (this);

   this->pairs = XOBJ_REALLOC(this->pairs, size * sizeof (this->pairs[0]),
                            file, line);

   if (this->pairs == NULL)
   {
      return ERR_RNA_ALLOC;
   }

   if (size > this->pairs_size)
   {
      memset ((this->pairs + this->pairs_size), INT_MAX,
              sizeof (this->pairs[0]) * (size - this->pairs_size));
   }

   this->pairs_size = size;

   return 0;
}

/** @brief Read the list of base pairs from a Vienna string.
 *
 * Initialise the base pairs component of an Rna object. The data source is a
 * structure string in Vienna notation which defines size of the and pairs in
 * the list. Position counting starts by 1. Non-paired positions will get a
 * value of @c NOT_PAIRED. Additionally the Vienna string itself is copied to
 * the data object. If compiled with enabled memory checking, @c file and
 * @c line should point to the position where the function was called. Both
 * parameters are automatically set by using the macro
 * @c RNA_INIT_PAIRLIST_VIENNA.\n
 * Returns 0 on success, @c ERR_RNA_ALLOC on memory problems. For problems
 * parsing the Vienna string, the following error codes are resrved:\n
 * @list
 * @c ERR_RNA_VIENNA_FORMAT if the string contains unknown symbols
 * @c ERR_RNA_VIENNA_MMC    if there are more closing then opening paring
 *                          partners
 * @c ERR_RNA_VIENNA_MMO    if there are more opening then closing paring
 *                          partners
 *
 * @params[in] vienna Structure string.
 * @params[in] length Length of @c vienna.
 * @params[in] this Rna object to store structure.
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
int
rna_init_pairlist_vienna (const char* vienna,
                          const unsigned long length,
                          Rna* this,
                          const char* file, const int line)
{
   int error = 0;
   unsigned long i;
   unsigned long open_pair;
   unsigned long* p_stack = NULL;

   assert (this);
   /* assert (this->vienna == NULL); */
   assert (this->pairs == NULL);
   assert (this->pairs_size == 0);

   /* allocate list */
   error = rna_allocate_pairlist (length, this, file, line);

   if (!error)
   {
      this->size = length;
      p_stack = XCALLOC (this->size, sizeof (p_stack[0]));
      if (p_stack == NULL)
      {
         error = ERR_RNA_ALLOC;
      }
   }

   /* transform structure string to base pair list */
   i = 0;
   open_pair = 0;
   while ((!error) && (i < this->size))
   {
      switch (vienna[i])
      {
         case VIENNA_OPEN :
            p_stack[open_pair] = i;
            open_pair++;
            break;     
            
         case VIENNA_CLOSE :
            if (open_pair == 0)
            {
               THROW_ERROR_MSG ("Mismatched nucleotide (closing base pair "
                                "partner) found at position %lu in structure "
                                "string (\'%s\')", i, vienna);
               error = ERR_RNA_VIENNA_MMC;               
            }
            else
            {
               open_pair--;
               this->pairs[i] = p_stack[open_pair];
               this->pairs[p_stack[open_pair]] = i;
            }
            break;
            
         case VIENNA_UNPAIRED : ;
            break;     
            
         default  : 
            THROW_ERROR_MSG ("Non-valid symbol found in structure string "
                             "(\'%s\'). Allowed characters: \'%c\', \'%c\', "
                             "\'%c\'. Found \'%c\' at position %lu.",
                             vienna,
                             VIENNA_OPEN,
                             VIENNA_CLOSE,
                             VIENNA_UNPAIRED,
                             vienna[i], i);
            error = ERR_RNA_VIENNA_FORMAT;
            break;
      }

      i++;
   }

   if (open_pair != 0)
   {
      THROW_ERROR_MSG ("Mismatched nucleotide (opening base pair "
                       "partner) found in structure string (\'%s\')", vienna);
      error = ERR_RNA_VIENNA_MMO;
   }
   
   XFREE (p_stack);

   if (error)
   {
      XFREE (this->pairs);
      this->pairs = NULL;
      this->pairs_size = 0;
   }

   return error;
}

/** @brief Initialise an Rna data object with a sequence and a structure.
 *
 * Store sequence and structure in an Rna data object. Does the same as using
 * @c rna_init_sequence together with @c rna_init_pairlist_vienna. Since there
 * is only one @c length parameter, we assume that sequence and structure are
 * of equal length. If compiled with enabled memory checking, @c file and
 * @c line should point to the position where the function was called. Both
 * parameters are automatically set by using the macro
 * @c RNA_INIT_SEQUENCE_STRUCTURE.\n
 * Returns 0 on success or an error code from the utilised functions.
 *
 * @param[in] seq Sequence.
 * @param[in] vienna Structure string in Vienna notation.
 * @param[in] length Length of the RNA.
 * @param[in] sigma Alphabet.
 * @param[in] this Rna data object.
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
int
rna_init_sequence_structure (const char* seq,
                             const char* vienna,
                             const unsigned long length,
                             Alphabet* sigma,
                             Rna* this,
                             const char* file, const int line)
{
   int error = 0;

   /* allocate and store sequence */
   error = rna_init_sequence (seq, length, sigma, this, file, line);

   /* allocate and store structure */
   if (!error)
   {
      error = rna_init_pairlist_vienna (vienna, length, this, file, line);
   }

   return error;
}

/** @brief Read a RNA structure and sequence from file
 *
 * Read in a file in connect format and store its contents in a Rna object.\n
 * Returns 0 on success, ... else.
 * GFILE_READ_ERROR on problems reading next line
 * ERR_RNA_ALLOC info component can't be stored
 * @param[in] this Rna object.
 * @param[in] path Path and file name.
 */
static unsigned long
rna_get_first_line_ct (char** buffer, size_t* length,
                       unsigned long* line_no, GFile* gfile)
{
   unsigned long lwp = 0;       /* last whitespace position */
   unsigned long col = 0;
   char* endptr;                /* needed for strtoul */
   bool ct_line = true;         /* is the first line already in ct format? */
   unsigned long n_bases = 1;   /* no. of bases */
   unsigned long no;            /* just the number in a col */
   unsigned long max_no = 0;    /* store max.no. found */
   int error = 0;
   unsigned long line_length;
   unsigned long i;

   /* Basic idea: Verify that first line is either a comment or the start of
      the base list */
   /* After this function gfile should be at the starting position of the list*/

   line_length = gfile_getline (&error, buffer, length, gfile);
   if (error)
   {
      return 0;
   }

   /* We run up to the '\0' at the end of the line, translate it to ' ' and
      fetch the last number in a ct line correctly. After the loop we restore
      the '\0'. */
   for (i = 0; i < line_length; i++)
   {
      /* transform \0 in line into ' ' */
      if ((*buffer)[i] == '\0')
      {
         (*buffer)[i] = ' ';
      }

      if ((ct_line == true) && ((*buffer)[i] == ' '))
      {
         if ((i > 0) && ((*buffer)[i - 1] != ' '))
         {
            col++;

            if (col != 2)
            {
               no = strtoul((*buffer + lwp), &endptr, 10);

               /* we have read a number, if pointers differ AND the the whole
                  part of the buffer we are looking at was read */
               if ((*buffer == endptr) || (((*buffer + i) - endptr) != 0))
               {
                  ct_line = false;

                  if (col == 1)
                  {
                     THROW_WARN_MSG ("File \"%s\" in ct format: Header does "
                                     "not start with no. of bases",
                                     str_get(gfile_get_path (gfile)));
                  }
               }
               else if (col == 1)
               {
                  n_bases = no;
               }
               else
               {
                  if (no > max_no)
                  {
                     max_no = no;
                  }
               }
            }
            else                /* col 2: Only non-numerical col */
            {
               /* we do not check the alphabet here */
               /* but the nucleotide col must have width 1 */
               if ((*buffer)[i - 2] != ' ')
               {
                  ct_line = false;
               }
            }

            /* remember position of last whitespace */
            lwp = i + 1;
         }
      }
   }
   if (line_length > 0)
   {
      (*buffer)[line_length - 1] = '\0';
   }

   *line_no = 1;

   if ((ct_line == true) && (col == 6))
   {
      if (max_no > n_bases)
      {
         n_bases = max_no;
      }

      THROW_WARN_MSG ("File \"%s\" in ct format: No header found in first line",
                      str_get(gfile_get_path (gfile)));

      *line_no = 0;

      error = gfile_rewind (gfile);

      if (error)
      {
         return 0;
      }
   }

   return n_bases;
}

/* read a ct line from file */
enum {
   Seq_pos = 0,                     /* position in sequence */
   p5_con,                          /* partner in 5' direction */
   p3_con,                          /* partner in 3' direction */
   partner,                         /* pairing partner of a base */
   Pos_seq,                         /* again position in sequence */
   N_ct_nos
};

static __inline__ unsigned long
s_rna_scan_line_ct (unsigned long cols[N_ct_nos],
                    char* base,
                    unsigned long* line_no,
                    int* error,
                    char** buf,
                    size_t* buf_size,
                    GFile* gfile)
{
   unsigned long i;
   unsigned long read;
   unsigned long col = 0;
   unsigned long st = 0;        /* stored number */
   unsigned long no;            /* number read */
   unsigned long lwp = 0;       /* pos of last ws */
   char* endptr;

   /* try to get a new line */
   read = gfile_getline (error, buf, buf_size, gfile);
   (*line_no)++;

   /* we do not check errors here because if(error) then read = 0 */
   while (read == 1)            /* only '\n' read */
   {
      THROW_WARN_MSG (CT_FILE_ERR "Empty line.",
                      str_get(gfile_get_path (gfile)),
                      *line_no);
      read = gfile_getline (error, buf, buf_size, gfile);

      (*line_no)++;
   }

   /* skip all ws at beginning of line */
   i = 0;
   while ((i < read) && ((*buf)[i] == ' '))
   {
      i++;
   }

   if (i > 0)
   {
      THROW_WARN_MSG (CT_FILE_ERR "Leading whitespaces.",
                      str_get(gfile_get_path (gfile)),
                      *line_no);      
   }

   /* changing the terminal '\0' to ' ' saves us 1 if-statement in the loop.
      Since we do not need the buf in the calling function, we do not care
      restoring '\0'.*/
   if (read > 0)
   {
      (*buf)[read - 1] = ' ';
   }
   else
   {
      /* In case we have read EOF (read == 0), we have not read a line but
         counted. Testing here instead of right after reading saves us 1
         if-statemnt. */
      (*line_no)--;
   }

   /* While reading we just store values read. We check for
      - the right symbols: No. or one letter chars
      - first and last col (seq.pos.) must contin same values
      here. */
   for (; i < read; i++)
   {
      if ((*buf)[i] == ' ' /* || (*buf)[i] == '\0' */)
      {
         /* ws found, let's check if we are behind a symbol */
         if ((*buf)[i - 1] != ' ')
         {
            /* col 2 is the only non-numerical */
            if (col != 1)
            {
               no = strtoul((*buf + lwp), &endptr, 10);

               /* we have read a number, if pointers differ AND the the whole
                  part of the buffer we are looking at was read */
               if ((*buf == endptr) || (((*buf + i) - endptr) != 0))
               {
                  (*buf)[i] = '\0';
                  THROW_ERROR_MSG (CT_FILE_ERR "Column %lu does not contain "
                                   "a number: \"%s\".",
                                   str_get(gfile_get_path (gfile)),
                                   *line_no,
                                   (col + 1),
                                   (*buf + lwp));
                  *error = ERR_RNA_CT_NAN;
                  return 0;
               }

               cols[st] = no;
               st++;
            }
            else /* store the base */
            {
               if ((*buf)[i - 2] != ' ')
               {
                  (*buf)[i] = '\0';
                  THROW_ERROR_MSG (CT_FILE_ERR "Column %lu is not a single "
                                   "letter nucleotide: \"%s\".",
                                   str_get(gfile_get_path (gfile)),
                                   *line_no,
                                   (col + 1),
                                   (*buf + lwp));
                  *error = ERR_RNA_CT_NN;
                  return 0;
               }

               *base = (*buf)[i - 1];
            }

            col++;
         }
         lwp = i + 1;
      }
   }

   /* do line syntax/ semantic check */
   if (cols[Seq_pos] != cols[Pos_seq])
   {
      THROW_ERROR_MSG (CT_FILE_ERR "First and last column, supposed to "
                       "redundantly describe the sequence position, do not "
                       "match: %lu != %lu.\n",
                       str_get(gfile_get_path (gfile)),
                       *line_no,
                       cols[Seq_pos],
                       cols[Pos_seq]);
      *error = ERR_RNA_CT_SM;
      col = 0;
   }
   if (cols[Seq_pos] == 0)
   {
      THROW_ERROR_MSG (CT_FILE_ERR "Sequence position is \'0\', bases are "
                       "counted beginning with \'1\'.\n",
                       str_get(gfile_get_path (gfile)),
                       *line_no);
      *error = ERR_RNA_CT_SM;
      col = 0;
   }
   if (cols[p5_con] == cols[Seq_pos])
   {
      THROW_ERROR_MSG (CT_FILE_ERR "5' connection (\'%lu\') points to the "
                       "sequence position (\'%lu\') of the current base.\n",
                       str_get(gfile_get_path (gfile)),
                       *line_no,
                       cols[p5_con],
                       cols[Seq_pos]);
      *error = ERR_RNA_CT_SM;
      col = 0;
   }
   if (cols[p3_con] == cols[Seq_pos])
   {
      THROW_ERROR_MSG (CT_FILE_ERR "3' connection (\'%lu\') points to the "
                       "sequence position (\'%lu\') of the current base.\n",
                       str_get(gfile_get_path (gfile)),
                       *line_no,
                       cols[p3_con],
                       cols[Seq_pos]);
      *error = ERR_RNA_CT_SM;
      col = 0;
   }

   return col;
}

/* move to rna_new_from_file_ct later */
/* reformat all msgs to "ct-file '.ct', line n: ..." */
int
rna_read_from_file_ct (Rna* this, GFile* gfile)
{
   int error = 0;
   char* line_buffer = NULL;
   size_t lb_size = 0;
   unsigned long n;
   unsigned long line_no = 0;       /* current line number */
   unsigned long ct_cols[N_ct_nos];
   unsigned long n_bases = 0;
   char base = 0;

   assert(this);
   assert(this->info == NULL);
   assert(gfile);

   /* fetch first line */
   n = rna_get_first_line_ct (&line_buffer, &lb_size, &line_no, gfile);
   if (!n)
   {
      error = GFILE_READ_ERROR;
   }

   /* store first line as info component */
   if (!error)
   {
      this->info = STR_NEW_CSTR(line_buffer);
      if (this->info == NULL)
      {
         error = ERR_RNA_ALLOC;
      }
   }

   /* init. allocation of rna structue */
   if (!error)
   {
      error = rna_allocate_pairlist (n, this, __FILE__, __LINE__);
   }
   if (!error)
   {
      error = rna_alloc_sequence (n, this, __FILE__, __LINE__);
   }

   /* read file */
   /* read a line into an array, store vals from array in Rna */
   while ((!error) && (s_rna_scan_line_ct (ct_cols,
                                           &base,
                                           &line_no,
                                           &error,
                                           &line_buffer,
                                           &lb_size,
                                           gfile) == (N_ct_nos + 1)))
   {
      n_bases++;

      /* check whether we have to reallocate */
      if ((n < ct_cols[Seq_pos]) || (n < ct_cols[partner]))
      {
         /* We are generous and give the Rna storage for 1000 new basese. Of
            course there could be a n. in that line telling us we need more
            space but usually we don't... */
         n = (n + 1000 > ct_cols[partner] ? n + 1000 : ct_cols[partner]);
         error = rna_realloc_sequence (n, this, __FILE__, __LINE__);
         if (!error)
         {
            error = rna_realloc_pairlist (n, this, __FILE__, __LINE__);
         }
      }

      if (!error)
      {
         /* store lines as Rna */
         this->seq[ct_cols[Seq_pos] - 1] = base;

         if (ct_cols[partner])
         {
            if (this->pairs[ct_cols[Seq_pos] - 1] == NOT_PAIRED)
            {
               this->pairs[ct_cols[Seq_pos] - 1] = ct_cols[partner] - 1;

               /* assure that partner is either unpaired or paired to us */
               if ((this->pairs[ct_cols[partner] - 1] != NOT_PAIRED)
                  &&(this->pairs[ct_cols[partner] - 1] != ct_cols[Seq_pos] - 1))
               {
                  THROW_ERROR_MSG (CT_FILE_ERR "Base %lu pairs with base %lu "
                                   "which is already paired with base %lu.",
                                   str_get(gfile_get_path(gfile)),
                                   line_no,
                                   ct_cols[Seq_pos] - 1,
                                   ct_cols[partner] - 1,
                                   this->pairs[ct_cols[partner] - 1]);
                  error = ERR_RNA_CT_SM;
               }
            }
            else
            {
               THROW_ERROR_MSG (CT_FILE_ERR "Base %lu should pair with base "
                                "%lu but is already paired to base %lu.\n",
                                str_get(gfile_get_path (gfile)),
                                line_no,
                                ct_cols[Seq_pos] - 1,
                                ct_cols[partner] - 1,
                                this->pairs[ct_cols[Seq_pos] - 1]);
               error = ERR_RNA_CT_SM;
            }
         }
      }
   }

   /* reallocate to true size */
   if (n > n_bases)
   {
      if (!error)
      {
         error = rna_realloc_sequence (n_bases, this, __FILE__, __LINE__);
      }
      if (!error)
      {
         error = rna_realloc_pairlist (n_bases, this, __FILE__, __LINE__);
      }
   }
   XFREE(line_buffer);

   return error;
}

int
rna_read_from_file (Rna* this, const char* path, const unsigned long length,
                    const char* fext, const unsigned long fext_len,
                    const GFileType ctype)
{
   int error = 0;
   const char* exts[] = { "ct" };
   unsigned long n_exts = sizeof (exts) / sizeof (&exts);
   unsigned long entry;
   GFile* gfile;
   unsigned long i;
   int __attribute__((unused)) (*read_funcs[1])(Rna*, GFile*) = {
      &rna_read_from_file_ct
   };

   assert (this);
   assert ((fext && fext_len) || (!fext && !fext_len));

   /* determine format */
   if (fext_len == 0)
   {
      entry = gfile_ext_from_list (path, length, exts, n_exts);

      if (entry == n_exts)
      {
         THROW_ERROR_MSG ("Extension of file \"%s\" not recognised.", path);
         return ERR_RNA_EXT_NR;
      }
   }
   else
   {
      for (i = 0; i < n_exts; i++)
      {
         if (strncmp (fext, exts[i], fext_len) == 0)
         {
            entry = i;
            i = n_exts;
         }
      }

      if (i == n_exts)
      {
         THROW_ERROR_MSG ("Extension \"%s\" not recognised.", fext);
         return ERR_RNA_EXT_NR;
      }   
   }

   /* open file */
   gfile = GFILE_OPEN(path, length, ctype, "r");
   if (gfile == NULL)
   {
      return ERR_RNA_FILE_OPEN;
   }

   /* read file */
   error = read_funcs[entry] (this, gfile);

   /* close file */
   if (!error)
   {
      error = gfile_close(gfile);
   }
   else
   {
      gfile_close(gfile);
   }

   return error;
}

/** @brief Decompose an already stored secondary structue.
 *
 * Initialise the secondary structure component of an Rna data object. Instead
 * of a simple list of pairs, this creates a set of structural components the
 * RNA structure consists of. If compiled with enabled memory checking,
 * @c file and @c line should point to the position where the function was
 * called. Both parameters are automatically set by using the macro
 * @c RNA_SECSTRUCT_INIT.\n
 * Returns 0 on success, @c ERR_RNA_ALLOC if memory allocation failed.
 *
 * @param[in] this Rna object.
 * @param[in] file Fill with name of calling file.
 * @param[in] line Fill with calling line.
 */
int
rna_secstruct_init (Rna* this, const char* file, const int line)
{
   int error = 0;

   assert (this);
   assert (this->pairs);
   assert (this->structure == NULL);

   this->structure = secstruct_new (file, line);

   if (this->structure == NULL)
   {
      error = ERR_RNA_ALLOC;
   }
   else
   {
      error = secstruct_find_interactions (this->pairs, this->size,
                                           this->structure);
   }

/*    if (!error) */
/*    { */
/*       mprintf ("Stacked base piars:\n"); */
/*       secstruct_fprintf_stacks (stdout, this->structure); */
/*       mprintf ("\nHairpin loops:\n"); */
/*       secstruct_fprintf_hairpins (stdout, this->structure); */
/*       mprintf ("\nBulge loops:\n"); */
/*       secstruct_fprintf_bulges (stdout, this->structure); */
/*       mprintf ("\nInternal loops:\n"); */
/*       secstruct_fprintf_internals (stdout, this->structure); */
/*       mprintf ("\nExternal loop:\n"); */
/*       secstruct_fprintf_external (stdout, this->structure); */
/*       mprintf ("\nMulti loops:\n"); */
/*       secstruct_fprintf_multiloops (stdout, this->structure); */
/*    } */

   return error;
}

/********************************   Altering   ********************************/

/** @brief Set a certain base at a certain position in an Rna sequence.
 *
 * Set a given base at a certain position in the sequence component of an Rna
 * object. Position counting starts at 0.\n
 *
 * @param[in] base The base type to set.
 * @param[in] pos Position where to place base.
 * @param[in] this Rna data object.
 */
void
rna_set_sequence_base (const char base, const unsigned long pos, Rna* this)
{
   assert (this);
   assert (this->seq);
   assert (this->size > pos);

   this->seq[pos] = base;

}

/** @brief Set a pairing partner for a base.
 *
 * Set a certain sequence position to be the pairing partner of another base.
 * Ony one base is modified here. To set a complete base pair you need to call
 * this function twice. We only check the existence of the base to be modified. 
 * We do not check whether the base to be modified is already paired to
 * somebody else.\n
 * Position counting starts at 0.\n
 *
 * @param[in] to The base type be modified.
 * @param[in] from Pairing position.
 * @param[in] this Rna data object.
 */
void
rna_set_pairing_base (const unsigned long to, const unsigned long from,
                      Rna* this)
{
   assert (this);
   assert (this->pairs);
   assert (this->size > to);

   this->pairs[to] = from;
}

/** @brief Copy a given sequence into an Rna object.
 *
 * Set the sequence component of an Rna object. Since the input sequence is
 * copied, memory for it has to be allocated before calling. The sequence is
 * stored verbatim, this means that it is not transformed into numbers or vice
 * versa automatically. If @c len exceeds the allocated memory, only a fitting
 * prefix of the sequence is copied. We do not check whether sequence consists
 * only of symbols from a certain RNA alphabet.
 *
 * @param[in] sequence Sequence to be copied.
 * @param[in] len Length of the sequence.
 * @param[in] this Rna object.
 */
void
rna_set_sequence (const char* sequence, const unsigned long len, Rna* this)
{
   assert (this);
   assert (sequence);

   memcpy (this->seq, sequence, (len < this->size ? len : this->size));
}

/** @brief Transform Rna sequence to numbers.
 *
 * Transform the sequence component of an Rna object to numbers.\n
 * Returns 0 on success, @c ERR_RNA_NO_BASE if the sequence contains invalid
 * bases. In case of error, the sequence is unchanged.
 *
 * @param[in] sigma Alphabet to use for transformation.
 * @param[in] this Rna data object.
 */
int
rna_transform_sequence_2_no (const Alphabet* sigma, Rna* this)
{
   unsigned long i;
   char tmp;

   assert (sigma);
   assert (this);
   assert (this->seq);

   for (i = 0; i < this->size; i++)
   {
      tmp = alphabet_base_2_no (this->seq[i], sigma);

      if (tmp == CHAR_UNDEF)
      {
         if (i > 0)
         {
            for (; i > 0; i--)
            {
               this->seq[i - 1] = alphabet_no_2_base (this->seq[i - 1], sigma);
            }
         }
         return ERR_RNA_NO_BASE;
      }

      this->seq[i] = tmp;
   }

   return 0;
}

/** @brief Transform Rna sequence from numbers to bases.
 *
 * Transform the sequence component of an Rna object to bases.\n
 * Returns 0 on success, @c ERR_RNA_NO_BASE if the sequence contains invalid
 * bases. In case of error, the sequence is unchanged.
 *
 * @param[in] sigma Alphabet to use for transformation.
 * @param[in] this Rna data object.
 */
int
rna_transform_sequence_2_bases (const Alphabet* sigma, Rna* this)
{
   unsigned long i;
   char tmp;

   assert (sigma);
   assert (this);
   assert (this->seq);

   for (i = 0; i < this->size; i++)
   {
      tmp = alphabet_no_2_base (this->seq[i], sigma);

      if (tmp == CHAR_UNDEF)
      {
         if (i > 0)
         {
            for (; i > 0; i--)
            {
               this->seq[i - 1] = alphabet_base_2_no (this->seq[i - 1], sigma);
            }
         }
         return ERR_RNA_NO_BASE;
      }

      this->seq[i] = tmp;
   }

   return 0;
}


/*********************************   Access   *********************************/

/** @brief Get the no. of nucleotides of an Rna object.
 *
 * Returns the no. of nucleotides of a given Rna object.
 *
 * @param[in] this Rna data object.
 */
unsigned long
rna_get_size (const Rna* this)
{
   assert (this);

   return this->size;
}

/** @brief Get the list of pairs of an Rna object.
 *
 * Returns the pairing list of an Rna object. If no list was set before,
 * returns @c NULL.
 *
 * @params[in] this Rna data object.
 */
const unsigned long*
rna_get_pairlist (const Rna* this)
{
   assert (this);

   return this->pairs;
}

/** @brief Get the sequence of an Rna object.
 *
 * Returns the sequence component of an Rna object as a cstring. If no sequence
 * was set before, returns @c NULL.
 *
 * @params[in] this Rna data object.
 */
char*
rna_get_sequence (const Rna* this)
{
   assert (this);

   return this->seq;
}

/** @brief Get the info component of an Rna object.
 *
 * Returns the info component of an Rna object as a cstring. If no sequence
 * was set before, returns @c NULL.
 *
 * @params[in] this Rna data object.
 */
const Str*
rna_get_info (const Rna* this)
{
   assert (this);

   return this->info;
}

/** @brief Get a base from a certain position in an Rna sequence.
 *
 * Get a base from a certain position in the sequence component of an Rna
 * object. Position counting starts at 0.\n
 *
 * @param[in] pos Position where to place base.
 * @param[in] this Rna data object.
 */
char
rna_get_sequence_base (const unsigned long pos, const Rna* this)
{
   assert (this);
   assert (this->seq);
   assert (this->size > pos);

   return this->seq[pos];

}

/** @brief Get the pairing partner of a base.
 *
 * Returns the position of a pairing partner of a position in the RNA sequence.
 * Counting starts at 1, @c NOT_PAIRED is returned for unpaired positions.\n
 *
 * @params[in] pos Position to be evaluated.
 * @params[in] this Rna data object.
 */
unsigned long
rna_base_pairs_with (const unsigned long pos, const Rna* this)
{
   assert (this);
   assert (this->pairs);
   assert (this->size > pos);

   return this->pairs[pos];
}

/** @brief Check all base pairs of a structure to be valid.
 *
 * For a given secondary structure and sequence, test all base pairs to be
 * allowed. This should be used to verify if a RNA secondary structure can be
 * evaluated with a certain scoring function.\n
 * Returns the size of the RNA if all base pairs are allowed or the index of
 * the 5' base of the first base pair not allowed.
 *
 * @param[in] this RNA structure to be checked.
 */
unsigned long
rna_validate_basepairs (bool (*validate_basepair) (const char, const char,
                                                   void*),
                        void* scores, const Rna* this)
{
   unsigned long k;

   assert (validate_basepair);
   assert (this);
   assert (this->seq);
   assert (this->pairs);

   for (k = 0; k < this->size; k++)
   {
      if (this->pairs[k] != NOT_PAIRED)
      {
         if (validate_basepair (this->seq[k], this->seq[this->pairs[k]], scores)
             == false)
         {
            return k;
         }  
      }
   }

   return this->size;
}

unsigned long
rna_validate_basepairs_nn_scores (NN_scores* scores, const Rna* this)
{
   return rna_validate_basepairs (nn_scores_is_allowed_basepair, scores, this);
}

int
rna_secstruct_calculate_DG (const NN_scores* scores, const Rna* this)
{
   assert (this);

   return secstruct_calculate_DG (this->seq, scores, this->structure);
}

unsigned long
rna_secstruct_get_noof_stacks (const Rna* this)
{
   return secstruct_get_noof_stacks (this->structure);
}

/** @brief Get the 5' and 3' base of a stem base pair 
 *
 * @param[in] i     5' base container.
 * @param[in] j     3' base container.
 * @param[in] stack no. of stack.
 * @param[in] this  Rna object.
 */
void
rna_secstruct_get_i_geometry_stack (unsigned long* i, unsigned long* j,
                                const unsigned long stack,
                                const Rna* this)
{
   assert (this);

   secstruct_get_i_geometry_stack (i, j, stack, this->structure);
}

unsigned long
rna_secstruct_get_noof_hairpins (const Rna* this)
{
   return secstruct_get_noof_hairpins (this->structure);
}

void
rna_secstruct_get_geometry_hairpin (unsigned long* start,
                                unsigned long* end,
                                unsigned long* size,
                                const unsigned long i,
                                const Rna* this)
{
   assert (this);

   secstruct_get_geometry_hairpin(start, end, size, i, this->structure);
}

/** @brief get the start base of the ith hairpin loop of a 2D structure
 *
 * @param[in] i Index of loop.
 * @param[in] this Rna object.
 */
unsigned long
rna_secstruct_get_i_start_hairpin (const unsigned long i, const Rna* this)
{
   assert (this);

   return secstruct_get_i_start_hairpin (i, this->structure);
}

/** @brief get the last base of the ith hairpin loop of a 2D structure
 *
 * @param[in] i Index of loop.
 * @param[in] this Rna object.
 */
unsigned long
rna_secstruct_get_i_end_hairpin (const unsigned long i, const Rna* this)
{
   assert (this);

   return secstruct_get_i_end_hairpin (i, this->structure);
}

/** @brief get the no. of unpaired bases of the ith hairpin loop
 *
 * @param[in] i Index of loop.
 * @param[in] this Rna object.
 */
unsigned long
rna_secstruct_get_i_size_hairpin (const unsigned long i, const Rna* this)
{
   assert (this);

   return secstruct_get_i_size_hairpin (i, this->structure);
}

/** @brief Get no. of stems in an external loop.
 *
 * @param[in] this Rna object.
 */
unsigned long
rna_secstruct_get_noof_stems_extloop (const Rna* this)
{
   assert (this);

   return secstruct_get_noof_stems_extloop (this->structure);
}

/** @brief Get the 5' base of the ith stem of the external loop of a structure.
 *
 * @param[in] i Index of stem.
 * @param[in] this Rna object.
 */
unsigned long
rna_secstruct_get_i_5p_stem_extloop (const unsigned long i, const Rna* this)
{
   assert (this);

   return secstruct_get_i_5p_stem_extloop (i, this->structure);
}

/** @brief Get the 3' base of the ith stem of the external loop of a structure.
 *
 * @param[in] i Index of stem.
 * @param[in] this Rna object.
 */
unsigned long
rna_secstruct_get_i_3p_stem_extloop (const unsigned long i, const Rna* this)
{
   assert (this);

   return secstruct_get_i_3p_stem_extloop (i, this->structure);
}

/** @brief Get the 5' and 3' base position of the ith stem in an external loop.
 *
 * @param[out] p5 5' psoition.
 * @param[out] p3 3' position.
 * @param[in]  i    Stem.
 * @param[in]  this Rna object.
 */
void
rna_secstruct_get_i_stem_extloop (unsigned long* p5,
                                  unsigned long* p3,
                                  const unsigned long i,
                                  const Rna* this)
{
   assert (this);

   secstruct_get_i_stem_extloop (p5, p3, i, this->structure);
}

/** @brief Get the number of 5' dangles in an external loop.
 *
 * @param[in] this Rna object.
 */
unsigned long
rna_secstruct_get_noof_5pdangles_extloop (const Rna* this)
{
   assert (this);

   return secstruct_get_noof_5pdangles_extloop (this->structure);
}

/** @brief Get the number of 3' dangles in an external loop.
 *
 * @param[in] this Rna object.
 */
unsigned long
rna_secstruct_get_noof_3pdangles_extloop (const Rna* this)
{
   assert (this);

   return secstruct_get_noof_3pdangles_extloop (this->structure);
}

/** @brief Get the ith 5', 3' and free base of a 5' dangling end.
 *
 * @param[out] p5 Container for the 5' base position.
 * @param[out] p3 Container for the 3' base position.
 * @param[out] fb Container for the free base.
 * @param[in]  i  ith 5' dangling end.
 * @param[in]  this Rna object.
 */
void
rna_secstruct_get_i_5pdangle_extloop (unsigned long* p5,
                                      unsigned long* p3,
                                      unsigned long* fb,
                                      const unsigned long i,
                                      const Rna* this)
{
   assert (this);

   secstruct_get_i_5pdangle_extloop (p5,  p3, fb, i, this->structure);
}

/** @brief Get the ith 5', 3' and free base of a 3' dangling end.
 *
 * @param[out] p5 Container for the 5' base position.
 * @param[out] p3 Container for the 3' base position.
 * @param[out] fb Container for the free base.
 * @param[in]  i  ith 3' dangling end.
 * @param[in]  this Rna object.
 */
void
rna_secstruct_get_i_3pdangle_extloop (unsigned long* p5,
                                      unsigned long* p3,
                                      unsigned long* fb,
                                      const unsigned long i,
                                      const Rna* this)
{
   assert (this);

   secstruct_get_i_3pdangle_extloop (p5,  p3, fb, i, this->structure);
}

/** @brief Get no. of multiloops in a secondary structure 
 *
 * @param[in] this Rna object.
 */
unsigned long
rna_secstruct_get_noof_multiloops (const Rna* this)
{
   assert (this);

   return secstruct_get_noof_multiloops (this->structure);
}

/** @brief Get no. of stems of a certain multiloop 
 *
 * @param[in] i index of multiloop.
 * @param[in] this Rna object.
 */
unsigned long
rna_secstruct_get_i_noof_stems_multiloop (unsigned long i, const Rna* this)
{
   assert (this);

   return secstruct_get_i_noof_stems_multiloop (i, this->structure);
}

/** @brief Get the 5' and 3' base position of a stem in an multiloop.
 *
 * @param[out] p5 5' psoition.
 * @param[out] p3 3' position.
 * @param[in]  i    Stem.
 * @param[in]  j    Multiloop.
 * @param[in]  this Rna object.
 */
void
rna_secstruct_get_i_stem_multiloop (unsigned long* p5, unsigned long* p3,
                                    const unsigned long i,
                                    const unsigned long j,
                                    const Rna* this)
{
   assert (this);

   secstruct_get_i_stem_multiloop (p5, p3, i, j, this->structure);
}

/** @brief Get the no. of 5' dangling ends involved in a certain multiloop. 
 *
 * @param[in] i Index of multiloop.
 * @param[in] this Rna object.
 */
unsigned long
rna_secstruct_get_i_noof_5pdangles_multiloop (unsigned long i, const Rna* this)
{
   assert (this);

   return secstruct_get_i_noof_5pdangles_multiloop (i, this->structure);
}

/** @brief Get the ith 5', 3' and free base of a 5' dangling end.
 *
 * @param[out] p5 Container for the 5' base position.
 * @param[out] p3 Container for the 3' base position.
 * @param[out] fb Container for the free base.
 * @param[in]  i  ith 5' dangling end.
 * @param[in]  j  jth multiloop.
 * @param[in]  this Rna object.
 */
void
rna_secstruct_get_i_5pdangle_multiloop (unsigned long* p5,
                                        unsigned long* p3,
                                        unsigned long* fb,
                                        const unsigned long i,
                                        const unsigned long j,
                                        const Rna* this)
{
   assert (this);

   secstruct_get_i_5pdangle_multiloop (p5, p3, fb, i, j, this->structure);
}

/** @brief Get the no. of 3' dangling ends involved in a certain multiloop. 
 *
 * @param[in] i Index of multiloop.
 * @param[in] this Rna object.
 */
unsigned long
rna_secstruct_get_i_noof_3pdangles_multiloop (unsigned long i,
                                              const Rna* this)
{
   assert (this);
   
   return secstruct_get_i_noof_3pdangles_multiloop (i, this->structure);
}

/** @brief Get the ith 5', 3' and free base of a 3' dangling end.
 *
 * @param[out] p5 Container for the 5' base position.
 * @param[out] p3 Container for the 3' base position.
 * @param[out] fb Container for the free base.
 * @param[in]  i  ith 3' dangling end.
 * @param[in]  j  jth multiloop.
 * @param[in]  this Secondary structure.
 */
void
rna_secstruct_get_i_3pdangle_multiloop (unsigned long* p5,
                                        unsigned long* p3,
                                        unsigned long* fb,
                                        const unsigned long i,
                                        const unsigned long j,
                                        const Rna* this)
{
   assert (this);

   secstruct_get_i_3pdangle_multiloop (p5, p3, fb, i, j, this->structure);
}

/** @brief get the no. of bulge loops of a 2D structure
 *
 * @param[in] this Rna object.
 */
unsigned long
rna_secstruct_get_noof_bulges (const Rna* this)
{
   assert (this);

   return secstruct_get_noof_bulges (this->structure);
}

/** @brief Retrieve bulge loop geometry in one go.
 *
 * @param[in/out] i1   Base i of opening pair.
 * @param[in/out] j1   Base j of opening pair.
 * @param[in/out] i2   Base i of closing pair.
 * @param[in/out] j2   Base j of closing pair
 * @param[in/out] size Size of the loop.
 * @param[in]     i    No. of bulge.
 * @param[in]     this Rna object.
 */
void
rna_secstruct_get_geometry_bulge (unsigned long* i1,
                                  unsigned long* j1,
                                  unsigned long* i2,
                                  unsigned long* j2,
                                  unsigned long* size,
                                  const unsigned long i,
                                  const Rna* this)
{
   assert (this);

   secstruct_get_geometry_bulge (i1, j1, i2, j2, size, i, this->structure);
}

/** @brief get the no. of internal loops of a 2D structure
 *
 * @param[in] this Rna object.
 */
unsigned long
rna_secstruct_get_noof_internals (const Rna* this)
{
   assert (this);

   return secstruct_get_noof_internals (this->structure);
}

/** @brief get the geometry of a certain internal loop.
 *
 * @param[out] i1    Base i of opening pair.
 * @param[out] j1    Base j of opening pair.
 * @param[out] i2    Base i of closing pair.
 * @param[out] j2    Base j of closing pair.
 * @param[out] size1 Size of first loop.
 * @param[out] size2 Size of 2nd loop.
 * @param[in] i      Index of loop.
 * @param[in] this   Rna object.
 */
void
rna_secstruct_get_geometry_internal (unsigned long* i1,
                                     unsigned long* j1,
                                     unsigned long* i2,
                                     unsigned long* j2,
                                     unsigned long* size1,
                                     unsigned long* size2,
                                     const unsigned long i,
                                     const Rna* this)
{
   assert (this);

   secstruct_get_geometry_internal (i1, j1, i2, j2, size1, size2, i,
                                    this->structure);
}
