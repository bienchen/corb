/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbrna/nn_scores.h
 *
 *  @brief Neares neighbour Model for evaluating RNA secondary structures
 *
 *  Module: nn_scores
 *
 *  Library: crbrna
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-07-03
 *
 *
 *  Revision History:
 *         - 2008Jul03 bienert: created
 *
 */


#ifdef __cplusplus
extern "C" {
#endif

#ifndef NN_SCORES_H
#define NN_SCORES_H


typedef struct NN_scores NN_scores;
   

/**********************   Constructors and destructors   **********************/

NN_scores*
nn_scores_new (const char*, const int);

#define NN_SCORES_NEW nn_scores_new (__FILE__, __LINE__)

NN_scores*
nn_socres_new_init (Alphabet*, const char*, const int);

#define NN_SCORES_NEW_INIT(A) nn_socres_new_init (A, __FILE__, __LINE__)

void
nn_scores_delete (NN_scores*);


/*********************************   Access   *********************************/

void
nn_scores_get_allowed_basepair (unsigned, char*, char*, const NN_scores*);


/*********************************    Size    *********************************/

unsigned long
nn_scores_no_allowed_basepairs (const NN_scores*);

/*********************************   Output   *********************************/

void
nn_scores_fprintf_bp_allowed (FILE*, const NN_scores*, const Alphabet*);

void
nn_scores_fprintf_bp_idx (FILE*, const NN_scores*, const Alphabet*);

void
nn_scores_fprintf_G_stack (FILE*, const NN_scores*, const Alphabet*);

void
nn_scores_fprintf_mm_G_stack (FILE*, const NN_scores*, const Alphabet*);


/******************************   Miscellaneous   *****************************/

unsigned long
nn_scores_bp_2_idx (const char, const char, const NN_scores*);


#endif /* NN_SCORES_H */

#ifdef __cplusplus
}
#endif
