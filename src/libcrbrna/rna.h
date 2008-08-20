/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbrna/rna.h
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

#ifdef __cplusplus
extern "C" {
#endif

#ifndef RNA_H
#define RNA_H


/* unpaired nucleotide */
#define NOT_PAIRED 0
#define VIENNA_OPEN  '('
#define VIENNA_CLOSE ')'
#define VIENNA_UNPAIRED '.'

/* error numbers */
enum rna_retvals{
   ERR_RNA_ALLOC = 1,      /* (re)allocation problems */
   ERR_RNA_VIENNA_FORMAT,  /* Wrong format of Vienna structure string */
   ERR_RNA_VIENNA_MMC,     /* Mismatched closing pairing partner found */
   ERR_RNA_VIENNA_MMO     /* Mismatched opening pairing partner found*/
};

typedef struct Rna Rna;


/**********************   Constructors and destructors   **********************/

Rna*
rna_new (const char*, const int);

#define RNA_NEW rna_new (__FILE__, __LINE__)

void
rna_delete (Rna*);
   
int
rna_init_pairlist_vienna (const char*, const unsigned long, Rna*,
                          const char*, const int);

#define RNA_INIT_PAIRLIST_VIENNA(A, B, C) \
   rna_init_pairlist_vienna (A, B, C, __FILE__, __LINE__);

int
rna_allocate_sequence (const unsigned long, Rna*, const char*, const int);

#define RNA_ALLOCATE_SEQUENCE(A, B) \
   rna_allocate_sequence (A, B, __FILE__, __LINE__)

/********************************   Altering   ********************************/


/*********************************   Access   *********************************/

unsigned long
rna_base_pairs_with (unsigned long, Rna*);


#endif /* RNA_H */

#ifdef __cplusplus
}
#endif
