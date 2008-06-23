/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbapps/fold.h
 *
 *  @brief RNA Secondary Structure Prediction Tool
 *
 *  Module: fold
 *
 *  Library: libcrbapps
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-06-17
 *
 *
 *  Revision History:
 *         - 2008Jun17 bienert: created
 *
 */


#ifdef __cplusplus
extern "C" {
#endif

#ifndef FOLD_H
#define FOLD_H

int
fold_main(const char*);
   
#endif /* FOLD_H */

#ifdef __cplusplus
}
#endif
