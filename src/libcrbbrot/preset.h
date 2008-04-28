/*
 * Copyright (C) 2008 Stefan Bienert
 * Copyright (C) 2008 Center for Bioinformatics, University of Hamburg
 *
 * See COPYING file in the top level directory of this tree for licence.
 */

/*
 ****   Documentation header   ***
 *
 *  @file libcrbbrot/preset.h
 *
 *  @brief Storing fixed sites for a sequence matrix
 *
 *  Module: preset
 *
 *  Library: libcrbbrot
 *
 *  Project: CoRB - Collection of RNAanalysis Binaries
 *
 *  @author Stefan Bienert
 *
 *  @date 2008-03-04
 *
 *
 *  Revision History:
 *         - 2008Mar04 bienert: created
 *
 */


#ifdef __cplusplus
extern "C" {
#endif

#ifndef PRESET_H
#define PRESET_H

enum preset_retvals{
   ERR_PA_ALLOC = 1,      /* (re)allocation problems */
};

typedef struct Preset Preset;

typedef struct PresetArray PresetArray;   

/**********************   Constructors and destructors   **********************/

Preset*
preset_new (const char*, const int);

#define PRESET_NEW preset_new (__FILE__, __LINE__)

Preset*
preset_new_preset (const char, const unsigned long, const char*, const int);

#define PRESET_NEW_PRESET(BASE, POS) preset_new_preset (BASE,     \
                                                        POS,      \
                                                        __FILE__, \
                                                        __LINE__)

PresetArray*
presetarray_new (const char*, const int);

#define PRESETARRAY_NEW presetarray_new (__FILE__, __LINE__)

PresetArray*
presetarray_new_size (const unsigned long, const char*, const int);

#define PRESETARRAY_NEW_SIZE(SIZE) presetarray_new_size (SIZE,     \
                                                         __FILE__, \
                                                         __LINE__)

void
preset_delete (Preset*);

void
presetarray_delete (PresetArray*);


/********************************   Altering   ********************************/

int
presetarray_add (const char, const unsigned long, PresetArray*);


/*********************************   Access   *********************************/

char
presetarray_get_ith_base (const unsigned long, const PresetArray*);

unsigned long
presetarray_get_ith_pos (const unsigned long, const PresetArray*);


/**********************************   Size   **********************************/

unsigned long
presetarray_get_length (const PresetArray*);

#endif /* PRESET_H */

#ifdef __cplusplus
}
#endif
