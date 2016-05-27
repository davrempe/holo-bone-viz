#ifndef ____RGL___DEFS___H___
#define ____RGL___DEFS___H___

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef __cplusplus 
#define BEGIN_C_DECL extern "C" {
#define END_C_DECL }
#else
#define BEGIN_C_DECL 
#define END_C_DECL 
#endif

BEGIN_C_DECL

/**
 * A more descriptive type for a general pointer.
 */

typedef void * ptr_t;

/**
 * A more descriptive type for a byte (8 bits).
 */

typedef unsigned char byte_t;

/**
 * A more descriptive type for a word (32 bits).
 */

typedef unsigned int word_t;

/**
 * A more descriptive type for a string.
 */

typedef char * string_t;

/**
 * A more descriptive type for a chunk of bytes.
 */

typedef byte_t * chunk_t;

/**
 * A more descriptive type for a boolean type.
 */

typedef int boolean_t_new;

#define TRUE 1
#define FALSE 0

#define NIL 0x0

END_C_DECL

#endif



