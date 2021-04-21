/*******************************************************************************************
 *
 *  Non-renetrent general sequence format reader
 *
 *  Author:  Gene Myers
 *  Date  :  April, 2021
 *
 ********************************************************************************************/

#ifndef _READER
#define _READER

#include <stdint.h>
#include <zlib.h>

#include "DB.h"

typedef struct
  { char   *seq;     // DNA sequence
    uint64  len;     // length of DNA sequence
  } Entry;

typedef Entry *(*Fetcher)();

Fetcher Open_File(char *path);

#endif // _READER
