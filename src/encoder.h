#ifndef MSPACK_ENCODER_H
#define MSPACK_ENCODER_H
#include "lossy.h"

int MSEncode(const char* input, const char* output, MSOptions& options, unsigned char is_gz);
int MZMLEncode(const char* input, const char* output, MSOptions& options, unsigned char is_gz);
#endif
