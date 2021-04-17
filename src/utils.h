#ifndef MSPREDICT_UTIL_H
#define MSPREDICT_UTIL_H

#include<algorithm>
#include "head.h"

#define FRAC_USE_BITS 8
#define FRAC_USE_MASK ((1 << FRAC_USE_BITS) - 1)

template <typename F>
uint32_t bucket_idx(F prev_float) {
  if (sizeof(F) == 8) {
    uint32_t frac_bits = 52;
    uint64_t val = *(uint64_t*)&prev_float;
    int64_t a = val >> frac_bits;
    a -= 1027;
    if (a < 0) {
      return 0;
    }
    a <<= FRAC_USE_BITS;
    int64_t b = (val >> (frac_bits - FRAC_USE_BITS)) & FRAC_USE_MASK; //top x bits from fraction
    return a + b;
  }
  else {
    uint32_t frac_bits = 23;
    uint32_t val = *(uint32_t*)&prev_float;
    int a = val >> frac_bits;
    a -= 131;
    if (a < 0) {
      return 0;
    }
    a <<= FRAC_USE_BITS;
    int b = (val >> (frac_bits - FRAC_USE_BITS)) & FRAC_USE_MASK; //top x bits from fraction
    return a + b;
  }
}

#ifdef __GNUC__
#define mspredict_ntohll __builtin_bswap64
#else
uint64_t mspredict_ntohll (uint64_t* in64) {
  uint8_t* in = (uint8_t*)in64;
  uint64_t out = 0;
  out += in[0]; out <<= 8;
  out += in[1]; out <<= 8;
  out += in[2]; out <<= 8;
  out += in[3]; out <<= 8;
  out += in[4]; out <<= 8;
  out += in[5]; out <<= 8;
  out += in[6]; out <<= 8;
  out += in[7];
  return out;
}

#endif
#endif

