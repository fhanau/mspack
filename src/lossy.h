#ifndef MSPREDICT_LOSSY_H
#define MSPREDICT_LOSSY_H

#include <tuple>
#include <cmath>
#include <limits>
#include <stdexcept>

#define MZ_DEFAULT_ERROR 0.0001
#define INT_DEFAULT_ERROR 0.01

//MAX_LOSS is the maximum absolute loss.
#define MAX_LOSS 1.0
#define MAX_REL_LOSS 0.5
#define LOSSY_FIXED_PO2_CUTOFF (1.0 / 65536.0)
#define LOSSY_MZ_CUTOFF (1.0 / 1024.0 / 1024.0 / 1024.0)
#define LOSSY_LOG_CUTOFF (1.0 / 32768.0)
#define LOSSY_REL_CUTOFF_SINGLE (1.0 / 1024.0 / 1024.0 / 8.0)
#define LOSSY_REL_CUTOFF_DOUBLE (1.0 / 1024.0 / 1024.0 / 1024.0)

typedef enum LossyMode {
  lossless = 0,
  truncate_to_int = 1,
  truncate_bits_ae = 2,
  truncate_bits_re = 3,
  fixed_ae = 4,
  log_re = 5
} LossyMode;

struct MSOptions {
  LossyMode mz_lossy_mode;
  double mz_lossy_error;
  LossyMode int_lossy_mode;
  double int_lossy_error;
  int no_buckets;
  unsigned block_size;
  //For debug purposes
  int scans_only;
} __attribute__((packed));

template<typename F, typename U> std::tuple<F, U> apply_loss(U val, LossyMode mode, double lossy_error) {
  F val_fp = reinterpret_cast<F&>(val);

  if(!std::isfinite(val_fp) || val_fp < 0.0){
    //We can not guarantee that this will not overlap with a value from a transform, throw an error for now.
    throw std::out_of_range("value out of range");
  }
  if(val_fp == 0.0){
    return std::make_tuple(val_fp, val);
  }
  if(mode == lossless){
    return std::make_tuple(val_fp, val);
  }
  if(mode == truncate_to_int){
    //If the value is beyond the range where the fp format  has an accuracy of at least one, it still works:
    //We cast to an int that is the same value (no additional loss) and cast back to the fp value we had before.
    if(val_fp > std::numeric_limits<U>::max()){
      throw std::out_of_range("value too large");
    }
    val = (U)floor(val_fp + 0.5);
    val_fp = (F)val;
    return std::make_tuple(val_fp, val);
  }
  if(mode == truncate_bits_re || mode == truncate_bits_ae){
    int frac_bits = sizeof(U) == 8 ? 52 : 23;
    int bits_to_trunc = frac_bits + ilogb(lossy_error) + 1;
    if(mode == truncate_bits_ae){
      bits_to_trunc -= ilogb(val_fp);
      if(bits_to_trunc > frac_bits) {
          bits_to_trunc = frac_bits;
      }
    }
    if(bits_to_trunc > 0){
      U mask = (((U)1) << bits_to_trunc) - 1;
      //This allows us to save an extra bit.
      U mask2 = ((U)1) << (bits_to_trunc - 1);
      if((val & mask) >= mask2){
        val = val & (~mask);
        val += (((U)1) << bits_to_trunc);
      }
      else {
        val = val & (~mask);
      }
      val_fp = reinterpret_cast<F&>(val);
    }
    return std::make_tuple(val_fp, val);
  }
  if (mode == fixed_ae) {
    double scaling_factor = 0.5 / lossy_error;
    if (val_fp * scaling_factor > std::numeric_limits<U>::max()) {
      throw std::out_of_range("value too large");
    }
    val = (U)floor(((double)val_fp * scaling_factor) + 0.5);
    val_fp = val / scaling_factor;
    return std::make_tuple(val_fp, val);
  }
  if(mode == log_re){
    if(val_fp < 1.0){
      //Should be exceedingly rare in practice, so it won't really hurt the compression ratio
      int num_bits = sizeof(U) == 8 ? 64 : 32;

      val |= (((U)1) << (num_bits - 1));
      return std::make_tuple(val_fp, val);
    }
    double base = 1.0 + lossy_error * 2.0;
    double k = 1.0 + lossy_error;
    double log_k = 1.0 - (log(k) / log(base));
    U idx = (U)floor(log((double)val_fp) / log(base) + log_k);
    val_fp = pow(base, idx);
    return std::make_tuple(val_fp, idx + 1);
  }
  throw std::runtime_error("unknown mode");
}

template<typename F, typename U> F apply_lossF(U val, LossyMode mode, double lossy_error){
  F f;
  U u;
  std::tie(f, u) = apply_loss<F,U>(val, mode, lossy_error);
  return f;
}

template<typename F, typename U> U apply_lossU(U val, LossyMode mode, double lossy_error){
  F f;
  U u;
  std::tie(f, u) = apply_loss<F,U>(val, mode, lossy_error);
  return u;
}

template<typename F, typename U> U decode_lossy(U val, LossyMode mode, double lossy_error) {
  F val_fp = reinterpret_cast<F&>(val);
  int num_bits = sizeof(U) == 8 ? 64 : 32;

  if(mode == log_re && (val >> (num_bits - 1)) > 0){
    return val & ((((U)1) << (num_bits - 1)) - 1);
  }
  if(mode == lossless){
    return val;
  }
  if(mode == truncate_to_int){
    val_fp = (F)val;
    return reinterpret_cast<U&>(val_fp);
  }
  if(mode == truncate_bits_re || mode == truncate_bits_ae){
    return val;
  }
  if(mode == fixed_ae){
    double scaling_factor = 0.5/lossy_error;
    val_fp = val / scaling_factor;
    return reinterpret_cast<U&>(val_fp);
  }
  if(mode == log_re){
    if(val == 0){
      return val;
    }
    double base = 1.0 + lossy_error * 2.0;
    val_fp = pow(base, val - 1);
    return reinterpret_cast<U&>(val_fp);
  }
  throw std::runtime_error("unknown mode");
}

#endif
