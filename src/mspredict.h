#ifndef MSPREDICT_H
#define MSPREDICT_H

#include <vector>
#include "utils.h"
#include "lossy.h"

static void write_to_vector4(uint32_t val, std::vector<std::vector<unsigned char> >& out_vec, unsigned offset){
  for(int k = 0; k < 4; k++){
    uint8_t c = val & 255;
    val >>= 8;
    out_vec[k + offset].push_back(c);
  }
}

static void write_to_vector8(uint64_t val, std::vector<std::vector<unsigned char> >& out_vec, unsigned offset){
  for(int k = 0; k < 8; k++){
    uint8_t c = val & 255;
    val >>= 8;
    out_vec[k + offset].push_back(c);
  }
}

static void write_to_vector(uint64_t val, std::vector<std::vector<unsigned char> >& out_vec, unsigned offset, unsigned size){
  if(size == 4){
    write_to_vector4(val, out_vec, offset);
    return;
  }
  write_to_vector8(val, out_vec, offset);
}

template <typename FA, typename FM, typename U>
void process_generic_bucket(FA** otherdata, FM** bucketData, size_t bucket_elements_scan[],
  size_t bucket_elements, std::vector<std::vector<unsigned char> >& main_output, uint64_t& all_integer_elements,
  bool is_mz, unsigned offset, LossyMode lossy_mode, double lossy_error){
  FA max_element = 1.0;
  all_integer_elements = 0;

  for (size_t i = 0; i < bucket_elements; i++) {
    for (size_t j = 0; j < bucket_elements_scan[i]; j++) {
      if(otherdata){
        max_element = std::max(max_element, otherdata[i][j]);
      }
      else {
        max_element = std::max(max_element, (FA)bucketData[i][j]);
      }
    }
  }

  //Not necessary for intensities since we have the adjusted stream.
  if(is_mz){
    U max_element2 = reinterpret_cast<U&>(max_element);
    std::tie(max_element, max_element2) = apply_loss<FM, U>(max_element2, lossy_mode, lossy_error);
  }

  uint64_t bucket_cnt = bucket_idx<FA>(max_element) + 1;
  //printf("Buckets: %llu\n", bucket_cnt);

  std::vector<std::vector<U> > buckets(bucket_cnt);
  for (size_t i = 0; i < bucket_elements; i++) {
    FM prevF = 0.0;
    U prevI = 0;
    for (size_t j = 0; j < bucket_elements_scan[i]; j++) {
      U val = *(U*)&bucketData[i][j];
      unsigned idx;
      FM lossyVal;
      std::tie(lossyVal, val) = apply_loss<FM, U>(val, lossy_mode, lossy_error);
      if(!otherdata){
        idx = bucket_idx<FM>(prevF);
      }
      if(otherdata){
        if(j == bucket_elements_scan[i] - 1){
          idx = 0;
        }
        else {
          idx = bucket_idx<FA>(otherdata[i][j]);
        }
      }
      if(!is_mz){
        buckets[idx].push_back(val);
      }
      else {
        buckets[idx].push_back(val - prevI);
      }
      prevI = val;
      prevF = lossyVal;
    }
  }

  if(!otherdata){
    for(uint64_t i = 0; i < bucket_cnt; i++){
      uint64_t val = buckets[i].size();
      write_to_vector8(val, main_output, offset);
    }
    offset += 8;
  }
#ifdef DEBUG_PRINT
  fprintf(stderr, "Buckets: %lu\n", bucket_cnt);
#endif

  //Do not use vector push_back directly to improve performance and memory usage.
  size_t total_size = 0;
  for(uint64_t i = 0; i < bucket_cnt; i++){
    total_size += buckets[i].size();
  }
  for (size_t j = 0; j < sizeof(U); j++) {
    main_output[offset + j].resize(main_output[offset + j].size() + total_size);
  }
  size_t pos = 0;
  for(uint64_t i = 0; i < bucket_cnt; i++){
    uint64_t sz = buckets[i].size();

    for(uint64_t j = 0; j < sz; j++){
      U val = buckets[i][j];
      for(size_t k = 0; k < sizeof(U); k++){
        uint8_t c = val & 255;
        val >>= 8;
        main_output[k + offset][pos] = c;
      }
      pos++;
    }
  }
}

template <typename F, typename U>
void process_generic_seq(F** bucketData, size_t bucket_elements_scan[],
  size_t bucket_elements, std::vector<std::vector<unsigned char> >& main_output, uint64_t& all_integer_elements,
  bool is_mz, unsigned offset, LossyMode lossy_mode, double lossy_error){
  all_integer_elements = 0;

  for (size_t i = 0; i < bucket_elements; i++) {
    U prev_val = 0;
    for (size_t j = 0; j < bucket_elements_scan[i]; j++) {
      U val;
      memcpy(&val, &bucketData[i][j], sizeof(F));
      val = apply_lossU<F, U>(val, lossy_mode, lossy_error);
      if(is_mz) {
        U new_val = val - prev_val;
        prev_val = val;
        val = new_val;
      }
      write_to_vector(val, main_output, offset, sizeof(U));
    }
  }
}

//Mid-level API
template <typename F, typename U>
void process_mz_bucket(F** mzBucketData, size_t bucket_elements_scan[],
  size_t bucket_elements, std::vector<std::vector<unsigned char> >& main_output, uint64_t& all_integer_elements,
  unsigned offset, LossyMode lossy_mode, double lossy_error){
    process_generic_bucket<F, F, U>(0, mzBucketData, bucket_elements_scan,
      bucket_elements, main_output, all_integer_elements,
      true, offset, lossy_mode, lossy_error);
}

template <typename FA, typename FM, typename U>
void process_int_bucket(FA** mzBucketData, FM** intBucketData, size_t bucket_elements_scan[],
  size_t bucket_elements, std::vector<std::vector<unsigned char> >& main_output, uint64_t& all_integer_elements,
  unsigned offset, LossyMode lossy_mode, double lossy_error){
    process_generic_bucket<FA, FM, U>(mzBucketData, intBucketData, bucket_elements_scan,
      bucket_elements, main_output, all_integer_elements,
      false, offset, lossy_mode, lossy_error);
}

template <typename F, typename U>
void mspredict_mz_seq(F** bucketData, size_t* bucket_elements_scan,
  size_t bucket_elements, std::vector<std::vector<unsigned char> >& main_output, uint64_t& all_integer_elements,
  unsigned offset, LossyMode lossy_mode, double lossy_error) {
  process_generic_seq<F, U>(bucketData, bucket_elements_scan,
  bucket_elements, main_output, all_integer_elements,
  true, offset, lossy_mode, lossy_error);
}

template <typename F, typename U>
void mspredict_int_seq(F** bucketData, size_t bucket_elements_scan[],
  size_t bucket_elements, std::vector<std::vector<unsigned char> >& main_output, uint64_t& all_integer_elements,
  unsigned offset, LossyMode lossy_mode, double lossy_error) {
  process_generic_seq<F, U>(bucketData, bucket_elements_scan,
  bucket_elements, main_output, all_integer_elements,
  false, offset, lossy_mode, lossy_error);
}

//High-level API

//Move options struct and LossyMode struct here

template <typename F, typename U, typename F2, typename U2>
int mspredict_both_bucket(F** mzBucketData, F2** intBucketData, size_t* bucket_elements_scan,
  size_t bucket_elements, unsigned char** out, size_t* out_size,
  LossyMode mz_lossy_mode, double mz_lossy_error, LossyMode int_lossy_mode, double int_lossy_error) {
  std::vector<std::vector<unsigned char>> main_output(8 + sizeof(F) + sizeof(F2));
  uint64_t all_integer_elements = 0;
  process_mz_bucket<F, U>(mzBucketData, bucket_elements_scan,
  bucket_elements, main_output, all_integer_elements,
  true, 0, mz_lossy_mode, mz_lossy_error);
  //TODO: If not lossless, do create a secondary array of mz actual value after lossy
  process_generic_bucket<F, F2, U2>(mzBucketData, intBucketData, bucket_elements_scan,
  bucket_elements, main_output, all_integer_elements,
  false, 8 + sizeof(F), int_lossy_mode, int_lossy_error);
  size_t num_elem = main_output[0].size();
  *out_size = num_elem * sizeof(F);

  //TODO: Include num_buckets, num_bucket_elem, all_integer_elements in output
  *out = (unsigned char*)malloc((*out_size) + 1);
  if (!(*out)) {
    return 1;
  }
  (*out)[0] = (all_integer_elements != 0);
  size_t pos = 1;
  for(int i = 0; i < 8 + sizeof(F) + sizeof(F2); i++) {
    memcpy(&((*out)[pos]), main_output[i].data(), num_elem);
    pos += num_elem;
  }
  return 0;
}

template <typename F, typename U>
int mspredict_generic_seq(F** bucketData, size_t* bucket_elements_scan,
  size_t bucket_elements, unsigned char is_mz, unsigned char** out, size_t* out_size,
  LossyMode lossy_mode, double lossy_error) {
  std::vector<std::vector<unsigned char>> main_output(sizeof(F));
  uint64_t all_integer_elements = 0;
  process_generic_seq<F, U>(bucketData, bucket_elements_scan,
  bucket_elements, main_output, all_integer_elements,
  is_mz, 0, lossy_mode, lossy_error);
  size_t num_elem = main_output[0].size();
  *out_size = num_elem * sizeof(F);

  //To keep
  *out = (unsigned char*)malloc((*out_size) + 1);
  if (!(*out)) {
    return 1;
  }
  //TODO: Also include number of elements for decoding
  (*out)[0] = (all_integer_elements != 0);
  size_t pos = 1;
  for(int i = 0; i < sizeof(F); i++) {
    memcpy(&((*out)[pos]), main_output[i].data(), num_elem);
    pos += num_elem;
  }
  return 0;
}

template <typename F, typename U>
int mspredict_mz_seq(F** bucketData, size_t* bucket_elements_scan,
  size_t bucket_elements, unsigned char** out, size_t* out_size,
  LossyMode lossy_mode, double lossy_error) {
  return mspredict_generic_seq<F, U>(bucketData, bucket_elements_scan,
  bucket_elements, out, out_size, true, lossy_mode, lossy_error);
}

template <typename F, typename U>
int mspredict_int_seq(F** bucketData, size_t* bucket_elements_scan,
  size_t bucket_elements, unsigned char** out, size_t* out_size,
  LossyMode lossy_mode, double lossy_error) {
  return mspredict_generic_seq<F, U>(bucketData, bucket_elements_scan,
    bucket_elements, out, out_size, false, lossy_mode, lossy_error);
}

#endif
