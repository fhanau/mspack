#ifndef MSPACK_H
#define MSPACK_H

#include <cstdlib>
#include <vector>

//TODO: Verify using this interface works. It should work if linked with the
//encoding code since it defines all allowed templates of the functions below.

template <typename F, typename U>
void process_mz_bucket(F** mzBucketData, size_t bucket_elements_scan[],
  size_t bucket_elements, std::vector<std::vector<unsigned char> >& main_output, uint64_t& all_integer_elements,
  unsigned offset, LossyMode lossy_mode, double lossy_error);

template <typename FA, typename FM, typename U>
void process_int_bucket(FA** mzBucketData, FM** intBucketData, size_t bucket_elements_scan[],
  size_t bucket_elements, std::vector<std::vector<unsigned char> >& main_output, uint64_t& all_integer_elements,
  unsigned offset, LossyMode lossy_mode, double lossy_error);

template <typename F, typename U>
void process_generic_seq(F** bucketData, size_t bucket_elements_scan[],
  size_t bucket_elements, std::vector<std::vector<unsigned char> >& main_output, uint64_t& all_integer_elements,
  bool is_mz, unsigned offset, LossyMode lossy_mode, double lossy_error);

//TODO: This should work since
template <typename FA, typename FM, typename UA, typename UM>
void process_ms_data(FA** mzBucketData, FM** intBucketData, size_t bucket_elements_scan[],
  size_t bucket_elements,
  FA** mzSeqData, FM** intSeqData, size_t seq_elements_scan[], size_t seq_elements,
  std::vector<std::vector<unsigned char> >& main_output, uint64_t& all_integer_elements1, uint64_t& all_integer_elements2,
  unsigned offsets[4], MSOptions options){
    uint64_t all_integer_elements = 0;
    process_mz_bucket<FA, UA>(mzBucketData, bucket_elements_scan,
      bucket_elements, main_output, all_integer_elements,
      offsets[0], options.mz_lossy_mode, options.mz_lossy_error);

    FA** lossyMzBucketData = (FA**)malloc(bucket_elements * sizeof(FA*));
    for(size_t i = 0; i < bucket_elements; i++){
      lossyMzBucketData[i] = (FA*)malloc(bucket_elements_scan[i] * sizeof(FA));
      for (size_t j = 0; j < bucket_elements_scan[i]; j++) {
        UA val = reinterpret_cast<UA&>(mzBucketData[i][j]);
        lossyMzBucketData[i][j] = apply_lossF<FA, UA>(val, options.mz_lossy_mode, options.mz_lossy_error);
      }
    }

    process_int_bucket<FA, FM, UM>(lossyMzBucketData, intBucketData, bucket_elements_scan,
      bucket_elements, main_output, all_integer_elements1,
      offsets[1], options.int_lossy_mode, options.int_lossy_error);

    for(size_t i = 0; i < bucket_elements; i++){
      free(lossyMzBucketData[i]);
    }
    free(lossyMzBucketData);

    process_generic_seq<FA, UA>(mzSeqData, seq_elements_scan,
      seq_elements, main_output, all_integer_elements,
      true, offsets[2], options.mz_lossy_mode, options.mz_lossy_error);
    process_generic_seq<FM, UM>(intSeqData, seq_elements_scan,
      seq_elements, main_output, all_integer_elements2,
      false, offsets[3], options.int_lossy_mode, options.int_lossy_error);
}

int process_ms_data(void** mzBucketData, void** intBucketData, size_t bucket_elements_scan[],
  size_t bucket_elements,
  void** mzSeqData, void** intSeqData, size_t seq_elements_scan[], size_t seq_elements,
  unsigned char** out_buf, size_t* out_size, uint64_t& all_integer_elements1, uint64_t& all_integer_elements2,
  MSOptions options, unsigned mz_depth, unsigned int_depth){
  std::vector<std::vector<unsigned char> > main_output(8 + (mz_depth + int_depth) * 2);
  unsigned offsets[4] = {0};
  offsets[1] = mz_depth;
  offsets[2] = mz_depth + int_depth;
  offsets[3] = 2 * mz_depth + int_depth;
  if(mz_depth == 4 && int_depth == 4){
    process_ms_data<float, float, unsigned, unsigned>((float**)mzBucketData, (float**)intBucketData, bucket_elements_scan,
    bucket_elements, (float**)mzSeqData, (float**)intSeqData, seq_elements_scan, seq_elements,
    main_output, all_integer_elements1, all_integer_elements2,
    offsets, options);
  }
  else if(mz_depth == 8 && int_depth == 8){
    process_ms_data<double, double, size_t, size_t>((double**)mzBucketData, (double**)intBucketData, bucket_elements_scan,
    bucket_elements, (double**)mzSeqData, (double**)intSeqData, seq_elements_scan, seq_elements,
    main_output, all_integer_elements1, all_integer_elements2,
    offsets, options);
  }
  else if(mz_depth == 8 && int_depth == 4){
    process_ms_data<double, float, size_t, unsigned>((double**)mzBucketData, (float**)intBucketData, bucket_elements_scan,
    bucket_elements, (double**)mzSeqData, (float**)intSeqData, seq_elements_scan, seq_elements,
    main_output, all_integer_elements1, all_integer_elements2,
    offsets, options);
  }
  else if(mz_depth == 4 && int_depth == 8){
    process_ms_data<float, double, unsigned, size_t>((float**)mzBucketData, (double**)intBucketData, bucket_elements_scan,
    bucket_elements, (float**)mzSeqData, (double**)intSeqData, seq_elements_scan, seq_elements,
    main_output, all_integer_elements1, all_integer_elements2,
    offsets, options);
  }
  else {
    return EXIT_FAILURE;
  }
  *out_size = 0;
  for(size_t i = 0; i < main_output.size(); i++){
    (*out_buf) = (unsigned char*)realloc(*out_buf, (*out_size) + main_output[i].size());
    memcpy((*out_buf) + (*out_size), main_output[i].data(), main_output[i].size());
    (*out_size) += main_output[i].size();
  }
  return EXIT_SUCCESS;
}

#endif
