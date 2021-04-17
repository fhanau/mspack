#ifndef MSPREDICT_DECODER_H
#define MSPREDICT_DECODER_H

static uint32_t assemble_uint4(unsigned char* base, uint64_t step) {
  uint32_t val = base[0];
  val += (1 << 8) * base[step];
  val += (1 << 16) * base[2 * step];
  val += (1 << 24) * base[3 * step];
  return val;
}

static uint64_t assemble_uint8_8(unsigned char* base, uint64_t step) {
  uint64_t val = base[0];
  val += (1UL << 8) * base[step];
  val += (1UL << 16) * base[2 * step];
  val += (1UL << 24) * base[3 * step];
  val += (1UL << 32) * base[4 * step];
  val += (1UL << 40) * base[5 * step];
  val += (1UL << 48) * base[6 * step];
  val += (1UL << 56) * base[7 * step];
  return val;
}

static uint64_t assemble_uint(unsigned char* base, uint64_t step, uint32_t size) {
  if (size == 4) {
    return assemble_uint4(base, step);
  }
  return assemble_uint8_8(base, step);
}

template <typename F, typename F2, typename U, typename U2>
void assembleSeqScan(U* mzData, U2* intData, unsigned char* mz_start,
  unsigned char* int_start, size_t peaksCount, size_t num_elem, MSOptions& options) {
  U mz = 0;
  for (size_t i = 0; i < peaksCount; i++) {
    mz += assemble_uint(mz_start + i, num_elem, sizeof(U));
    U2 intensity = assemble_uint(int_start + i, num_elem, sizeof(U2));

    mzData[i] = decode_lossy<F, U>(mz, options.mz_lossy_mode, options.mz_lossy_error);
    intData[i] = decode_lossy<F2, U2>(intensity, options.int_lossy_mode, options.int_lossy_error);
  }
}

template <typename F, typename F2, typename U, typename U2>
void assembleBucketScan(U* mzData, U2* intData, unsigned char** bucket_start,
  uint64_t* bucket_used, size_t peaksCount, size_t num_elem, MSOptions& options) {
  U mz = 0;
  F prev_float = 0.0;
  for (size_t i = 0; i < peaksCount; i++) {
    unsigned bucket = bucket_idx<F>(prev_float);

    unsigned char* bucket_search = bucket_start[bucket];

    unsigned num_used = bucket_used[bucket];

    mz += assemble_uint(bucket_search + num_used, num_elem, sizeof(U));
    U mz_decoded = decode_lossy<F, U>(mz, options.mz_lossy_mode, options.mz_lossy_error);
    bucket_used[bucket]++;

    //Look up mz
    F pos = *(F*)&mz_decoded;
    bucket = bucket_idx<F>(pos);
    num_used = bucket_used[bucket];
    if (i == peaksCount - 1) {
      bucket = 0;
      num_used = bucket_used[0] - 1;
    }
    unsigned char* int_bucket_search = bucket_start[bucket] + num_elem * sizeof(U);

    U2 intensity = assemble_uint(int_bucket_search + num_used, num_elem, sizeof(U2));

    prev_float = pos;

    mzData[i] = mz_decoded;
    intData[i] = decode_lossy<F2, U2>(intensity, options.int_lossy_mode, options.int_lossy_error);
  }

}

#endif
