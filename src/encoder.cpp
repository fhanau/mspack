#include <arpa/inet.h>
#include <assert.h>
#include "base64.h"
#include "encoder.h"
#include "head.h"
#include "mspredict.h"

int verify_loss(MSOptions& options, int mz_bit_depth, int int_bit_depth) {
  if (options.mz_lossy_mode == fixed_ae && (mz_bit_depth == 4 || options.mz_lossy_error < LOSSY_FIXED_PO2_CUTOFF)) {
    double lossy_error = options.mz_lossy_error;
    unsigned long le_integral = *(unsigned long*)&lossy_error;
    le_integral &= 0xFFF0000000000000UL;
    lossy_error = *(double*)&le_integral;
    options.mz_lossy_error = lossy_error;
  }
  if (options.mz_lossy_mode == fixed_ae || options.mz_lossy_mode == truncate_bits_ae) {
    if (options.mz_lossy_error < LOSSY_FIXED_PO2_CUTOFF && mz_bit_depth == 4) {
      return 1;
    }
    if (options.mz_lossy_error < LOSSY_MZ_CUTOFF) {
      return 1;
    }
  }
  if (options.int_lossy_mode == log_re || options.int_lossy_mode == truncate_bits_re) {
    if (options.int_lossy_mode == log_re && options.int_lossy_error < LOSSY_LOG_CUTOFF) {
      return 1;
    }
    if (options.int_lossy_error < LOSSY_REL_CUTOFF_SINGLE && int_bit_depth == 4) {
      return 1;
    }
    if (options.int_lossy_error < LOSSY_REL_CUTOFF_DOUBLE) {
      return 1;
    }
  }
  return 0;
}

template <typename F, typename U>
int MSEncode2(const char* output, XMLDocument* doc,
  XMLElement* msRun, XMLElement* scan, XMLElement* peaks, int precision, MSOptions& options, unsigned char is_gz);

int MSEncode(const char* input, const char* output, MSOptions& options, unsigned char is_gz){
  XMLDocument doc;
  XMLError error = doc.LoadFile(input);
  if(error != XML_SUCCESS){
    fprintf(stderr, "Error while reading %s\n", input);
    return EXIT_FAILURE;
  }

  XMLElement* msRun = doc.RootElement()->FirstChildElement("msRun");
  if(!msRun){
    fprintf(stderr, "XML error\n");
    return EXIT_FAILURE;
  }
  XMLElement * scan = msRun -> FirstChildElement("scan");
  XMLElement * peaks = scan -> FirstChildElement("peaks");

  // check for scan precision
  int precision = 0;
  peaks->QueryIntAttribute("precision", &precision);

  if (verify_loss(options, precision / 8, precision / 8) != 0) {
    fprintf(stderr, "Illegal lossy options\n");
    return 1;
  }

  if(precision == 32){
    return MSEncode2<float, int>(output, &doc, msRun, scan, peaks, precision, options, is_gz);
  }
  else if(precision == 64){
    fprintf(stderr, "64-bit MZXML is currently not supported.\n");
    //return MSEncode2<double, long>(output, &doc, msRun, scan, peaks, precision, options);
  }
  return EXIT_FAILURE;
}

typedef enum {
  TBD,
  MZ,
  INT,
  ION
} BinDataType;

int parseCv(XMLElement* cvParam, BinDataType& type, unsigned& precision){
  unsigned char uncompressed = 0;
  while(cvParam){
    const char* name;
    cvParam->QueryStringAttribute("name", &name);
    const char* accession;
    cvParam->QueryStringAttribute("accession", &accession);
    if(strcmp(accession, "MS:1000514") == 0 && strcmp(name, "m/z array") == 0){
      type = MZ;
    }
    if(strcmp(accession, "MS:1000515") == 0 && strcmp(name, "intensity array") == 0){
      type = INT;
    }
    if(strcmp(name, "non-standard data array") == 0){
      type = ION;
    }
    if(strcmp(accession, "MS:1000523") == 0 && strcmp(name, "64-bit float") == 0){
      precision = 64;
    }
    if(strcmp(accession, "MS:1000521") == 0 && strcmp(name, "32-bit float") == 0){
      precision = 32;
    }
    if(strcmp(accession, "MS:1000576") == 0 && strcmp(name, "no compression") == 0){
      uncompressed = 1;
    }
    cvParam = cvParam->NextSiblingElement();
  }
  if(uncompressed != 1){
    fprintf(stderr, "Could not verify that data array is uncompressed\n");
    return 1;
  }
  if(type == TBD){
    fprintf(stderr, "Could not parse binary array type\n");
    return 1;
  }
  if(precision == 0){
    fprintf(stderr, "Could not parse binary array precision\n");
    return 1;
  }
  return 0;
}

int verify_increasing_values(void* bin, uint64_t length, unsigned bit_depth){
  if(bit_depth == 4){
    float* data = (float*)bin;
    for(uint64_t i = 0; i < length - 1; i++){
      if(data[i] > data[i + 1]) {
        return 1;
      }
    }
  }
  if(bit_depth == 8){
    double* data = (double*)bin;
    for(uint64_t i = 0; i < length - 1; i++){
      if(data[i] > data[i + 1]) {
        return 1;
      }
    }
  }
  return 0;
}

int queryPrecision(XMLElement* currentSpectrum, BinDataType query){
  XMLElement* cvParam;

  XMLElement* binaryDataArrayList = currentSpectrum->FirstChildElement("binaryDataArrayList");
  int arrayCount;
  binaryDataArrayList->QueryIntAttribute("count", &arrayCount);
  XMLElement* binaryDataArray = binaryDataArrayList->FirstChildElement("binaryDataArray");
  for(int j = 0; j < arrayCount; j++){
    BinDataType type = TBD;
    unsigned precision = 0;
    cvParam = binaryDataArray->FirstChildElement("cvParam");
    if(parseCv(cvParam, type, precision) != 0){
      fprintf(stderr, "Failed to parse CV information\n");
      return -1;
    }
    if(type == query){
      return precision;
    }
    binaryDataArray = binaryDataArray->NextSiblingElement();
  }
  return -1;
}

template <typename F, typename F2, typename U, typename U2>
int MZMLEncode2(FILE* error_split,unsigned start_scan, unsigned scanCount, XMLElement*& currentSpectrum, uint64_t mz_precision, uint64_t int_precision, MSOptions& options);

template <typename F, typename F2, typename U, typename U2>
int MZMLCompressNScans(FILE* error_split, unsigned start_scan, unsigned scanCount, XMLElement*& currentSpectrum, uint64_t mz_precision, uint64_t int_precision, MSOptions& options);

int MZMLEncode(const char* input, const char* output, MSOptions& options, unsigned char is_gz){
  XMLDocument doc;
  XMLError error = doc.LoadFile(input);
  if(error != XML_SUCCESS){
    fprintf(stderr, "Error while reading %s\n", input);
    return EXIT_FAILURE;
  }

  XMLElement* mzML = doc.RootElement();
  //Support both indexedML and non-indexedML files, although we cannot maintain the indicies at this time.
  XMLElement* spectrumList = mzML->FirstChildElement("mzML");
  if(spectrumList){
    fprintf(stderr, "warning: mspredict currently has limited support for indexedML. Due to whitespace issues, processed files may have invalid indexedML indicies.\n");
    spectrumList = spectrumList->FirstChildElement("run")->FirstChildElement("spectrumList");
  }
  else {
    spectrumList = mzML->FirstChildElement("run")->FirstChildElement("spectrumList");
  }
  if(!spectrumList){
    fprintf(stderr, "XML error\n");
    return 1;
  }
  //Get count
  unsigned scanCount;
  spectrumList->QueryUnsignedAttribute("count", &scanCount);

  XMLElement* currentSpectrum = spectrumList->FirstChildElement("spectrum");
  uint64_t mz_bit_depth = queryPrecision(currentSpectrum, MZ);
  uint64_t int_bit_depth = queryPrecision(currentSpectrum, INT);
  if(mz_bit_depth == 0 || int_bit_depth == 0){
    fprintf(stderr, "Could not parse precision\n");
    return 1;
  }

  uint64_t mz_precision = mz_bit_depth / 8;
  uint64_t int_precision = int_bit_depth / 8;
  if (verify_loss(options, mz_precision, int_precision) != 0) {
    fprintf(stderr, "Illegal lossy options\n");
    return 1;
  }
  unsigned start_scan = 0;
  int return_code = 1;

  FILE* error_split = fopen(output, "wb");
  uint64_t xml_offset = 0;

  fwrite(&xml_offset, 1, sizeof(xml_offset), error_split);
  fwrite(&options, 1, sizeof(MSOptions), error_split);
  fwrite(&mz_precision, 1, sizeof(uint64_t), error_split);
  fwrite(&int_precision, 1, sizeof(uint64_t), error_split);

  while(start_scan < scanCount) {
    if(mz_precision == 4 && int_precision == 4){
      return_code = MZMLCompressNScans<float, float, uint32_t, uint32_t>(error_split, start_scan, scanCount, currentSpectrum, mz_precision, int_precision, options);
      //return_code = MZMLEncode2<float, float, uint32_t, uint32_t>(error_split, start_scan, scanCount, currentSpectrum, mz_precision, int_precision, options);
    }
    else if(mz_precision == 8 && int_precision == 4){
      return_code = MZMLCompressNScans<double, float, uint64_t, uint32_t>(error_split,start_scan, scanCount, currentSpectrum, mz_precision, int_precision, options);
    }
    else if(mz_precision == 4 && int_precision == 8){
      return_code = MZMLCompressNScans<float, double, uint32_t, uint64_t>(error_split, start_scan, scanCount, currentSpectrum, mz_precision, int_precision, options);
    }
    else if(mz_precision == 8 && int_precision == 8){
      return_code = MZMLCompressNScans<double, double, uint64_t, uint64_t>(error_split, start_scan, scanCount, currentSpectrum, mz_precision, int_precision, options);
    }
    else {
      return 1;
    }
    if (return_code != 0) {
      return return_code;
    }
    start_scan += options.block_size;
#ifdef DEBUG_PRINT
    fprintf(stderr, "%d\n", start_scan);
#endif
    if (options.block_size == 0) {
      break;
    }
  }
  xml_offset = ftell(error_split);
  fprintf(stderr, "xml offset: %zu\n", xml_offset);
  if(!options.scans_only){
    XMLPrinter printIO(0, true, 0);
    doc.Accept(&printIO);
    fwrite(printIO.CStr(), 1, printIO.CStrSize(), error_split);
  }
  rewind(error_split);
  fwrite(&xml_offset, 1, sizeof(xml_offset), error_split);
  fclose(error_split);

  char cmd [PATH_MAX];
  if (is_gz) {
    strcpy(cmd, "gzip -6f -S " EXT_GZ " ");
    strcat(cmd, output);
  }
  else {
    strcpy(cmd, "bsc e ");
    strcat(cmd, output);
    strcat(cmd, " ");
    strcat(cmd, output);
    strcat(cmd, EXT_BSC " -T");
  }
  fprintf(stderr, "%s\n", cmd);
  system(cmd);

  //Remove intermediate file
  if (!is_gz) {
    remove(output);
  }
  return return_code;
}

template <typename F>
struct ms_chunk {
  F** scans;
  size_t* scan_sizes;
  size_t size;
};

template <typename F>void initMSChunk(size_t size, ms_chunk<F>* chunk){
  chunk->scans = (F**)calloc(size, sizeof(F*));
  chunk->scan_sizes = (size_t*)malloc(size * sizeof(size_t));
  chunk->size = 0;
}

template <typename F>void freeMSChunk(ms_chunk<F>* chunk){
  free(chunk->scan_sizes);
  for(size_t i = 0; i < chunk->size; i++) {
    if(chunk->scans[i] != 0) {
      free(chunk->scans[i]);
    }
  }
  free(chunk->scans);
}

template <typename F, typename F2, typename U, typename U2>
int MZMLReadNScans(XMLElement*& currentSpectrum, unsigned num_blocks,
  unsigned mz_precision, unsigned int_precision, ms_chunk<F>* mz_bucket_scans,
  ms_chunk<F>* mz_seq_scans, ms_chunk<F2>* int_bucket_scans,
  ms_chunk<F2>* int_seq_scans) {
  for (unsigned idx = 0; idx < num_blocks; idx++) {
    XMLElement* cvParam = currentSpectrum->FirstChildElement("cvParam");
    unsigned msLevel = 0;
    while(cvParam){
      const char* name = 0;
      cvParam->QueryStringAttribute("name", &name);
      if(strcmp(name, "ms level") == 0){
        cvParam->QueryUnsignedAttribute("value", &msLevel);
        break;
      }
      cvParam = cvParam->NextSiblingElement();
    }

    unsigned peaksCount;
    currentSpectrum->QueryUnsignedAttribute("defaultArrayLength", &peaksCount);
    if(!peaksCount){
      currentSpectrum = currentSpectrum->NextSiblingElement();
      continue;
    }

    ms_chunk<F>* mz_scans;
    ms_chunk<F2>* int_scans;
    if (msLevel == 1) {
      mz_scans = mz_bucket_scans;
      int_scans = int_bucket_scans;
    }
    else {
      mz_scans = mz_seq_scans;
      int_scans = int_seq_scans;
    }
    mz_scans->scans[mz_scans->size] = (F*)malloc(peaksCount * mz_precision);
    mz_scans->scan_sizes[mz_scans->size] = peaksCount;
    mz_scans->size++;
    int_scans->scans[int_scans->size] = (F2*)malloc(peaksCount * int_precision);
    int_scans->scan_sizes[int_scans->size] = peaksCount;
    int_scans->size++;

    XMLElement* binaryDataArrayList = currentSpectrum->FirstChildElement("binaryDataArrayList");
    int arrayCount;
    binaryDataArrayList->QueryIntAttribute("count", &arrayCount);
    XMLElement* binaryDataArray = binaryDataArrayList->FirstChildElement("binaryDataArray");
    for(int j = 0; j < arrayCount; j++){
      int encodedLength;
      binaryDataArray->QueryIntAttribute("encodedLength", &encodedLength);

      BinDataType type = TBD;
      unsigned precision = 0;
      cvParam = binaryDataArray->FirstChildElement("cvParam");
      if(parseCv(cvParam, type, precision) != 0){
        fprintf(stderr, "Failed to parse CV information\n");
        return 1;
      }
      if((type == MZ && mz_precision != precision / 8) || (type == INT && int_precision != precision / 8)){
        fprintf(stderr, "Mismatching precision\n");
        return 1;
      }

      XMLElement* base64 = binaryDataArray->FirstChildElement("binary");
      const char* payload = base64->GetText();
      int dataLen = 0;
      unsigned char * bin = unbase64(payload, encodedLength, &dataLen);
      if(type == MZ){
        assert((unsigned)dataLen / mz_precision == peaksCount);
        memcpy(mz_scans->scans[mz_scans->size - 1], bin, dataLen);
        if (verify_increasing_values(mz_scans->scans[mz_scans->size - 1], peaksCount, mz_precision)){
          fprintf(stderr, "m/z values are not increasing, can not process file\n");
          return 1;
        }
      }
      if(type == INT){
        assert((unsigned)dataLen / int_precision == peaksCount);
        memcpy(int_scans->scans[int_scans->size - 1], bin, dataLen);
      }
      if(type == MZ || type == INT){
        base64->SetText("");
      }

      free(bin);
      binaryDataArray = binaryDataArray->NextSiblingElement();
    }
    currentSpectrum = currentSpectrum->NextSiblingElement();
  }
  return 0;
}

template <typename F, typename F2, typename U, typename U2> int MZMLCompressNScans(FILE* error_split, unsigned start_scan, unsigned scanCount, XMLElement*& currentSpectrum, uint64_t mz_precision, uint64_t int_precision, MSOptions& options) {
  size_t num_blocks = scanCount - start_scan;
  if (options.block_size && num_blocks > options.block_size) {
    num_blocks = options.block_size;
  }
  ms_chunk<F> mz_bucket_scans;
  ms_chunk<F> mz_seq_scans;
  ms_chunk<F2> int_bucket_scans;
  ms_chunk<F2> int_seq_scans;

  initMSChunk(num_blocks, &mz_bucket_scans);
  initMSChunk(num_blocks, &mz_seq_scans);
  initMSChunk(num_blocks, &int_bucket_scans);
  initMSChunk(num_blocks, &int_seq_scans);

  int error = MZMLReadNScans<F, F2, U, U2>(currentSpectrum, num_blocks, mz_precision, int_precision, &mz_bucket_scans, &mz_seq_scans, &int_bucket_scans, &int_seq_scans);
  if (error) {
    return error;
  }

  std::vector<std::vector<unsigned char> > main_output(8 + (mz_precision + int_precision) * 2);

  uint64_t all_integer_intensities1 = 1;
  uint64_t all_integer_intensities2 = 1;
  uint64_t all_integer_elements = 0;

  //Use variables for sizes shared between m/z and intensities
  size_t bucket_elements = mz_bucket_scans.size;
  size_t* bucket_elements_scan = mz_bucket_scans.scan_sizes;
  if(!options.no_buckets){
  process_mz_bucket<F, U>(mz_bucket_scans.scans, bucket_elements_scan, bucket_elements, main_output, all_integer_elements, 0, options.mz_lossy_mode, options.mz_lossy_error);

  F** lossyMzBucketData = mz_bucket_scans.scans;
  if(options.mz_lossy_mode != lossless){
    lossyMzBucketData = (F**)malloc(bucket_elements * sizeof(F*));
    for(size_t i = 0; i < bucket_elements; i++){
      lossyMzBucketData[i] = (F*)malloc(bucket_elements_scan[i] * sizeof(F));
      for (size_t j = 0; j < bucket_elements_scan[i]; j++) {
        U val = reinterpret_cast<U&>(mz_bucket_scans.scans[i][j]);
        lossyMzBucketData[i][j] = apply_lossF<F, U>(val, options.mz_lossy_mode, options.mz_lossy_error);

        #if 0
                if(mz_bucket_scans.scans[i][j] > 1.0){
                  static int cnt2 = 0;
                  cnt2++;
                  static double maxerror2 = 0.0;
                  static double error_sum2 = 0.0;
                  double original = mz_bucket_scans.scans[i][j];
                  double newV = lossyMzBucketData[i][j];
                  double error = abs(newV - original);
                  error_sum2 += error;
                  if(error > maxerror2){
                    maxerror2 = error;
                  }

                  if(cnt2 % 100000 == 0 && cnt2 > 200000){
                    fprintf(stderr, "max absolute error: %.12f\n", maxerror2);
                    fprintf(stderr, "avg absolute error: %.12f\n", error_sum2 / cnt2);
                    fprintf(stderr, "count: %d\n", cnt2);
                  }
                }
        #endif

      }
    }
  }
  #if 0
    if(options.int_lossy_mode != lossless){
      for(size_t i = 0; i < bucket_elements; i++){
        for (size_t j = 0; j < bucket_elements_scan[i]; j++) {
          U2 val = reinterpret_cast<U2&>(int_bucket_scans.scans[i][j]);
          F2 newV = apply_lossF<F2, U2>(val, options.int_lossy_mode, options.int_lossy_error);
          if(int_bucket_scans.scans[i][j] > 1.0){
            static int cnt = 0;
            cnt++;
            static double maxerror = 0.0;
            static double error_sum = 0.0;
            F2 original = int_bucket_scans.scans[i][j];
            double error = abs((newV - original)/original);
            error_sum += error;
            if(error > maxerror){
              maxerror = error;
            }

            if(cnt % 100000 == 0 && cnt > 200000){
              fprintf(stderr, "max relative error: %.12f\n", maxerror);
              fprintf(stderr, "avg relative error: %.12f\n", error_sum / cnt);
              fprintf(stderr, "count: %d\n", cnt);
            }
          }
        }
      }
    }
  #endif

  process_int_bucket<F, F2, U2>(lossyMzBucketData, int_bucket_scans.scans, bucket_elements_scan, bucket_elements, main_output, all_integer_intensities1, 8 + mz_precision, options.int_lossy_mode, options.int_lossy_error);
  if(options.mz_lossy_mode != lossless && bucket_elements){
    for(size_t i = 0; i < bucket_elements; i++){
      if (lossyMzBucketData[i]) {
        free(lossyMzBucketData[i]);
      }
    }
    free(lossyMzBucketData);
  }
  }
  else {
    //The stream offset needs to be nonzero to produce the right bucket size
    process_generic_seq<F, U>(mz_bucket_scans.scans, bucket_elements_scan, bucket_elements, main_output, all_integer_elements, true, 8, options.mz_lossy_mode, options.mz_lossy_error);
    process_generic_seq<F2, U2>(int_bucket_scans.scans, bucket_elements_scan, bucket_elements, main_output, all_integer_intensities1, false, 8 + mz_precision, options.int_lossy_mode, options.int_lossy_error);
  }
  freeMSChunk(&mz_bucket_scans);
  freeMSChunk(&int_bucket_scans);

  process_generic_seq<F, U>(mz_seq_scans.scans, mz_seq_scans.scan_sizes, mz_seq_scans.size, main_output, all_integer_elements, true, 8 + mz_precision + int_precision, options.mz_lossy_mode, options.mz_lossy_error);
  process_generic_seq<F2, U2>(int_seq_scans.scans, int_seq_scans.scan_sizes, int_seq_scans.size, main_output, all_integer_intensities2, false, 8 + mz_precision * 2 + int_precision, options.int_lossy_mode, options.int_lossy_error);
  all_integer_elements = all_integer_intensities1 * 2 + all_integer_intensities2;

  freeMSChunk(&mz_seq_scans);
  freeMSChunk(&int_seq_scans);

  uint64_t num_buckets = main_output[0].size();
  uint64_t num_bucket_elem = main_output[8].size();
  uint64_t num_seq_elem = main_output[8 + mz_precision + int_precision].size();

  fwrite(&num_blocks, 1, sizeof(uint64_t), error_split);
  fwrite(&num_buckets, 1, sizeof(uint64_t), error_split);
  fwrite(&num_bucket_elem, 1, sizeof(uint64_t), error_split);
  fwrite(&num_seq_elem, 1, sizeof(uint64_t), error_split);
  fwrite(&all_integer_elements, 1, sizeof(uint64_t), error_split);

  for(uint64_t i = 0; i < main_output.size(); i++){
    fwrite(main_output[i].data(), 1, main_output[i].size(), error_split);
  }
  return 0;
}

template <typename F, typename F2, typename U, typename U2>
int MZMLEncode2(FILE* error_split, unsigned start_scan, unsigned scanCount, XMLElement*& currentSpectrum, uint64_t mz_precision, uint64_t int_precision, MSOptions& options){
  size_t num_blocks = 0;

  F** mzBucketData = (F**)malloc(scanCount * sizeof(F*));
  F2** intBucketData = (F2**)malloc(scanCount * sizeof(F2*));
  F** mzSeqData = (F**)malloc(scanCount * sizeof(F*));
  F2** intSeqData = (F2**)malloc(scanCount * sizeof(F2*));
  size_t bucket_elements = 0;
  size_t seq_elements = 0;
  size_t bucket_elements_scan[scanCount];
  size_t seq_elements_scan[scanCount];

  for(unsigned i = start_scan; i < scanCount; i++){
    XMLElement* cvParam = currentSpectrum->FirstChildElement("cvParam");
    unsigned msLevel = 0;
    while(cvParam){
      const char* name;
      cvParam->QueryStringAttribute("name", &name);
      if(strcmp(name, "ms level") == 0){
        cvParam->QueryUnsignedAttribute("value", &msLevel);
        break;
      }
      cvParam = cvParam->NextSiblingElement();
    }

    unsigned peaksCount;
    currentSpectrum->QueryUnsignedAttribute("defaultArrayLength", &peaksCount);
    if(!peaksCount){
      currentSpectrum = currentSpectrum->NextSiblingElement();
      //TODO
    num_blocks++;
    if (num_blocks == options.block_size) {
      break;
    }
      continue;
    }

    F* append_here_mz;
    F2* append_here_int;
    if(msLevel == 1){
      mzBucketData[bucket_elements] = (F*)malloc(peaksCount * mz_precision);
      intBucketData[bucket_elements] = (F2*)malloc(peaksCount * int_precision);
      append_here_mz = mzBucketData[bucket_elements];
      append_here_int = intBucketData[bucket_elements];
      bucket_elements_scan[bucket_elements] = peaksCount;
      bucket_elements++;
    }
    else {
      mzSeqData[seq_elements] = (F*)malloc(peaksCount * mz_precision);
      intSeqData[seq_elements] = (F2*)malloc(peaksCount * int_precision);
      append_here_mz = mzSeqData[seq_elements];
      append_here_int = intSeqData[seq_elements];
      seq_elements_scan[seq_elements] = peaksCount;
      seq_elements++;
    }

    XMLElement* binaryDataArrayList = currentSpectrum->FirstChildElement("binaryDataArrayList");
    int arrayCount;
    binaryDataArrayList->QueryIntAttribute("count", &arrayCount);
    XMLElement* binaryDataArray = binaryDataArrayList->FirstChildElement("binaryDataArray");
    for(int j = 0; j < arrayCount; j++){
      int encodedLength;
      binaryDataArray->QueryIntAttribute("encodedLength", &encodedLength);

      BinDataType type = TBD;
      unsigned precision = 0;
      cvParam = binaryDataArray->FirstChildElement("cvParam");
      if(parseCv(cvParam, type, precision) != 0){
        fprintf(stderr, "Failed to parse CV information\n");
        return 1;
      }
      if((type == MZ && mz_precision != precision / 8) || (type == INT && int_precision != precision / 8)){
        fprintf(stderr, "Mismatching precision\n");
        return 1;
      }

      XMLElement* base64 = binaryDataArray->FirstChildElement("binary");
      const char* payload = base64->GetText();
      int dataLen = 0;
      unsigned char * bin = unbase64(payload, encodedLength, &dataLen);
      if(type == MZ){
        assert((unsigned)dataLen / mz_precision == peaksCount);
        memcpy(append_here_mz, bin, dataLen);
        if (verify_increasing_values(append_here_mz, peaksCount, mz_precision)){
          fprintf(stderr, "m/z values are not increasing, can not process file\n");
          return 1;
        }
      }
      if(type == INT){
        assert((unsigned)dataLen / int_precision == peaksCount);
        memcpy(append_here_int, bin, dataLen);
      }
      if(type == MZ || type == INT){
        base64->SetText("");
      }

      free(bin);
      binaryDataArray = binaryDataArray->NextSiblingElement();
    }
    currentSpectrum = currentSpectrum->NextSiblingElement();
    num_blocks++;
    if (num_blocks == options.block_size) {
      break;
    }
  }

  std::vector<std::vector<unsigned char> > main_output(8 + (mz_precision + int_precision) * 2);

  uint64_t all_integer_intensities1 = 1;
  uint64_t all_integer_intensities2 = 1;
  uint64_t all_integer_elements = 0;

  if(!options.no_buckets){
  process_mz_bucket<F, U>(mzBucketData, bucket_elements_scan, bucket_elements, main_output, all_integer_elements, 0, options.mz_lossy_mode, options.mz_lossy_error);

  F** lossyMzBucketData = mzBucketData;
  if(options.mz_lossy_mode != lossless){
    lossyMzBucketData = (F**)malloc(bucket_elements * sizeof(F*));
    for(size_t i = 0; i < bucket_elements; i++){
      lossyMzBucketData[i] = (F*)malloc(bucket_elements_scan[i] * sizeof(F));
      for (size_t j = 0; j < bucket_elements_scan[i]; j++) {
        U val = reinterpret_cast<U&>(mzBucketData[i][j]);
        lossyMzBucketData[i][j] = apply_lossF<F, U>(val, options.mz_lossy_mode, options.mz_lossy_error);

#if 0
        if(mzBucketData[i][j] > 1.0){
          static int cnt2 = 0;
          cnt2++;
          static double maxerror2 = 0.0;
          static double error_sum2 = 0.0;
          double original = mzBucketData[i][j];
          double newV = lossyMzBucketData[i][j];
          double error = abs(newV - original);
          error_sum2 += error;
          if(error > maxerror2){
            maxerror2 = error;
          }

          if(cnt2 % 10000 == 0){
            fprintf(stderr, "max absolute error: %.12f\n", maxerror2);
            fprintf(stderr, "avg absolute error: %.12f\n", error_sum2 / cnt2);
          }
        }
#endif
      }
    }
  }

#if 0
  if(options.int_lossy_mode != lossless){
    for(size_t i = 0; i < bucket_elements; i++){
      for (size_t j = 0; j < bucket_elements_scan[i]; j++) {
        U2 val = reinterpret_cast<U2&>(intBucketData[i][j]);
        F2 newV = apply_lossF<F2, U2>(val, options.int_lossy_mode, options.int_lossy_error);
        if(intBucketData[i][j] > 1.0){
          static int cnt = 0;
          cnt++;
          static double maxerror = 0.0;
          static double error_sum = 0.0;
          F2 original = intBucketData[i][j];
          double error = abs((newV - original)/original);
          error_sum += error;
          if(error > maxerror){
            maxerror = error;
          }

          if(cnt % 100000 == 0){
            fprintf(stderr, "max relative error: %.12f\n", maxerror);
            fprintf(stderr, "avg relative error: %.12f\n", error_sum / cnt);
          }
        }
      }
    }
  }
#endif

  process_int_bucket<F, F2, U2>(lossyMzBucketData, intBucketData, bucket_elements_scan, bucket_elements, main_output, all_integer_intensities1, 8 + mz_precision, options.int_lossy_mode, options.int_lossy_error);
  if(options.mz_lossy_mode != lossless){
    for(size_t i = 0; i < bucket_elements; i++){
      free(lossyMzBucketData[i]);
    }
    free(lossyMzBucketData);
  }
  }
  else {
    //The stream offset needs to be nonzero to produce the right bucket size
    process_generic_seq<F, U>(mzBucketData, bucket_elements_scan, bucket_elements, main_output, all_integer_elements, true, 8, options.mz_lossy_mode, options.mz_lossy_error);
    process_generic_seq<F2, U2>(intBucketData, bucket_elements_scan, bucket_elements, main_output, all_integer_intensities1, false, 8 + mz_precision, options.int_lossy_mode, options.int_lossy_error);
  }

  process_generic_seq<F, U>(mzSeqData, seq_elements_scan, seq_elements, main_output, all_integer_elements, true, 8 + mz_precision + int_precision, options.mz_lossy_mode, options.mz_lossy_error);
  process_generic_seq<F2, U2>(intSeqData, seq_elements_scan, seq_elements, main_output, all_integer_intensities2, false, 8 + mz_precision * 2 + int_precision, options.int_lossy_mode, options.int_lossy_error);
  all_integer_elements = all_integer_intensities1 * 2 + all_integer_intensities2;

  free(mzBucketData);
  free(mzSeqData);
  free(intBucketData);
  free(intSeqData);

  uint64_t num_buckets = main_output[0].size();
  uint64_t num_bucket_elem = main_output[8].size();
  uint64_t num_seq_elem = main_output[8 + mz_precision + int_precision].size();

  if (num_blocks == 0) {
    num_blocks = scanCount;
  }
  fwrite(&num_blocks, 1, sizeof(uint64_t), error_split);
  fwrite(&options, 1, sizeof(MSOptions), error_split);
  fwrite(&num_buckets, 1, sizeof(uint64_t), error_split);
  fwrite(&num_bucket_elem, 1, sizeof(uint64_t), error_split);
  fwrite(&num_seq_elem, 1, sizeof(uint64_t), error_split);
  fwrite(&all_integer_elements, 1, sizeof(uint64_t), error_split);

  for(uint64_t i = 0; i < main_output.size(); i++){
    fwrite(main_output[i].data(), 1, main_output[i].size(), error_split);
  }
  return 0;
}

template <typename F, typename U>
int MSEncode2(const char* output, XMLDocument* doc, XMLElement* msRun,
  XMLElement* scan, XMLElement* peaks, int precision, MSOptions& options, unsigned char is_gz){
  int value_size = precision / 8;
  int scanCount = 0;
  msRun->QueryIntAttribute("scanCount",&scanCount);

  FILE* error_split = fopen(output, "wb");

  U ** mzArrays = (U**)malloc(sizeof(U*) * scanCount);
  U ** intArrays = (U**)malloc(sizeof(U*) * scanCount);

  F** mzBucketData = (F**)malloc(scanCount * sizeof(F*));
  F** intBucketData = (F**)malloc(scanCount * sizeof(F*));
  F** mzSeqData = (F**)malloc(scanCount * sizeof(F*));
  F** intSeqData = (F**)malloc(scanCount * sizeof(F*));
  size_t bucket_elements = 0;
  size_t seq_elements = 0;
  size_t bucket_elements_scan[scanCount];
  size_t seq_elements_scan[scanCount];

for(int i = 0; i < scanCount; i++)
{
  peaks = scan -> FirstChildElement("peaks");
  int peaksCount;
  scan -> QueryIntAttribute("peaksCount",&peaksCount);
  int msLevel;
  scan -> QueryIntAttribute("msLevel",&msLevel);
  if(!peaksCount){
    scan = scan -> NextSiblingElement();
    continue;
  }

  F* append_here_mz, *append_here_int;
  if(msLevel == 1){
    mzBucketData[bucket_elements] = (F*)malloc(peaksCount * value_size);
    intBucketData[bucket_elements] = (F*)malloc(peaksCount * value_size);
    append_here_mz = mzBucketData[bucket_elements];
    append_here_int = intBucketData[bucket_elements];
    bucket_elements_scan[bucket_elements] = peaksCount;
    bucket_elements++;
  }
  else {
    mzSeqData[seq_elements] = (F*)malloc(peaksCount * value_size);
    intSeqData[seq_elements] = (F*)malloc(peaksCount * value_size);
    append_here_mz = mzSeqData[seq_elements];
    append_here_int = intSeqData[seq_elements];
    seq_elements_scan[seq_elements] = peaksCount;
    seq_elements++;
  }

  const char * base64char = peaks->GetText();
  int len = strlen(base64char);
  int flen = 0;
  unsigned char * bin = unbase64(base64char, len, &flen);// credit to ...
  for(int j = 0; j < peaksCount; j++) {
      U tmp, tmp2;
      memcpy(&tmp, bin + (2 * j) * value_size, value_size);
      memcpy(&tmp2, bin + (2 * j + 1) * value_size, value_size);

      if (value_size == 8) {
        tmp = mspredict_ntohll(tmp);
        tmp2 = mspredict_ntohll(tmp2);
      }
      else {
        tmp = ntohl(tmp);
        tmp2 = ntohl(tmp2);
      }

      memcpy(&append_here_mz[j], &tmp, value_size);
      memcpy(&append_here_int[j], &tmp2, value_size);
  }
  if (verify_increasing_values(append_here_mz, peaksCount, value_size)) {
    fprintf(stderr, "m/z values are not increasing, can not process file\n");
    return 1;
  }
  free(bin);

  peaks -> SetText("");
  scan = scan -> NextSiblingElement();
}
  std::vector<std::vector<unsigned char> > main_output(8 + 4 * value_size);

  uint64_t all_integer_intensities1 = 1;
  uint64_t all_integer_intensities2 = 1;
  uint64_t all_integer_elements = 0;

  if(!options.no_buckets){
  process_mz_bucket<F, U>(mzBucketData, bucket_elements_scan, bucket_elements, main_output, all_integer_elements, 0, options.mz_lossy_mode, options.mz_lossy_error);

  F** lossyMzBucketData = mzBucketData;
  if(options.mz_lossy_mode != lossless){
    lossyMzBucketData = (F**)malloc(bucket_elements * sizeof(F*));
    for(size_t i = 0; i < bucket_elements; i++){
      lossyMzBucketData[i] = (F*)malloc(bucket_elements_scan[i] * sizeof(F));
      for (size_t j = 0; j < bucket_elements_scan[i]; j++) {
        U val = reinterpret_cast<U&>(mzBucketData[i][j]);
        lossyMzBucketData[i][j] = apply_lossF<F, U>(val, options.mz_lossy_mode, options.mz_lossy_error);
        #if 0
                if(mzBucketData[i][j] > 1.0){
                  static int cnt2 = 0;
                  cnt2++;
                  static double maxerror2 = 0.0;
                  static double error_sum2 = 0.0;
                  double original = mzBucketData[i][j];
                  double newV = lossyMzBucketData[i][j];
                  double error = abs(newV - original);
                  error_sum2 += error;
                  if(error > maxerror2){
                    maxerror2 = error;
                  }

                  if(cnt2 % 100000 == 0){
                    fprintf(stderr, "max absolute error: %.12f\n", maxerror2);
                    fprintf(stderr, "avg absolute error: %.12f\n", error_sum2 / cnt2);
                  }
                }
        #endif
      }
    }
  }

#if 0
if(options.int_lossy_mode != lossless){
  for(size_t i = 0; i < bucket_elements; i++){
    for (size_t j = 0; j < bucket_elements_scan[i]; j++) {
      U val = reinterpret_cast<U&>(intBucketData[i][j]);
      F newV = apply_lossF<F, U>(val, options.int_lossy_mode, options.int_lossy_error);
      if(intBucketData[i][j] > 1.0){
        static int cnt = 0;
        cnt++;
        static double maxerror = 0.0;
        static double error_sum = 0.0;
        F original = intBucketData[i][j];
        double error = abs(((double)newV - (double)original)/(double)original);
        error_sum += error;
        if(error > maxerror){
          maxerror = error;
        }

        if(cnt % 100000 == 0){
          fprintf(stderr, "max relative error: %.12f\n", maxerror);
          fprintf(stderr, "avg relative error: %.12f\n", error_sum / cnt);
        }
      }
    }
  }
}
#endif
  process_int_bucket<F, F, U>(lossyMzBucketData, intBucketData, bucket_elements_scan, bucket_elements, main_output, all_integer_intensities1, 8 + value_size, options.int_lossy_mode, options.int_lossy_error);
  if(options.mz_lossy_mode != lossless) {
    for(size_t i = 0; i < bucket_elements; i++){
      free(lossyMzBucketData[i]);
    }
    free(lossyMzBucketData);
  }
  }
  else {
    //The stream offset needs to be nonzero to produce the right bucket size
    process_generic_seq<F, U>(mzBucketData, bucket_elements_scan, bucket_elements, main_output, all_integer_elements, true, 8, options.mz_lossy_mode, options.mz_lossy_error);
    process_generic_seq<F, U>(intBucketData, bucket_elements_scan, bucket_elements, main_output, all_integer_intensities1, false, 8 + value_size, options.int_lossy_mode, options.int_lossy_error);
  }

  process_generic_seq<F, U>(mzSeqData, seq_elements_scan, seq_elements, main_output, all_integer_elements, true, 8 + 2 * value_size, options.mz_lossy_mode, options.mz_lossy_error);
  process_generic_seq<F, U>(intSeqData, seq_elements_scan, seq_elements, main_output, all_integer_intensities2, false, 8 + 3 * value_size, options.int_lossy_mode, options.int_lossy_error);
  all_integer_elements = all_integer_intensities1 * 2 + all_integer_intensities2;

  //TODO: Free the new structures
  for(int i = 0; i < scanCount; i++){
    free(mzArrays[i]);
    free(intArrays[i]);
  }

  free(mzArrays);
  free(intArrays);

  if(!options.scans_only){
    XMLPrinter printIO(0, true, 0);
    doc->Accept(&printIO);
    fwrite(printIO.CStr(), 1, printIO.CStrSize(), error_split);
  }

  fwrite(&options, 1, sizeof(MSOptions), error_split);
  uint64_t num_buckets = main_output[0].size();
  uint64_t num_bucket_elem = main_output[8].size();
  uint64_t num_seq_elem = main_output[8 + 2 * value_size].size();

  fwrite(&num_buckets, 1, sizeof(uint64_t), error_split);
  fwrite(&num_bucket_elem, 1, sizeof(uint64_t), error_split);
  fwrite(&num_seq_elem, 1, sizeof(uint64_t), error_split);
  fwrite(&all_integer_elements, 1, sizeof(uint64_t), error_split);

  for(uint64_t i = 0; i < main_output.size(); i++){
    fwrite(main_output[i].data(), 1, main_output[i].size(), error_split);
  }
  fclose(error_split);

  char cmd [PATH_MAX];
  if (is_gz) {
    strcpy(cmd, "gzip -6f -S " EXT_GZ " ");
    strcat(cmd, output);
  }
  else {
    strcpy(cmd, "bsc e ");
    strcat(cmd, output);
    strcat(cmd, " ");
    strcat(cmd, output);
    strcat(cmd, EXT_BSC " -T");
  }
  fprintf(stderr, "%s\n", cmd);
  system(cmd);

  //Remove intermediate file
  if (!is_gz) {
    remove(output);
  }
  return 0;
}
