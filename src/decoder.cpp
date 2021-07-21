#include <arpa/inet.h>
#include <cmath>
#include "head.h"
#include "utils.h"
#include "lossy.h"
#include "decoder.h"

//When using diff to verify that decoded output is equivalent, enable this to
//include newlines at a cost of larger files.
#define XML_DEBUG 0

void LossyInfo(MSOptions& options) {
  fprintf(stderr, "Decoding mass spectrometry data\n");
  fprintf(stderr, "m/z loss: %f\n", options.mz_lossy_error);
  fprintf(stderr, "intensity loss: %f\n", options.int_lossy_error);

  if(options.int_lossy_mode == lossless){
    fprintf(stderr, "Lossless intensities\n");
  }
  else if (options.int_lossy_mode == log_re){
    fprintf(stderr, "Log transform intensities\n");
  }
  //TODO: Mode, format (mzml, mzxml,)
}

int MSDecode(const char* input, const char* output, unsigned char is_gz){
  char cmd [PATH_MAX];
  if (is_gz) {
    strcpy(cmd, "gzip -dfk -S " EXT_GZ " ");
    strcat(cmd, input);
    strcat(cmd, EXT_GZ);
  }
  else {
    strcpy(cmd, "bsc d ");
    strcat(cmd, input);
    strcat(cmd, EXT_BSC);
    strcat(cmd, " ");
    strcat(cmd, input);
    strcat(cmd, " -T");
  }
#ifndef NDEBUG
  fprintf(stderr, "%s\n", cmd);
#endif
  system(cmd);
  XMLDocument doc;
  XMLError error = doc.LoadFile(input);
  if(error != XML_SUCCESS){
    fprintf(stderr, "Error while reading %s%s\n", input, is_gz ? EXT_GZ : EXT_BSC);
    return EXIT_FAILURE;
  }
  XMLElement * msRun = doc.RootElement()->FirstChildElement("msRun");
  if(!msRun){
    fprintf(stderr, "XML error\n");
    return EXIT_FAILURE;
  }
  XMLElement * scan = msRun -> FirstChildElement("scan");
  XMLElement * peaks = scan -> FirstChildElement("peaks");

  int scanCount = 0;
  msRun->QueryIntAttribute("scanCount",&scanCount);

  int precision = 0;
  peaks->QueryIntAttribute("precision", &precision);

  int value_size = precision / 8;
  int pair_size = value_size * 2;

  FILE* fp = fopen(input, "rb");
  if(!fp){
    fprintf(stderr, "Error: Could not open file\n");
    return EXIT_FAILURE;
  }

  fseek(fp, 0, SEEK_END);
  size_t raw_size = ftell(fp);
  rewind(fp);
  unsigned char* raw_data = (unsigned char*)malloc(raw_size);
  fread(raw_data, 1, raw_size, fp);
  uint64_t xml_size = strlen((char*)raw_data) + 1;
  unsigned char* sz_data = raw_data + xml_size;

  uint64_t num_buckets = 0;
  uint64_t num_bucket_elem = 0;
  uint64_t num_seq_elem = 0;
  MSOptions options;

  //Read information on lossy
  memcpy(&options, sz_data, sizeof(MSOptions));
  sz_data += sizeof(MSOptions);
  if(options.scans_only) {
    //metadata is required to reconstruct file
    fprintf(stderr, "Error: Invalid data or invalid options structure for mzXML or invalid file format\n");
    fclose(fp);
    return EXIT_FAILURE;
  }

  //Acquire the stream sizes
  memcpy(&num_buckets, sz_data, sizeof(uint64_t));
  memcpy(&num_bucket_elem, sz_data + 8, sizeof(uint64_t));
  memcpy(&num_seq_elem, sz_data + 16, sizeof(uint64_t));
  //all_integer_elements is unused

  LossyInfo(options);

  unsigned char* payload = sz_data + 4 * sizeof(uint64_t);
  unsigned char* mz_1_start = payload + num_buckets * 8;
  unsigned char* int_1_start = mz_1_start + num_bucket_elem * value_size;
  unsigned char* mz_2_start = int_1_start + num_bucket_elem * value_size;
  unsigned char* int_2_start = mz_2_start + num_seq_elem * value_size;
  //fprintf(stderr, "pred_file_sz: %llu\n", xml_size + 4 * sizeof(uint64_t) + (num_bucket_elem + num_seq_elem) * value_size * 2);

  uint64_t ms1_read = 0;
  uint64_t ms2_read = 0;

  uint64_t* bucket_used = (uint64_t*)calloc(num_buckets, sizeof(uint64_t));
  unsigned char** bucket_start = (unsigned char**)malloc(sizeof(unsigned char*) * num_buckets);

  if (num_buckets) {
    bucket_start[0] = mz_1_start;
    for(uint64_t i = 1; i < num_buckets; i++){
      bucket_start[i] = bucket_start[i - 1] + assemble_uint8_8(payload + i - 1, num_buckets);
    }
  }

  for(int i = 0; i < scanCount; i++)
  {
    peaks = scan -> FirstChildElement("peaks");
    int peaksCount;
    scan -> QueryIntAttribute("peaksCount",&peaksCount);
    if(!peaksCount){
      scan = scan->NextSiblingElement("scan");
      continue;
    }
    unsigned* fpdata = (unsigned*)malloc(pair_size * peaksCount);

    int msLevel;
    scan -> QueryIntAttribute("msLevel",&msLevel);
    if(msLevel == 1 && options.no_buckets){
      unsigned mz = 0;
        for (int j = 0; j < peaksCount; j++) {
          //Assemble a single data pair
          mz += assemble_uint(mz_1_start + ms1_read, num_bucket_elem, value_size);
          unsigned intensity = assemble_uint(int_1_start + ms1_read, num_bucket_elem, value_size);
          ms1_read++;

          unsigned mz_tmp = decode_lossy<float, unsigned>(mz, options.mz_lossy_mode, options.mz_lossy_error);
          fpdata[j * 2] = htonl(mz_tmp);
          intensity = decode_lossy<float, unsigned>(intensity, options.int_lossy_mode, options.int_lossy_error);
          fpdata[j * 2 + 1] = htonl(intensity);
        }
    }
    else if(msLevel == 1){
      unsigned mz = 0;
      float prev_float = 0.0;
      for (int j = 0; j < peaksCount; j++) {
        unsigned bucket = bucket_idx<float>(prev_float);
        unsigned char* bucket_search = bucket_start[bucket];

        unsigned num_used = bucket_used[bucket];

        mz += assemble_uint(bucket_search + num_used, num_bucket_elem, value_size);
        unsigned mz_decoded = decode_lossy<float, unsigned>(mz, options.mz_lossy_mode, options.mz_lossy_error);
        bucket_used[bucket]++;

        //Look up mz
        prev_float = *(float*)&mz_decoded;
        bucket = bucket_idx<float>(prev_float);
        num_used = bucket_used[bucket];
        if(j == peaksCount - 1){
          bucket = 0;
          num_used = bucket_used[0] - 1;
        }
        unsigned char* int_bucket_search = bucket_start[bucket] + num_bucket_elem * value_size;

        unsigned intensity = assemble_uint(int_bucket_search + num_used, num_bucket_elem, value_size);

        fpdata[j * 2] = htonl(mz_decoded);
        unsigned intensity2 = decode_lossy<float, unsigned>(intensity, options.int_lossy_mode, options.int_lossy_error);
        fpdata[j * 2 + 1] = htonl(intensity2);
      }
    }
    if(msLevel != 1){
      unsigned mz = 0;
        for (int j = 0; j < peaksCount; j++) {
          //Assemble a single data pair
          mz += assemble_uint(mz_2_start + ms2_read, num_seq_elem, value_size);
          unsigned intensity = assemble_uint(int_2_start + ms2_read, num_seq_elem, value_size);
          ms2_read++;

          unsigned mz_tmp = decode_lossy<float, unsigned>(mz, options.mz_lossy_mode, options.mz_lossy_error);
          fpdata[j * 2] = htonl(mz_tmp);
          intensity = decode_lossy<float, unsigned>(intensity, options.int_lossy_mode, options.int_lossy_error);
          fpdata[j * 2 + 1] = htonl(intensity);
        }
    }

    int base64_len = 0;
    char* recoded = base64(fpdata, peaksCount * pair_size, &base64_len);
    free(fpdata);
    peaks->SetText(recoded);
    free(recoded);
    scan = scan->NextSiblingElement("scan");
  }

  error = doc.SaveFile(output, XML_DEBUG == 0);
  if(error != XML_SUCCESS){
    fprintf(stderr, "Error while saving file\n");
    return EXIT_FAILURE;
  }
  //Delete intermediate file
  remove(input);
  return 0;
}

template <typename F, typename F2, typename U, typename U2>
int MZMLDecode2(XMLElement*& currentSpectrum, unsigned char*& sz_data,
  uint64_t mz_precision, uint64_t int_precision, MSOptions options);

int MZMLDecode(const char* input, const char* output, unsigned char is_gz){
  char cmd [PATH_MAX];
  if (is_gz) {
    strcpy(cmd, "gzip -dfk -S " EXT_GZ " ");
    strcat(cmd, input);
    strcat(cmd, EXT_GZ);
  }
  else {
    strcpy(cmd, "bsc d ");
    strcat(cmd, input);
    strcat(cmd, EXT_BSC);
    strcat(cmd, " ");
    strcat(cmd, input);
    strcat(cmd, " -T");
  }
#ifndef NDEBUG
  fprintf(stderr, "%s\n", cmd);
#endif
  system(cmd);
  XMLDocument doc;
  FILE* fp = fopen(input, "rb");
  if(!fp){
    fprintf(stderr, "Error: Could not open file\n");
    return EXIT_FAILURE;
  }

  uint64_t xml_offset = 0;
  fread(&xml_offset, 1, sizeof(xml_offset), fp);
  fseek(fp, 0, SEEK_END);
  uint64_t xml_size = ftell(fp) - xml_offset;
  fprintf(stderr, "xml_size: %llu\n", xml_size);
  char* xml_data = (char*)malloc(xml_size);
  fseek(fp, xml_offset, SEEK_SET);
  fread(xml_data, xml_size, 1, fp);

  XMLError error = doc.Parse(xml_data, xml_size);
  if(error != XML_SUCCESS){
    fprintf(stderr, "Error while reading %s%s\n", input, is_gz ? EXT_GZ : EXT_BSC);
    return EXIT_FAILURE;
  }
  free(xml_data);

  XMLElement* mzML = doc.RootElement();
  //Support both indexedML and non-indexedML files
  XMLElement* spectrumList = mzML->FirstChildElement("mzML");
  if(spectrumList){
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
  int scanCount;
  spectrumList->QueryIntAttribute("count", &scanCount);
  fprintf(stderr, "scanCount: %d\n", scanCount);

  XMLElement* currentSpectrum = spectrumList->FirstChildElement("spectrum");

  fseek(fp, sizeof(uint64_t), SEEK_SET);
  uint64_t raw_size = xml_offset - sizeof(uint64_t);
  unsigned char* raw_data = (unsigned char*)malloc(raw_size);
  fread(raw_data, 1, raw_size, fp);

  //Parse precision and options independently from scans
  MSOptions options;
  memcpy(&options, raw_data, sizeof(MSOptions));
  if(options.scans_only) {
    fprintf(stderr, "Error: Invalid data or invalid options structure for mzML or invalid file format\n");
    fclose(fp);
    //metadata is required to reconstruct file
    return EXIT_FAILURE;
  }
  unsigned char* sz_data = raw_data + sizeof(MSOptions);
  uint64_t mz_precision = 0;
  uint64_t int_precision = 0;
  memcpy(&mz_precision, sz_data, sizeof(uint64_t));
  memcpy(&int_precision, sz_data + 8, sizeof(uint64_t));
  fprintf(stderr, "mz_precision: %llu\n", mz_precision);
  fprintf(stderr, "int_precision: %llu\n", int_precision);
  int return_code = 0;
  while (return_code != 2) {
    if(mz_precision == 4 && int_precision == 4){
      return_code = MZMLDecode2<float, float, uint32_t, uint32_t>(currentSpectrum, sz_data, mz_precision, int_precision, options);
    }
    else if(mz_precision == 4 && int_precision == 8){
      return_code = MZMLDecode2<float, double, uint32_t, uint64_t>(currentSpectrum, sz_data, mz_precision, int_precision, options);
    }
    else if(mz_precision == 8 && int_precision == 4){
      return_code = MZMLDecode2<double, float, uint64_t, uint32_t>(currentSpectrum, sz_data, mz_precision, int_precision, options);
    }
    else if(mz_precision == 8 && int_precision == 8){
      return_code = MZMLDecode2<double, double, uint64_t, uint64_t>(currentSpectrum, sz_data, mz_precision, int_precision, options);
    }
    else {
      return EXIT_FAILURE;
    }
    if (return_code == EXIT_FAILURE) {
      return EXIT_FAILURE;
    }
  }

  if(doc.SaveFile(output, XML_DEBUG == 0) != XML_SUCCESS){
    fprintf(stderr, "Error while saving file\n");
    return EXIT_FAILURE;
  }
  free(raw_data);
  //Delete intermediate file
  remove(input);
  return 1;
}

template <typename F, typename F2, typename U, typename U2>
int MZMLDecode2(XMLElement*& currentSpectrum, unsigned char*& sz_data,
  uint64_t mz_precision, uint64_t int_precision, MSOptions options){
  uint64_t num_buckets = 0;
  uint64_t num_bucket_elem = 0;
  uint64_t num_seq_elem = 0;

  //Acquire the stream sizes
  uint64_t num_blocks = 0;
  memcpy(&num_blocks, sz_data + 16, sizeof(uint64_t));
  //fprintf(stderr, "num_blocks: %llu\n", num_blocks);
  sz_data += 8;
  memcpy(&num_buckets, sz_data + 16, sizeof(uint64_t));
  memcpy(&num_bucket_elem, sz_data + 24, sizeof(uint64_t));
  memcpy(&num_seq_elem, sz_data + 32, sizeof(uint64_t));
  //all_integer_elements is unused
#ifdef DEBUG_PRINT
  LossyInfo(options);
#endif
  unsigned char* payload = sz_data + 6 * sizeof(uint64_t);
  unsigned char* mz_1_start = payload + num_buckets * 8;
  unsigned char* int_1_start = mz_1_start + num_bucket_elem * mz_precision;
  unsigned char* mz_2_start = int_1_start + num_bucket_elem * int_precision;
  unsigned char* int_2_start = mz_2_start + num_seq_elem * mz_precision;
  sz_data = int_2_start + num_seq_elem * int_precision - 16;
#ifdef PRINT_DEBUG
  fprintf(stderr, "number of buckets: %llu\n", num_buckets);
#endif
  //fprintf(stderr, "pred_file_sz: %llu\n", xml_size + 6 * sizeof(uint64_t) + (num_bucket_elem + num_seq_elem) * (num_bucket_elem + num_seq_elem));

  uint64_t ms1_read = 0;
  uint64_t ms2_read = 0;

  uint64_t* bucket_used = (uint64_t*)calloc(num_buckets, sizeof(uint64_t));
  unsigned char** bucket_start = (unsigned char**)malloc(sizeof(unsigned char*) * num_buckets);

  if (num_buckets) {
    bucket_start[0] = mz_1_start;
    for(uint64_t i = 1; i < num_buckets; i++){
      bucket_start[i] = bucket_start[i - 1] + assemble_uint8_8(payload + i - 1, num_buckets);
    }
  }

  for(uint64_t i = 0; i < num_blocks; i++)
  {
    XMLElement* cvParam = currentSpectrum->FirstChildElement("cvParam");
    int msLevel = -1;
    while(cvParam){
      const char* name = 0;
      XMLError error = cvParam->QueryStringAttribute("name", &name);
      if (error != XML_SUCCESS) {
        cvParam = cvParam->NextSiblingElement();
      }
      if(strcmp(name, "ms level") == 0){
        cvParam->QueryIntAttribute("value", &msLevel);
        break;
      }
      cvParam = cvParam->NextSiblingElement();
    }

    int peaksCount;
    currentSpectrum->QueryIntAttribute("defaultArrayLength", &peaksCount);
    if(!peaksCount){
      currentSpectrum = currentSpectrum->NextSiblingElement();
      continue;
    }

    U* mzdata = (U*)malloc(mz_precision * peaksCount);
    U2* intdata = (U2*)malloc(int_precision * peaksCount);
    if (msLevel == 1 && options.no_buckets) {
      assembleSeqScan<F, F2, U, U2>(mzdata, intdata, mz_1_start + ms1_read, int_1_start + ms1_read, peaksCount, num_bucket_elem, options);
      ms1_read += peaksCount;
    }
    else if(msLevel == 1){
      assembleBucketScan<F, F2, U, U2>(mzdata, intdata, bucket_start, bucket_used,
        peaksCount, num_bucket_elem, options);
    }
    if(msLevel != 1){
      assembleSeqScan<F, F2, U, U2>(mzdata, intdata, mz_2_start + ms2_read, int_2_start + ms2_read, peaksCount, num_seq_elem, options);
      ms2_read += peaksCount;
    }

      XMLElement* binaryDataArrayList = currentSpectrum->FirstChildElement("binaryDataArrayList");
      int arrayCount;
      binaryDataArrayList->QueryIntAttribute("count", &arrayCount);
      XMLElement* binaryDataArray = binaryDataArrayList->FirstChildElement("binaryDataArray");

      for(int j = 0; j < arrayCount; j++){
        int encodedLength;
        binaryDataArray->QueryIntAttribute("encodedLength", &encodedLength);

        cvParam = binaryDataArray->FirstChildElement("cvParam");

        typedef enum {
          TBD,
          MZ,
          INT,
          ION
        } BinDataType;

        BinDataType type = TBD;
        unsigned precision = 0;
        while(cvParam){
          const char* name = 0;
          XMLError error = cvParam->QueryStringAttribute("name", &name);
          if (error != XML_SUCCESS) {
            cvParam = cvParam->NextSiblingElement();
            continue;
          }
          if(strcmp(name, "m/z array") == 0){
            type = MZ;
          }
          if(strcmp(name, "intensity array") == 0){
            type = INT;
          }
          if(strcmp(name, "non-standard data array") == 0){
            type = ION;
          }
          if(strcmp(name, "64-bit float") == 0){
            precision = 64;
          }
          if(strcmp(name, "32-bit float") == 0){
            precision = 32;
          }
          cvParam = cvParam->NextSiblingElement();
        }

        XMLElement* binary = binaryDataArray->FirstChildElement("binary");
        if(precision == mz_precision * 8 && type == MZ){
          int base64_len = 0;
          char* recoded = base64(mzdata, peaksCount * mz_precision, &base64_len);
          free(mzdata);
          binary->SetText(recoded);
          free(recoded);
        }
        if(precision == int_precision * 8 && type == INT){
          int base64_len = 0;
          char* recoded = base64(intdata, peaksCount * int_precision, &base64_len);
          free(intdata);
          binary->SetText(recoded);
          free(recoded);
        }
        binaryDataArray = binaryDataArray->NextSiblingElement();
      }

    currentSpectrum = currentSpectrum->NextSiblingElement();
  }

  free(bucket_used);
  free(bucket_start);

  if(currentSpectrum == NULL ) {
    return 2;
  }

  return EXIT_SUCCESS;
}
