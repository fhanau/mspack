#include <getopt.h>
#include <inttypes.h>
#include <sys/stat.h>
#include "encoder.h"
#include "head.h"
#include "mspredict.h"

typedef enum mspredict_mode {
  encodeMZML,
  decodeMZML,
  encodeMZXML,
  decodeMZXML,
  printMZXML,
  cmpMZXML,
  invalid
} mspredict_mode;

void help(const char* prog){
  fprintf(stderr,
    "mspredict – mass spectrometry compressor\n"
    "Usage: %s <command> [Options] input output\n"
    "Commands:\n"
    " --mzmle  Compress MZML file\n"
    " --mzmld  Decompress MZML file\n"
    " --mzxmle Compress MZXML file\n"
    " --mzxmld Decompress MZXML file\n"
//Less stable option, undocumented for now
//    " --cmp    Verify that two MZXML files contain the same data\n"
    "\n"
    "Lossy options:\n"
    " --mz-lossless            Lossless m/z compression (default)\n"
    " --mz-trunc-abs=<error>   Truncate bits, absolute error (default error 10^-4)\n"
    " --mz-fixed-abs=<error>   Fixed point transform, absolute error (default error 10^-4)\n"
    "\n"
    " --int-lossless           Lossless intensity compression (default)\n"
//    " --int-integers           Truncate intensities to integers\n"
    " --int-trunc-rel=<error>  Truncate bits, relative error (default error 10^-2)\n"
    " --int-log=<error>        Log transform, relative error (default error 10^-2)\n"
    "\n"
    " --block-size=n Maximum number of blocks to store sequentially\n"
    " --no-buckets   Disable bucket transform\n"
    " --scans-only   Only write scan data (testing only)\n",
#ifdef __DATE__
    __DATE__,
#endif
    prog
  );
}

int main(int argc,char* argv[]) {
  mspredict_mode mode = invalid;
  MSOptions options;
  options.mz_lossy_mode = lossless;
  options.mz_lossy_error = 0.0;
  options.int_lossy_mode = lossless;
  options.int_lossy_error = 0.0;
  options.scans_only = 0;
  options.block_size = 0;
  options.no_buckets = 0;
  int c;

#define MZ_OPT_OFFSET 10
  struct option getopt_options[] = {
    //TODO: Add support for printing values for easier debugging
    //{"printmzxml", required_argument, (int*)&mode, printMZXML},

    {"cmpmzxml",         no_argument, (int*)&mode, cmpMZXML},
    {"mzmle",            no_argument, (int*)&mode, encodeMZML},
    {"mzmld",            no_argument, (int*)&mode, decodeMZML},
    {"mzxmle",           no_argument, (int*)&mode, encodeMZXML},
    {"mzxmld",           no_argument, (int*)&mode, decodeMZXML},

    {"int-lossless",  no_argument,       (int*)&options.int_lossy_mode, lossless},
//int-integers has been disabled as it does not work on all files in its current form. If the option is enabled, this will cause an intentional exception.
    {"int-integers",  no_argument,       (int*)&options.int_lossy_mode, truncate_to_int},
    {"int-trunc-rel", optional_argument, (int*)&options.int_lossy_mode, truncate_bits_re},
    {"int-log",       optional_argument, (int*)&options.int_lossy_mode, log_re},

    {"mz-lossless",  no_argument,       (int*)&options.mz_lossy_mode, lossless},
    {"mz-trunc-abs", optional_argument, (int*)&options.mz_lossy_mode, truncate_bits_ae},
    {"mz-fixed-abs", optional_argument, (int*)&options.mz_lossy_mode, fixed_ae},

    {"scans-only",   no_argument, &options.scans_only, 1},
    {"block-size",   required_argument, (int*)&options.block_size, 1},
    {"no-buckets",   no_argument, &options.no_buckets, 1},
    {"help",         no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };

  while (true) {
    int option_index = 0;
    double error;
    c = getopt_long(argc, argv, "", getopt_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {
      case 0:
        if (getopt_options[option_index].has_arg == no_argument) {
          break;
        }
        if(getopt_options[option_index].has_arg == required_argument) {
          intmax_t block_size = strtoimax(optarg, 0, 10);
          if(block_size < 0 || block_size > UINT_MAX) {
            fprintf(stderr, "Error: Block size must not be negative\n");
            return EXIT_FAILURE;
          }
          options.block_size = block_size;
          break;
        }
        if (!optarg){
          if(option_index < MZ_OPT_OFFSET){
            options.int_lossy_error = INT_DEFAULT_ERROR;
            fprintf(stderr, "info: no error parameter specified, using default error of 0.01\n");
          }
          else {
            options.mz_lossy_error = MZ_DEFAULT_ERROR;
            fprintf(stderr, "info: no error parameter specified, using default error of 0.0001\n");
          }
          continue;
        }

        error = strtod(optarg, 0);
        if(error <= 0.0 || error >= MAX_LOSS){
          fprintf(stderr, "Error out of range: must be between 0 and 1.\n");
          return EXIT_FAILURE;
        }

        LossyMode mode;
        if(option_index < MZ_OPT_OFFSET){
          mode = options.int_lossy_mode;
          options.int_lossy_error = error;
        }
        else {
          mode = options.mz_lossy_mode;
          options.mz_lossy_error = error;
        }
        if((mode == truncate_bits_re || mode == log_re) && error > MAX_REL_LOSS){
          fprintf(stderr, "Error out of range: Relative error must be between 0 and 0.5.\n");
          return EXIT_FAILURE;
        }
        break;
      case 'h':
          help(argv[0]);
          return EXIT_SUCCESS;
        break;
      case ':':
        fprintf(stderr, "missing option argument\n");
        return EXIT_FAILURE;
        break;
      case '?':
        break;
      default:
        return EXIT_FAILURE;
    }
  }

  if(options.mz_lossy_mode == lossless || options.mz_lossy_mode == truncate_to_int){
    options.mz_lossy_error = 0.0;
  }
  if(options.int_lossy_mode == lossless || options.int_lossy_mode == truncate_to_int){
    options.int_lossy_error = 0.0;
  }

  if (mode == encodeMZXML && options.block_size != 0) {
    fprintf(stderr, "error: using mzxml with block size is not supported at this time.\n");
  }

  //TODO
  /*if(mode == printMZXML){
    return MSPrint(argv[optind], i, j);
    return EXIT_FAILURE;
  }*/
  int is_gz = 0;
  if ((argc - optind == 2) && (mode == encodeMZML || mode == encodeMZXML)) {
    int out_len = strnlen(argv[optind + 1], PATH_MAX);
    if (out_len < (EXT_LEN+1) || (memcmp(EXT_GZ, &(argv[optind + 1][out_len - EXT_LEN]), EXT_LEN) != 0 && memcmp(EXT_BSC, &(argv[optind + 1][out_len - EXT_LEN]), EXT_LEN) != 0)) {
      fprintf(stderr, "Encoder out file name must end with '%s' or '%s'\n", EXT_GZ, EXT_BSC);
      return EXIT_FAILURE;
    }
    is_gz = memcmp(EXT_GZ, &(argv[optind + 1][out_len - EXT_LEN]), EXT_LEN) == 0;
    //Let gzip/bsc add the extension later
    argv[optind + 1][out_len - EXT_LEN] = 0;
    struct stat stats;
    if(stat(argv[optind + 1], &stats) == 0) {
      fprintf(stderr, "Intermediate file name is used, pick a different outname\n");
      return EXIT_FAILURE;
    }
  }
  if ((argc - optind == 2) && (mode == decodeMZML || mode == decodeMZXML)) {
    int out_len = strnlen(argv[optind], PATH_MAX);
    if (out_len < (EXT_LEN+1) || (memcmp(EXT_GZ, &(argv[optind][out_len - EXT_LEN]), EXT_LEN) != 0 &&  memcmp(EXT_BSC, &(argv[optind][out_len - EXT_LEN]), EXT_LEN) != 0)) {
      fprintf(stderr, "Decoder in file name must end with '%s' or '%s'\n", EXT_GZ, EXT_BSC);
      return EXIT_FAILURE;
    }
    //Let gzip/bsc add the extension later
    is_gz = memcmp(EXT_GZ, &(argv[optind][out_len - EXT_LEN]), EXT_LEN) == 0;
    argv[optind][out_len - EXT_LEN] = 0;
    struct stat stats;
    if(stat(argv[optind], &stats) == 0) {
      fprintf(stderr, "Intermediate file name is used, try moving %s\n", argv[optind]);
      return EXIT_FAILURE;
    }
  }

  if (options.scans_only) {
    fprintf(stderr, "warning: scans-only is for experimental purposes only as the deta can not be decoded without metadata");
  }

  if(mode == encodeMZML && argc - optind == 2){
    return MZMLEncode(argv[optind], argv[optind + 1], options, is_gz);
  }
  if(mode == decodeMZML && argc - optind == 2){
    return MZMLDecode(argv[optind], argv[optind + 1], is_gz);
  }
  if(mode == encodeMZXML && argc - optind == 2){
    return MSEncode(argv[optind], argv[optind + 1], options, is_gz);
  }
  if(mode == decodeMZXML && argc - optind == 2){
    return MSDecode(argv[optind], argv[optind + 1], is_gz);
  }
  if(mode == cmpMZXML && argc - optind == 2){
    return MZXMLCompare(argv[optind], argv[optind + 1]);
  }

  help(argv[0]);
  return EXIT_SUCCESS;
}
