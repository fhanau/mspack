# mspack
mspack is a C++ program for lossless and lossy mass spectrometry data compression, achieving a high compression ratio without sacrificing performance. It includes an example implementation for mzXML and mzML as well as a format-agnostic API.

## Requirements
- C++11 compatible compiler, preferably gcc/clang.
- While a Unix environment is not required, compiling on Windows using other compilers might fail due to missing headers.
- [gzip](https://ftp.gnu.org/gnu/gzip/), should already be installed on a Unix system. On Linux, please make sure to use the latest gzip version (gzip 1.10). The gzip included in macOS works as well.
- [bsc](https://github.com/IlyaGrebnov/libbsc) for better, but slower compression (recommended)

## API
The API is defined in ```include/mspack.h```. The mspack binary represents an example implementation to showcase the capabilities of the program. Users are encouraged to use the API for advanced usage and to use the mspack transforms in other software.

## Binary Usage:
Compile by running make in src/ directory. This will place the binary in the top directory. Running the program without arguments will provide a help interface explaining the program usage, including lossy and advanced options. Use the .mgz file extension to use the gzip backend and .bsc to use bsc, which will improve compression significantly, but slow down compression.

mzML encode:
```mspack --mzmle (options) <in> <out><.mgz|.bsc>```

mzml decode:
```mspack --mzmld <in><.mgz|.bsc> <out>```

mzXML encode:
```mspack --mzxmle (options) <in> <out><.mgz|.bsc>```

mzXML decode:
```mspack --mzxmld <in><.mgz|.bsc> <out>```

## Example
The examples directory contains a test file which can be compressed using the following commands:
### Encode
```./mspack --mzmle examples/BSA1.mzml BSA.mgz```
### Encode file using bsc backend
```./mspack --mzmle examples/BSA1.mzml BSA.bsc```
### Encode lossily using default error
```./mspack --mzmle --mz-fixed-abs --int-log --examples/BSA1.mzml BSA-lossy.mgz```
### Decode
```./mspack --mzmld BSA.bsc BSA-decoded.mzml```

## Caveats
Due to the XML library we use, whitespace is not preserved in the decoded file. While the decoded file will have the same data as the original file, it will have different whitespace and therefore a different checksum. This can also cause warnings in mzml software if the indexedML format is used as indexing information becomes incorrect.

The current implementation depends on the mz values being increasing. This is sometimes not the case if a third data point is present, i.e. ion mobility in addition to mz and intensity.

The implementation can only use raw mz and intensity arrays. If the arrays are compressed with e.g. zlib or msnumpress, the file can be converted using msconvert so it can be used by mspack.

The mzXML implementation does not support all features as the format has been superseded by mzML. mzXML is limited to 32-bit files and does not support the block-based I/O feature.

### Licenses and Acknowledgements
mspack is available under the Apache License. It uses the tinyxml2 library by Lee Thomason, SHA1 by Dominik Reichl and base64 by William Sherif. Consult the source code for license information of these libraries. The mspack XML parsing code is partially based on the parsing code of [Masscomp](https://github.com/iochoa/MassComp) by Ruochen Yang, Xi Chen, and Idoia Ochoa.

### Authors
Felix Hanau and Idoia Ochoa, University of Illinois at Urbana-Champaign
