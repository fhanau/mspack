/*
FPMSComp: mass spectrometry data compression

Ruochen Yang
email: rcyang624@126.com
*/

#define _FILE_OFFSET_BITS 64
#define EXT_GZ ".mgz"
#define EXT_BSC ".bsc"
#define EXT_LEN 4

#include "tinyxml2.h"

using namespace tinyxml2;

char* base64(const void* binaryData, int len, int *flen);
unsigned char* unbase64(const char* ascii, int len, int *flen);

int MSPrint(const char* input, int x, int y);

int MSDecode(const char* input, const char* output, unsigned char is_gz);
int MZMLDecode(const char* input, const char* output, unsigned char is_gz);

int MZXMLCompare(const char* input1, const char* input2);
