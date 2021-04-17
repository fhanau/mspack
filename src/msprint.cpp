#include <arpa/inet.h>
#include <cmath>
#include "head.h"

int MSPrint(const char* input, int x, int y){
  XMLDocument doc;
  XMLError error = doc.LoadFile(input);
  if(error != XML_SUCCESS){
    fprintf(stderr, "Error while reading file %s\n", input);
    return EXIT_FAILURE;
  }
  XMLElement * mzXML = doc.RootElement();

  XMLElement * msRun = mzXML -> FirstChildElement("msRun");

  XMLElement * scan = msRun -> FirstChildElement("scan");
  XMLElement * peaks = scan -> FirstChildElement("peaks");

  int scanCount = 0;
  msRun->QueryIntAttribute("scanCount",&scanCount);

  // check for scan precision
  int precision = 0;

  peaks->QueryIntAttribute("precision", &precision);

  int doubleprecision = precision == 64;
  int value_size = precision / 8;
  double total_sum = 0.0;
  double log_sum = 0.0;
  unsigned globalPeaks = 0;

  for(int i = 0; i < scanCount; i++)
  {
    peaks = scan -> FirstChildElement("peaks");
    int peaksCount;
    scan -> QueryIntAttribute("peaksCount",&peaksCount);

    const char * base64char = peaks->GetText();
    int len = strlen(base64char);
    int flen = 0;
    unsigned char * bin = unbase64(base64char, len, &flen);// credit to ...
    for(int j = 0; j < peaksCount; j++) {
      if(!doubleprecision){
        unsigned tmp, tmp2;
        memcpy(&tmp, bin + (2 * j) * value_size, 4);
        memcpy(&tmp2, bin + (2 * j + 1) * value_size, 4);

        tmp = ntohl(tmp);
        tmp2 = ntohl(tmp2);
        float intensity = *(float*)&tmp;
        total_sum += intensity;
        log_sum += log2(intensity);
        globalPeaks++;
        if (i == x && j == y) {
          printf("m/z at %d %d: %f\n", x, y, intensity);
        }
      }
    }

    if(i < scanCount-1){
      scan = scan -> NextSiblingElement();
    }
  }
  printf("Peaks: %d Avg intensity: %f Avg log: %f\n", globalPeaks, total_sum / globalPeaks, log_sum / globalPeaks);
  return 0;
}
