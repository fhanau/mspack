#include <iostream>
#include "head.h"

int MZXMLCompare(const char* input1, const char* input2){
  // compare the original file and the decompressed file
  // note: in the compression and decompression processes,
  // only the mz-int pairs are extracted and compressed, so here we only need to check if the pairs are same or not
  // mz-int pairs are encoded with base 64
  XMLDocument ori_doc;
  XMLError error = ori_doc.LoadFile(input1);
  if(error != XML_SUCCESS){
    fprintf(stderr, "Error while reading %s\n", input1);
    return EXIT_FAILURE;
  }
  XMLElement * msRun_ori = ori_doc.RootElement()->FirstChildElement("msRun");
  XMLElement * scan_ori = msRun_ori -> FirstChildElement("scan");
  int scanCount_ori = 0;
  msRun_ori->QueryIntAttribute("scanCount",&scanCount_ori);

  XMLDocument dc_doc;
  error = dc_doc.LoadFile(input2);
  if(error != XML_SUCCESS){
    fprintf(stderr, "Error while reading %s\n", input2);
    return EXIT_FAILURE;
  }
  XMLElement * scan_dc = dc_doc.RootElement() -> FirstChildElement("msRun") ->
    FirstChildElement("scan");

  for(int i = 0; i < scanCount_ori; i++)
  {
    XMLElement* peaks_ori = scan_ori -> FirstChildElement("peaks");
    XMLElement* peaks_dc = scan_dc -> FirstChildElement("peaks");
    int peaksCount_ori;
    scan_ori -> QueryIntAttribute("peaksCount",&peaksCount_ori);
    if(peaksCount_ori == 0){
      scan_ori = scan_ori -> NextSiblingElement("scan");
      scan_dc = scan_dc -> NextSiblingElement("scan");
      continue;
    }
    const char * base64char_ori = peaks_ori -> GetText();
    const char * base64char_dc = peaks_dc -> GetText();
    if(!base64char_ori || !base64char_dc){
      std::cout << "Missing scans in at least one file\n";
      return EXIT_FAILURE;
    }
    if(strcmp(base64char_ori, base64char_dc) != 0)
    {
      std::cout << "wrong decompression at scanNum = " << i << std::endl;
      return EXIT_FAILURE;
    }
    // turn to next scan
    scan_ori = scan_ori -> NextSiblingElement("scan");
    scan_dc = scan_dc -> NextSiblingElement("scan");
  }
  std::cout << "No errors detected\n";
  return EXIT_SUCCESS;
}
