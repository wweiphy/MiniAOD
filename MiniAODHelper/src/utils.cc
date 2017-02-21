#include <fstream>
#include <string>

#include "MiniAOD/MiniAODHelper/interface/utils.h"


bool utils::fileExists(const std::string& fileName) {
  std::ifstream infile(fileName.c_str());
  return infile.good();
}
