//
// Created by Tsz-Chiu Au on 6/22/22.
//

#include "debug.h"

#include <sstream>
#include <filesystem>
#include <algorithm>


// source_location not available in Apple clang yet

//void __line__(const std::string& message, const std::source_location& location) {
//  std::cout << "--- " << location.line();
//  std::cout << " in " << std::filesystem::path(location.file_name()).filename();
//  if (!message.empty()) {
//    std::cout << " : " << message;
//  }
//  std::cout << " ---" << std::endl;
//}

void __line_impl__(const std::string& message, std::string filename, int line_no) {
  std::cout << "--- " << line_no;
  if (!filename.empty()) {

    std::cout << " in " << std::filesystem::path(filename).filename().string();
  }
  if (!message.empty()) {
    std::cout << " : " << message;
  }
  std::cout << " ---" << std::endl;
}


void __pp__() {
  std::cout << std::endl;
}

void __pp_impl__(const std::vector<bool>& vs) {
  std::cout << "[";
  bool isFirst = true;
  for(auto v : vs) {
    if (isFirst) {
      isFirst = false;
    } else {
      std::cout << ",";
    }
    __pp_impl__(v);
  }
  std::cout << "]";
}

void __pp_no_endl__() {
  // do nothing
}


