//
// Created by Tsz-Chiu Au on 6/22/22.
//

#include "stl.h"


bool advance_vector_counter(std::vector<int>& counter, const std::vector<int>& counter_size) {
  unsigned int i = 0;
  while(i < counter.size()) {
    counter[i]++;
    if (counter[i] < counter_size[i]) return true;
    counter[i]=0;
    i++;
  }
  return false;
}
