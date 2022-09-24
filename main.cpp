#include <iostream>
#include <string>
#include <random>
#include <memory>
#include <limits>
#include <chrono>

#include "shared.h"
#include "formation_graph.h"

// ********************************************************************************
//   The Main Function
// ********************************************************************************

int main(int argc, char** argv) {

  // std::random_device::result_type rand_seed = 2746795615;
  std::random_device::result_type rand_seed = 0;
  Shared::init(rand_seed, false);
  std::cout << "rand_seed = " << Shared::getInstance().getRandSeed() << std::endl;

  test_formation_graph_1();
  // test_formation_graph_2();

  return 0;
}
