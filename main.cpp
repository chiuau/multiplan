#include <iostream>

#include "shared.h"
#include "formation_graph.h"

// ********************************************************************************
//   The Main Function
// ********************************************************************************

int main(int argc, char** argv) {

  std::random_device::result_type rand_seed = 0;
  Shared::init(rand_seed);

  pure_formation_graph_expr_1();

  return 0;
}
