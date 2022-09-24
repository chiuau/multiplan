#include <sstream>
#include <iomanip>
#include <cctype>
#include <filesystem>
#include <algorithm>

#include "shared.h"


std::unique_ptr<Shared> Shared::instance = nullptr;
std::once_flag Shared::flag;

void Shared::init(std::random_device::result_type rand_seed, bool is_show_rand_seed) {

  // cout setting
  std::cout << std::boolalpha;

  // random seek
  std::call_once(Shared::flag, [&]() {
    Shared::instance.reset(new Shared);
    Shared::instance->setRandSeed(rand_seed);
    if (is_show_rand_seed) {
      std::cout << "Random seed = " << Shared::instance->rand_seed << std::endl;
    }
  });
}

void Shared::setRandSeed(std::random_device::result_type new_rand_seed) {
  if (new_rand_seed == 0) {
    std::random_device dev;
    new_rand_seed = dev();
  }
  rand_seed = new_rand_seed;

  rng = std::mt19937{rand_seed};
}

