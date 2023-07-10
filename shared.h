#ifndef MULTIPLAN_SHARED_H
#define MULTIPLAN_SHARED_H

#include <iostream>
#include <memory>
#include <mutex>
#include <random>
#include <cassert>


/* --------------------------------------------------------------------------------------------------
 * A singleton of global variables.
 * -------------------------------------------------------------------------------------------------- */

class Shared final {
public:

  ~Shared() = default;

  Shared(const Shared &) = delete;
  Shared(Shared &&) = delete;
  Shared &operator=(const Shared &) = delete;
  Shared &operator=(Shared &&) = delete;

  static void init(std::random_device::result_type rand_seed);

  static Shared& getInstance() {
    return *Shared::instance;
  }

  [[nodiscard]] std::random_device::result_type getRandSeed() const { return rand_seed; }
  std::mt19937& getRng() { return rng; }

private:

  Shared() = default;

  void setRandSeed(std::random_device::result_type rand_seed);

private:

  static std::unique_ptr<Shared> instance;
  static std::once_flag flag;

  std::random_device::result_type rand_seed = 0;
  std::mt19937 rng;
};


#endif //MULTIPLAN_SHARED_H
