//
// Created by Tsz-Chiu Au on 6/22/22.
//

#ifndef UTIL_STL_H
#define UTIL_STL_H

#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <unordered_map>
#include <memory>
#include <mutex>
#include <random>
#include <cassert>
#include <string_view>
#include <cmath>
// #include <experimental/source_location>
#include <utility>
#include <type_traits>


/* --------------------------------------------------------------------------------------------------
 * unique_cast(v) - Cast a unique_ptr to a different type
 *
 * Usage: std::unique_ptr<MyClass> obj = unique_cast<MyClass>(std::make_unique<MyOtherClass>());
 *
 * See: https://stackoverflow.com/questions/11002641/dynamic-casting-for-unique-ptr
 * -------------------------------------------------------------------------------------------------- */

template <class destinationT, typename sourceT>
std::unique_ptr<destinationT> unique_cast(std::unique_ptr<sourceT>&& source)
{
  if (!source)
    return std::unique_ptr<destinationT>();

  destinationT* dest_ptr = dynamic_cast<destinationT*>(source.get());
  if(dest_ptr) {
    source.release();
    return std::unique_ptr<destinationT>(dest_ptr);
  } else {
    return std::unique_ptr<destinationT>();
  }
}


/* --------------------------------------------------------------------------------------------------
 * A hash function for integer vectors.
 *
 * See: https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector/20511429
 * See: https://en.cppreference.com/w/cpp/utility/hash
 * See: https://en.cppreference.com/w/cpp/container/vector_bool/hash
 * See: https://www.boost.org/doc/libs/1_76_0/doc/html/hash/reference.html#header.boost.container_hash.hash_hpp
 * -------------------------------------------------------------------------------------------------- */

namespace std
{
  template<> struct hash<std::vector<int>> {
    std::size_t operator()(std::vector<int> const& vec) const {
      std::size_t seed = vec.size();
      for(auto& i : vec) {
        seed ^= static_cast<uint32_t>(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }
      return seed;
    }
  };

  // see https://stackoverflow.com/questions/32685540/why-cant-i-compile-an-unordered-map-with-a-pair-as-key
  template<class T1,class T2>
  struct hash<std::pair<T1,T2>> {
    std::size_t operator()(const std::pair<T1, T2> &p) const {
      auto h1 = std::hash<T1>{}(p.first);
      auto h2 = std::hash<T2>{}(p.second);
      // Mainly for demonstration purposes, i.e. works but is overly simple
      // In the real world, use sth. like boost.hash_combine
      return (h1 * 65536) ^ h2;
    }
  };

}


/* --------------------------------------------------------------------------------------------------
 * vector-based counter
 * -------------------------------------------------------------------------------------------------- */

bool advance_vector_counter(std::vector<int>& counter, const std::vector<int>& counter_size);


/* --------------------------------------------------------------------------------------------------
 * remove_front() - remove some elements from the front of a vector and store the result in a new vector.
 *
 * Usage: auto v2 = remove_front(v1, 3);
 *
 * -------------------------------------------------------------------------------------------------- */

template<typename T>
std::vector<T> remove_front(const std::vector<T>& vs, int n = 1) {
  assert(n >= 1);
  std::vector<T> result;
  for(int i=n; i<vs.size(); i++) {
    result.push_back(vs[i]);
  }
  return result;
}


#endif //UTIL_STL_H
