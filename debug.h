//
// Created by Tsz-Chiu Au on 6/22/22.
//

#ifndef UTIL_DEBUG_H
#define UTIL_DEBUG_H

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

// #define NDEBUG


/* --------------------------------------------------------------------------------------------------
 * __line__() - print the line number for debugging
 *
 * Usage: __line__()
 *        __line__("Here")
 *
 * See: https://stackoverflow.com/questions/597078/file-line-and-function-usage-in-c
 * see: https://stackoverflow.com/questions/3046889/optional-parameters-with-c-macros
 * see: https://gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
 * -------------------------------------------------------------------------------------------------- */

// void __line__(const std::string& message = "", const std::source_location& location = std::source_location::current());

void __line_impl__(const std::string& message, std::string filename, int line_no);
#define __line_0_args__() __line_impl__("", __FILE__, __LINE__)
#define __line_1_args__(message) __line_impl__(message, __FILE__, __LINE__)
#define __line_get_2nd_arg__(arg1, arg2, ...) arg2
#define __line_macro_chooser__(...) __line_get_2nd_arg__(__VA_ARGS__ __VA_OPT__(,) __line_1_args__, __line_0_args__)
#define __line__(...) __line_macro_chooser__(__VA_ARGS__)(__VA_ARGS__)


/* --------------------------------------------------------------------------------------------------
 * __pp__() - print something for debugging
 *
 * Usage: __pp__(v)
 *
 * -------------------------------------------------------------------------------------------------- */

template<typename T1, typename T2>
void __pp_impl__(const std::pair<T1, T2>& p);

template<typename T>
void __pp_impl__(const std::vector<T>& vs);

template<typename T>
void __pp_impl__(const std::list<T>& vs);

template<typename T>
void __pp_impl__(const std::set<T>& vs);

template<typename K, typename T>
void __pp_impl__(const std::unordered_map<K,T>& map);


template<typename T>
void __pp_impl__(T v) {
  std::cout << v;
}

template<typename T1, typename T2>
void __pp_impl__(const std::pair<T1, T2>& p) {
  std::cout << "(";
  __pp_impl__(p.first);
  std::cout << ",";
  __pp_impl__(p.second);
  std::cout << ")";
}

template<typename T>
void __pp_impl__(const std::vector<T>& vs) {
  std::cout << "[";
  bool isFirst = true;
  for(auto& v : vs) {
    if (isFirst) {
      isFirst = false;
    } else {
      std::cout << ",";
    }
    __pp_impl__(v);
  }
  std::cout << "]";
}

template<typename T>
void __pp_impl__(const std::list<T>& vs) {
  std::cout << "<";
  bool isFirst = true;
  for(auto& v : vs) {
    if (isFirst) {
      isFirst = false;
    } else {
      std::cout << ",";
    }
    __pp_impl__(v);
  }
  std::cout << ">";
}

void __pp_impl__(const std::vector<bool>& vs);

template<typename T>
void __pp_impl__(const std::set<T>& vs) {
  std::cout << "{";
  bool isFirst = true;
  for(auto& v : vs) {
    if (isFirst) {
      isFirst = false;
    } else {
      std::cout << ",";
    }
    __pp_impl__(v);
  }
  std::cout << "}";
}

template<typename K, typename T>
void __pp_impl__(const std::unordered_map<K,T>& map) {
  std::cout << "{";
  bool isFirst = true;
  for(auto& [k, v] : map) {
    if (isFirst) {
      isFirst = false;
    } else {
      std::cout << ",";
    }
    __pp_impl__(k);
    std::cout << "->";
    __pp_impl__(v);
  }
  std::cout << "}";
}



void __pp__();

template<typename T, typename... TS>
void __pp__(T&& v, TS... vs) {
  __pp_impl__(std::forward<T>(v));
  __pp__(vs...);
}

void __pp_no_endl__();

template<typename T, typename... TS>
void __pp_no_endl__(T&& v, TS... vs) {
  __pp_impl__(std::forward<T>(v));
  __pp_no_endl__(vs...);
}


// see https://www.scs.stanford.edu/~dm/blog/va-opt.html

#define PARENS ()
#define EXPAND(...) EXPAND4(EXPAND4(EXPAND4(EXPAND4(__VA_ARGS__))))
#define EXPAND4(...) EXPAND3(EXPAND3(EXPAND3(EXPAND3(__VA_ARGS__))))
#define EXPAND3(...) EXPAND2(EXPAND2(EXPAND2(EXPAND2(__VA_ARGS__))))
#define EXPAND2(...) EXPAND1(EXPAND1(EXPAND1(EXPAND1(__VA_ARGS__))))
#define EXPAND1(...) __VA_ARGS__
#define FOR_EACH(macro, ...)                                    \
  __VA_OPT__(EXPAND(FOR_EACH_HELPER(macro, __VA_ARGS__)))
#define FOR_EACH_HELPER(macro, a1, ...)                         \
  macro(a1)                                                     \
  __VA_OPT__(FOR_EACH_AGAIN PARENS (macro, __VA_ARGS__))
#define FOR_EACH_AGAIN() FOR_EACH_HELPER

#define __vv_impl__(...) __pp_no_endl__(#__VA_ARGS__, " = ", __VA_ARGS__, "  ");

#define __vv__(...) FOR_EACH(__vv_impl__, __VA_ARGS__); __pp__();



//* See https://stackoverflow.com/questions/40626433/c-how-to-specialize-a-template-using-vectort
//template<typename T>
//struct is_vector {
//  static constexpr bool value = false;
//};
//
//template<template<typename...> class C, typename U>
//struct is_vector<C<U>> {
//  static constexpr bool value = std::is_same<C<U>,std::vector<U>>::value;
//};



/* --------------------------------------------------------------------------------------------------
 * type_name(v) - Show the type information of a variable v
 *
 * Usage: std::cout << type_name<decltype(v)>() << std::endl;
 *
 * See: https://stackoverflow.com/questions/81870/is-it-possible-to-print-a-variables-type-in-standard-c
 * -------------------------------------------------------------------------------------------------- */

template <class T>
constexpr std::string_view type_name() {
  using namespace std;
#ifdef __clang__
  string_view p = __PRETTY_FUNCTION__;
  return string_view(p.data() + 34, p.size() - 34 - 1);
#elif defined(__GNUC__)
  string_view p = __PRETTY_FUNCTION__;
#  if __cplusplus < 201402
  return string_view(p.data() + 36, p.size() - 36 - 1);
#  else
  return string_view(p.data() + 49, p.find(';', 49) - 49);
#  endif
#elif defined(_MSC_VER)
  string_view p = __FUNCSIG__;
    return string_view(p.data() + 84, p.size() - 84 - 7);
#endif
}


#endif //UTIL_DEBUG_H
