// -*- mode: c++ -*-
#pragma once

#include <boost/functional/hash.hpp>

class AnnotatedRule;

namespace std
{
template <class U, class V>
struct hash<pair<U,V> > {
 public:

  inline size_t operator()(const pair<U,V> & p) const {
    size_t seed = 0;
    seed ^= std::hash<U>()(p.first) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    seed ^= std::hash<V>()(p.second) + 0x9e3779b9 + (seed<<6) + (seed>>2);

    return seed;
  }
};

// from somewhere on StackOverflow

template<typename... TTypes>
struct hash<std::tuple<TTypes...>>
{
  typedef std::tuple<TTypes...> Tuple;

  template<int N>
      size_t operator()(Tuple /*value*/) const { return 0; }

  template<int N, typename THead, typename... TTail>
      size_t operator()(Tuple value) const
  {
    constexpr int Index = N - sizeof...(TTail) - 1;
    return hash<THead>()(std::get<Index>(value)) ^ operator()<N, TTail...>(value);
  }

  size_t operator()(Tuple value) const
  {
    return operator()<sizeof...(TTypes), TTypes...>(value);
  }
};
};
