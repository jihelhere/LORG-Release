// -*- mode: c++ -*-
#pragma once

#include <boost/functional/hash.hpp>

namespace std
{
template <class U, class V>
class hash<pair<U,V> > {
 public:
  // from boost
  inline size_t operator()(const pair<U,V> & p) const {
    size_t seed = 0;
    seed ^= std::hash<U>()(p.first) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    seed ^= std::hash<V>()(p.second) + 0x9e3779b9 + (seed<<6) + (seed>>2);

    return seed;
  }
};

template <typename... Args>
class hash< tuple<Args...> >
{
 public:
  inline size_t operator()(const tuple<Args...>& t) const
  {

    // temporary fix
    return boost::hash_value(t);

    // following code doesn't compile!!
    // size_t seed = 0;

    // for(size_t i = 0; i < std::tuple_size<decltype(t)>(t); ++i)
    // {
    //   seed ^= std::hash<decltype(std::get<i>(t))>()(std::get<i>(t)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    // }

    // return seed;
  }
};

};
