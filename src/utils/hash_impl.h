// -*- mode: c++ -*-
#pragma once

#include <boost/functional/hash.hpp>

class AnnotatedRule;

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

// template<>
// class hash<const AnnotatedRule*>
// {
//  public:
//   inline size_t operator()(const AnnotatedRule* ar) const
//   {
//     return 0;
//   }
// };


template<unsigned int N>
struct my_hash
{

  template< typename... Args>
  inline size_t operator()(const std::tuple<Args...>& t) const
  {
    size_t seed = my_hash<N-1>()(t);

    auto n = std::get<N-1>(t);
    seed ^=  std::hash<decltype(n)>()(n) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    return seed;
  }
};

template<>
struct my_hash<1>
{

  template<typename... Args>
  inline size_t operator()(const std::tuple<Args...>& t) const
  {
    auto n = std::get<0>(t);
    return std::hash<decltype(n)>()(n) + 0x9e3779b9;
    //return 0;
  }
};


template <typename... Args>
class hash< tuple<Args...> >
{
 public:


  inline size_t operator()(const tuple<Args...>& t) const
  {

    // temporary fix
    //return boost::hash_value(t);

    // following code doesn't compile!!
    size_t seed = 0;

    seed = my_hash<sizeof...(Args)>()(t);


    // for(size_t i = 0; i < std::tuple_size<decltype(t)>(t); ++i)
    // {
    //   seed ^= std::hash<decltype(std::get<i>(t))>()(std::get<i>(t)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    // }

    return seed;
  }
};

};
