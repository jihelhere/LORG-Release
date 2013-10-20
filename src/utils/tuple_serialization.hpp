# pragma once




#include <boost/config.hpp>

#include <boost/serialization/utility.hpp>
#include <boost/serialization/collections_save_imp.hpp>
#include <boost/serialization/collections_load_imp.hpp>
#include <boost/serialization/split_free.hpp>

// shamelessly copied from
// https://sydius.me/2011/02/c0x-tuple-boost-serialization/

namespace boost {
namespace serialization {

template<uint N>
struct Serialize
{
  template<class Archive, typename... Args>
  static void serialize(Archive & ar, std::tuple<Args...> & t, const unsigned int version)
  {
    ar & std::get<N-1>(t);
    Serialize<N-1>::serialize(ar, t, version);
  }
};

template<>
struct Serialize<0>
{
  template<class Archive, typename... Args>
  static void serialize(Archive & /*ar*/, std::tuple<Args...> & /*t*/, const unsigned int /*version*/)
  {
  }
};

template<class Archive, typename... Args>
void serialize(Archive & ar, std::tuple<Args...> & t, const unsigned int version)
{
  Serialize<sizeof...(Args)>::serialize(ar, t, version);
}

}
}
