// -*- mode: c++ -*-
#pragma  once

#include <cassert>
#include <cmath>

// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>


namespace {
  void set_log(std::vector<std::vector<std::vector<double> > >& probabilities)
  {
    for(unsigned i = 0; i < probabilities.size(); ++i)
      for(unsigned j = 0; j < probabilities[i].size(); ++j)
	for(unsigned k = 0; k < probabilities[i][j].size(); ++k)
	  probabilities[i][j][k] = std::log(probabilities[i][j][k]);
  }

  void set_log(std::vector<std::vector<double> >& probabilities)
  {
    for(unsigned i = 0; i < probabilities.size(); ++i)
      for(unsigned j = 0; j < probabilities[i].size(); ++j)
	probabilities[i][j] = std::log(probabilities[i][j]);
  }

  void set_log(std::vector<double>& probabilities)
  {
    for(unsigned i = 0; i < probabilities.size(); ++i)
      probabilities[i] = std::log(probabilities[i]);
  }
}




template <class MyRule>
class C2f : public MyRule
{
public:
  C2f(const MyRule& r) : MyRule(r), finers(), coarser(NULL), logmode(false){};
  C2f() : MyRule(), finers(), coarser(NULL), logmode(false){};

  virtual ~C2f() {};

  // void set_finers(const std::vector< C2f*>& fs)
  // {
  //   finers = fs;
  // }

  void add_finer(C2f* f)
  {
    finers.push_back(f);
    finers.back()->coarser = this;
  }

  const std::vector<const C2f*>& get_finers() const
  {
    return finers;
  }

  const C2f * get(int index) const
  {
    if(index < 0) return this;
    assert(unsigned(index) < finers.size());
    return finers[index];
  }


  const C2f* get_coarser() const
  {
    assert(coarser);
    return coarser;
  }


  const C2f* get_coarser(unsigned i) const
  {
    const C2f* res = this;
    while(i-- && res)
      res = res->coarser;

    return res;
  }


  bool is_logmode() const {return logmode;}
  void set_logmode()
  {
    set_log(this->probabilities); logmode=true;
  };




private:
  std::vector< C2f*> finers;
  C2f* coarser;

  friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/)
    {
      ar & boost::serialization::base_object<MyRule>(*this);
      ar & finers;
      ar & coarser;
      ar & logmode;
    }

protected:
  bool logmode;

};
