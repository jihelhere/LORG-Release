// -*- mode: c++ -*-

#ifndef _ANNOTATEDRULE_H_
#define _ANNOTATEDRULE_H_

#include "MetaProduction.h"
#include <set>
#include <map>
#include <vector>

typedef std::map<int,std::vector<int> > Merge_map;
typedef std::map<int,std::set<int> > Merge_lookup;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#pragma clang diagnostic pop

class AnnotatedRule : public MetaProduction
{
public:
  AnnotatedRule() : MetaProduction() {};
  AnnotatedRule(int l) : MetaProduction(l) {};

  virtual ~AnnotatedRule() {};
  virtual bool is_lexical() const=0;
  virtual bool is_unary() const=0;
  virtual bool is_binary() const=0;

  static bool continue_with_cell (int split_size, int index,int, const std::map<int,std::set<int> > & );


 private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/)
    {
      ar & lhs;
    }

};


//why is this a method for Rule ?
// it doesn.t modify a rule
// it doesn't read anything from a rule
// TODO: move it to AnnotatedRule.cpp ?
inline
bool AnnotatedRule::continue_with_cell (int split_size, int index, int label, const std::map<int,std::set<int> >& full_set)
{

  //if this label contains a latent annotation set that will be merged
  //if this latent annotation index is part of an annotation set that will be merged....
  //but if the annotation is not the begin index of the annotation set
  return !((index % split_size) != 0 && full_set.count(label) &&  full_set.find(label)->second.count(index));
}

#endif /* _ANNOTATEDRULE_H_ */
