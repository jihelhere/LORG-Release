#include "compact_binary_rules.h"

#include <algorithm>
#include <stdexcept>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#include <boost/phoenix/bind/bind_member_function.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/phoenix/stl.hpp>
#pragma clang diagnostic pop
using namespace compact_binary_rules;

using namespace boost::phoenix;


template<class BinaryRule, typename info>
info transform(const BinaryRule&)
{
  throw std::runtime_error("compact_binary_rules::transform should be specialised");
}


template<class BinaryRule, typename info>
void build_vector_rhs1(const std::vector<const BinaryRule*>& pre, vector_rhs1<info>& result )
{
  result.rules.resize(pre.size());

  std::transform(pre.begin(), pre.end(), result.rules.begin(),
                 [](info e)
                 {return transform<BinaryRule,info>(*e);}
                 );

  result._begin = result.rules.begin();
  result._end = result.rules.end();
  result.rhs1 = pre.front()->get_rhs1();

  result.rules.shrink_to_fit();
}


  template<class BinaryRule, typename info>
    void build_vector_rhs0(const std::vector<std::vector<const BinaryRule*> >& pre,
                           vector_rhs0<info>& result)
{
  result.vrhs1.resize(pre.size());

  for(unsigned i = 0; i < pre.size(); ++i)
    build_vector_rhs1<BinaryRule,info>(pre[i],result.vrhs1[i]);

  result._begin = result.vrhs1.begin();
  result._end = result.vrhs1.end();
  result.rhs0 = pre.front().front()->get_rhs0();

  result.vrhs1.shrink_to_fit();
}


  // the first argument is just a vector of vector of vector of BinaryRule pointers
  template<class BinaryRule, typename info>
  void build_vector_brules(const std::vector<std::vector<std::vector<const BinaryRule*> > >& pre,
			     vector_brules<info>& result)
{
  result.vrhs0.resize(pre.size());

  for(unsigned i = 0; i < pre.size(); ++i)
    build_vector_rhs0<BinaryRule,info>(pre[i],result.vrhs0[i]);

  result._begin = result.vrhs0.begin();
  result._end = result.vrhs0.end();

  result.vrhs0.shrink_to_fit();
}

template<typename info>
template<class BinaryRule>
vector_brules<info> *
vector_brules<info>::convert(const std::vector<BinaryRule >& binary_rules)
{
  typedef std::vector<BinaryRule > vbr;
  typedef std::vector<const BinaryRule* > vbrp;
  typedef std::vector<vbrp > vvbrp;
  typedef std::vector<vvbrp> vvvbrp;

  //partition brules according to their rhs0
  vvbrp brulesranked;

  for(const auto& br : binary_rules)
  //   typename vbr::const_iterator brules_itr = binary_rules.begin();
  // brules_itr != binary_rules.end(); ++brules_itr)
  {

    int rhs0 = br.get_rhs0();

    auto br_itr = std::find_if(brulesranked.begin(),brulesranked.end(),
                               bind(&BinaryRule::get_rhs0,front(arg_names::arg1)) == rhs0
                               );

    if(br_itr != brulesranked.end())
      br_itr->push_back(&br);
    else
      brulesranked.push_back(vbrp(1,&br));
  }

  vvvbrp prebrulesrankedranked;
  prebrulesrankedranked.resize(brulesranked.size());

  //in each partition, repartition brules according to their rhs1
  for(unsigned i = 0; i < brulesranked.size();++i) {
    vbrp same_rhs0_vect = brulesranked[i];
    for(const auto& r : same_rhs0_vect)
    {
      int rhs1 = r->get_rhs1();

      auto itr = std::find_if(prebrulesrankedranked[i].begin(),
                              prebrulesrankedranked[i].end(),
                              bind(&BinaryRule::get_rhs1,front(arg_names::arg1)) == rhs1
                              );

      if(itr != prebrulesrankedranked[i].end())
	itr->push_back(r);
      else
	prebrulesrankedranked[i].push_back(vbrp(1,r));
    }
  }


  vector_brules<info> * res = new vector_brules<info>();

  build_vector_brules<BinaryRule,info>(prebrulesrankedranked,*res);

  return res;
}
