// -*- mode: c++ -*-
#ifndef _MAXRULEPROBABILITYKB_CPP_
#define _MAXRULEPROBABILITYKB_CPP_

#include "MaxRuleProbabilityKB.h"


unsigned MaxRuleProbabilityKB::size = 0;


std::ostream & MaxRuleProbabilityKB::operator>> (std::ostream & out) const
{
  for(auto& cand: candidates) { out << "cand:" << cand.probability << " "; }
  return out;
}

template<>
ParserCKYAllFactory::MaxParsing_Calculation
MaxRuleTreeLogProbaComputer<MaxRuleProbabilityKB>::calculation = ParserCKYAllFactory::Product;


#endif /* _MAXRULEPROBABILITYKB_CPP_ */
