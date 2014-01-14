// // -*- mode: c++ -*-
#ifndef _MAXRULEPROBABILITY_CPP_
#define _MAXRULEPROBABILITY_CPP_

#include "MaxRuleProbability1B.h"

#include <numeric>

template<>
ParserCKYAllFactory::MaxParsing_Calculation
MaxRuleTreeLogProbaComputer<MaxRuleProbability1B>::calculation = ParserCKYAllFactory::Product;

#endif /* _MAXRULEPROBABILITY_H_ */
