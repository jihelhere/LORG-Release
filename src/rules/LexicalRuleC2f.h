// -*- mode: c++ -*-

#ifndef _LEXICALRULEC2F_H_
#define _LEXICALRULEC2F_H_

#include "LexicalRule.h"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#include "C2f.h"
#pragma clang diagnostic pop
//#include <cmath>

typedef C2f<LexicalRule> LexicalRuleC2f;

// class LexicalRuleC2f : public C2f<LexicalRule>
// {
// public:
//   LexicalRuleC2f(const LexicalRule& r) : C2f<LexicalRule>(r) {}

//   void set_logmode()
//   {
//     for(unsigned i = 0; i < probabilities.size(); ++i)
//       probabilities[i] = std::log(probabilities[i]);
//     logmode = true;
//   }
// };

#endif /* _LEXICALRULEC2F_H_ */
