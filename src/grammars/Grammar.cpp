#include "Grammar.h"

#include <iostream>


template <typename Bin, typename Un, typename Lex>
std::ostream& operator<<(std::ostream& out, const Grammar<Bin,Un,Lex> & gram)
{
  for(auto& rule:  gram.binary_rules) {out << rule << std::endl;}
  for(auto& rule:   gram.unary_rules) {out << rule << std::endl;}
  for(auto& rule: gram.lexical_rules) {out << rule << std::endl;}
}
