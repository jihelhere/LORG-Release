// -*- mode: c++ -*-
#ifndef _PARSERCKYALLFACTORY_H_
#define _PARSERCKYALLFACTORY_H_

#include "utils/ConfigTable.h"
#include "ParserCKYAll.h"


namespace ParserCKYAllFactory {
  enum Parsing_Algorithm {MaxRule, Viterbi, MaxN, KMaxRule, MinDiv, Variational};
  enum MaxParsing_Calculation {Product, Sum, ProdSum};
  std::vector<ParserCKYAll *> create_parser(ConfigTable& config);
  Parsing_Algorithm string_to_pa(const std::string& s);
  MaxParsing_Calculation string_to_mpc(const std::string& s);
}

#endif /* _PARSERCKYALLFACTORY_H_ */
