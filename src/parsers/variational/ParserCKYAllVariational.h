// -*- mode: c++ -*-
#ifndef _PARSERCKYALLVARIATIONAL_H_
#define _PARSERCKYALLVARIATIONAL_H_

#include "parsers/maxrule/ParserCKYAllMaxRule.h"
#include "VariationalTypes.h"



class ParserCKYAllVariational : public ParserCKYAllMaxRule<VariationalTypes>
{
 private:
  unsigned k;

 public:
  ParserCKYAllVariational(std::vector<AGrammar*>& cgs,
                          const std::vector<double>& p, double b_t,
                          const annot_descendants_type& annot_descendants_,
                          bool accurate_, unsigned min_beam, int stubborn, unsigned k_);

  ~ParserCKYAllVariational() {};

  void extract_solution();

 private:
  void initialise_candidates();

  void extend_all_derivations();
};

#endif /* _PARSERCKYALLVARIATIONAL_H_ */
