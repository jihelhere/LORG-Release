// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVARKB_H_
#define _PARSERCKYALLMAXVARKB_H_

#include "MaxRuleProbabilityKB.h"
#include "ParserCKYAllMaxRule.h"


class ParserCKYAllMaxRuleKB : public ParserCKYAllMaxRule<MaxRuleKBTypes>
{
 private:
  unsigned k;

 public:
  ParserCKYAllMaxRuleKB(ParserCKYAllFactory::MaxParsing_Calculation c,
                        const std::vector<AGrammar*>& cgs,
                        const std::vector<double>& p, double b_t,
                        const annot_descendants_type& annot_descendants_,
                        bool accurate_, unsigned min_beam, int stubborn, unsigned k_);

  ~ParserCKYAllMaxRuleKB() {};

  void extract_solution();

  void simple_extract_solution();

 private:
  void initialise_candidates();

  void extend_all_derivations();

  double log_normalisation_factor;
  void set_log_normalisation_factor(double lnf) {log_normalisation_factor = lnf;};
};

#endif /* _PARSERCKYALLMAXVARKB_H_ */
