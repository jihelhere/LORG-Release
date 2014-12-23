// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVARONEBEST_H_
#define _PARSERCKYALLMAXVARONEBEST_H_

#include "ParserCKYAllMaxRule.h"
#include "MaxRuleProbability1B.h"
#include "Word.h"

#include "emptystruct.h"

#include "parsers/ParserCKYAllFactory.h"


class ParserCKYAllMaxRule1B : public ParserCKYAllMaxRule<MaxRule1BTypes>
{
public:
  ParserCKYAllMaxRule1B(ParserCKYAllFactory::MaxParsing_Calculation c,
                        const std::vector<AGrammar*>& cgs,
                        const std::vector<double>& p, double b_t,
                        const annot_descendants_type& annot_descendants_,
                        bool accurate_, unsigned min_beam, int stubborn);

  ~ParserCKYAllMaxRule1B() {};

  void extract_solution();

  void simple_extract_solution();

protected:

  void change_rules_reset() const;


  /**
     \brief Calculates the chart specific rule probabilities of the packed edges in the chart
     and uses this to select the best edge (max rule parsing)
   */
  void calculate_chart_specific_rule_probabilities_and_best_edge();

 private:
  double log_normalisation_factor;

  void set_log_normalisation_factor(double lnf);


};

#endif /* _PARSERCKYALLMAXVARONEBEST_H_ */
