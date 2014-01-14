// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVARONEBEST_CPP_
#define _PARSERCKYALLMAXVARONEBEST_CPP_

#include "ParserCKYAllMaxRule1B.h"
#include "parsers/ParserCKYAll.hpp"

ParserCKYAllMaxRule1B::ParserCKYAllMaxRule1B(ParserCKYAllFactory::MaxParsing_Calculation c,
                                             std::vector<AGrammar*>& cgs,
                                             const std::vector<double>& p, double b_t,
                                             const annot_descendants_type& annot_descendants_,
                                             bool accurate_, unsigned min_beam, int stubborn)
    : ParserCKYAllMaxRule<MaxRule1BTypes>(c, cgs, p, b_t, annot_descendants_, accurate_, min_beam, stubborn)
{
  // this is not in the super class because maxn parsing uses a
  //different mapping
  //create the coarse-to-fine map
  this->create_coarse_to_fine_mapping(this->grammars);

  Edge::set_unary_chains(this->grammars[this->grammars.size() - 1]->get_unary_decoding_paths());
}

void ParserCKYAllMaxRule1B::extract_solution()
{
  //   std::cout << *chart << std::endl;
  compute_inside_outside_probabilities();
  //   std::cout << *chart << std::endl;
  calculate_chart_specific_rule_probabilities_and_best_edge();
  //   std:e:cout << *chart << std::endl;

  // if (!chart->get_root().exists_edge(SymbolTable::instance_nt().get_label_id(LorgConstants::tree_root_name)))
  //   std::cout << "no axiom at root" << std::endl;
}

void ParserCKYAllMaxRule1B::calculate_chart_specific_rule_probabilities_and_best_edge()
{
  //std::cout << "calculate_chart_specific_rule_probabilities_and_best_edge: " << get_sentence_probability() << std::endl;

  double sentence_log_probability = std::log(get_sentence_probability());
  this->set_log_normalisation_factor(sentence_log_probability);
  ParserCKYAllMaxRule::calculate_maxrule_probabilities(log_normalisation_factor);
}

void ParserCKYAllMaxRule1B::simple_extract_solution()
{
  ParserCKYAllMaxRule::calculate_maxrule_probabilities(log_normalisation_factor);
}


void ParserCKYAllMaxRule1B::set_log_normalisation_factor(double lnf)
{
  log_normalisation_factor = lnf;
}

#endif /* _PARSERCKYALLMAXVARONEBEST_H_ */
