// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVARMULTIPLE_CPP_
#define _PARSERCKYALLMAXVARMULTIPLE_CPP_

#include "parsers/ParserCKYAll.hpp"
#include "ParserCKYAllMaxRuleMultiple.h"



ParserCKYAllMaxRuleMultiple::ParserCKYAllMaxRuleMultiple(ParserCKYAllFactory::MaxParsing_Calculation c,
                                                         std::vector<AGrammar*>& cgs,
                                                         const std::vector<double>& p, double b_t,
                                                         const std::vector< std::vector<AGrammar*> >& fgs,
                                                         const std::vector< annot_descendants_type >& all_annot_descendants_,
                                                         bool accurate_, unsigned min_beam, int stubborn, unsigned k_)
: ParserCKYAllMaxRule<MaxRuleMultipleTypes>(c, cgs, p, b_t, all_annot_descendants_[0], accurate_, min_beam, stubborn),
    fine_grammars(fgs), all_annot_descendants(all_annot_descendants_), nb_grammars(fgs.size() + 1), k(k_)
{

  // create a mapping of all grammars
  std::vector<AGrammar*> all_grammars(grammars);

  for (unsigned i = 0; i < fine_grammars.size(); ++i)
    all_grammars.insert(all_grammars.end(), fine_grammars[i].begin(), fine_grammars[i].end());

  create_coarse_to_fine_mapping(all_grammars);


  // create another mapping between final grammars starting with the last one
  // that's how MaxRuleProbability class knows where to find rules
  std::vector<AGrammar*>& lasts = fine_grammars.back();
  std::vector<AGrammar*> maxn_mapping(1,lasts.back());
  maxn_mapping.push_back(grammars.back());
  //  create_coarse_to_fine_mapping(maxn_mapping);
  for (unsigned i = 0; i < fine_grammars.size(); ++i) {
    maxn_mapping.push_back(fine_grammars[i][fine_grammars[i].size()-1]);
  }

  create_coarse_to_fine_mapping(maxn_mapping);

  //TODO calculate this properly for multiple grammars
  Edge::set_unary_chains(grammars.back()->get_unary_decoding_paths());
}


ParserCKYAllMaxRuleMultiple::~ParserCKYAllMaxRuleMultiple()
{
  for(unsigned i = 0; i < fine_grammars.size(); ++i)
    for(unsigned j = 0; j < fine_grammars[i].size(); ++j)
      delete fine_grammars[i][j];
}


void ParserCKYAllMaxRuleMultiple::change_rules_reset()
{
  this->chart->opencells_apply(
      [](Cell& cell)
      {
        // 0 means c2f
        // 1 means multiple grammar decoding
        cell.change_rules_resize(1,0);
      });
}


void ParserCKYAllMaxRuleMultiple::change_rules_load_backup(unsigned backup_idx, unsigned size) const
{
  //  std::cout << "change_rules_load_backup" << std::endl;
  function<void(Edge&)> replace_rules = std::bind(&Edge::replace_rule_probabilities, std::placeholders::_1, size);
  function<void(Edge&)> replace_annotations =
      [backup_idx](Edge& e)
      {
        e.get_annotations() = e.get_prob_model().get_annotations_backup()[backup_idx];
      };

  chart->opencells_apply(
    [&replace_rules, &replace_annotations](Cell&cell)
    {
      cell.apply_on_edges(replace_rules, replace_annotations);
    });
}

void ParserCKYAllMaxRuleMultiple::modify_backup(unsigned backup_idx) const
{
  function<void(Edge&)> modify = [backup_idx](Edge& e){e.get_prob_model().get_annotations_backup()[backup_idx] = e.get_annotations();};
  chart->opencells_apply([&modify](Cell&cell){cell.apply_on_edges(modify);});
}



void ParserCKYAllMaxRuleMultiple::precompute_all_backups()
{
  backup_annotations();

  for(unsigned i = 0; i < fine_grammars.size(); ++i) {

    //    std::cout << "fine grammar " << i << std::endl;

    //    std::cout << "changing rules" << std::endl;
    // change grammar rules in the chart and  pick the new baseline
    change_rules_reset();
    //    std::cout << "changed rules" << std::endl;

    //    std::cout << "changing annotation history" << std::endl;
    // update annotation history
    annot_descendants = all_annot_descendants[i+1];
    //    std::cout << "changed annotation history" << std::endl;

    //    std::cout << "before beam_c2f " << std::endl;
    beam_c2f(fine_grammars[i], annot_descendants);
    //    std::cout << "after beam_c2f " << std::endl;

    backup_annotations();
  }
}

void ParserCKYAllMaxRuleMultiple::multiple_inside_outside_specific()
{
  static int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  MaxRuleProbabilityMultiple::reset_log_normalisation_factor();

  for(unsigned i = 0; i < fine_grammars.size() + 1; ++i) {

    //    std::cout << "computation " << i << std::endl;

    //assume that all grammars have the same size
    unsigned distance = i==0 ? 0 : 1;
    change_rules_load_backup(i, distance);

    compute_inside_probabilities();
    //    std::cout << "sentence_prob: " << std::log(get_sentence_probability()) << std::endl;


    if(!chart->get_root().is_closed() && chart->get_root().exists_edge(start_symbol)) {
      chart->get_root().get_edge(start_symbol).get_annotations().reset_outside_probabilities(1.0, false);
      compute_outside_probabilities();

      MaxRuleProbabilityMultiple::set_log_normalisation_factor(std::log(get_sentence_probability()));

    }

    modify_backup(i);

  }
}


void ParserCKYAllMaxRuleMultiple::extract_solution()
{
  static int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  // first we backup all annotations for all grammars

  //std::cerr << "precompute_all_backups" << std::endl;
  precompute_all_backups();

  // actual computation
  // inside/outside + max rule score
  // all in one function
  //std::cerr << "multiple_inside_outside_specific" << std::endl;
  multiple_inside_outside_specific();


  if(!chart->get_root().is_closed() && chart->get_root().exists_edge(start_symbol)) {
    //std::cerr << "initialise_candidates" << std::endl;
    initialise_candidates();

    //std::cerr << "extend_all_derivations" << std::endl;
    extend_all_derivations();
    //std::cerr << "after extend_all_derivations" << std::endl;
  }


  //reset to first grammar for next parse
  annot_descendants = all_annot_descendants[0];
}

void ParserCKYAllMaxRuleMultiple::simple_extract_solution()
{
  static int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  if(!chart->get_root().is_closed() && chart->get_root().exists_edge(start_symbol)) {
    //    std::cerr << "calculate_best_edge" << std::endl;

    initialise_candidates();
    extend_all_derivations();
  }


  //reset to first grammar for next parse
  annot_descendants = all_annot_descendants[0];

  //throw std::runtime_error("not implemented yet");
}


void ParserCKYAllMaxRuleMultiple::calculate_maxrule_probabilities()
{

  //unsigned sent_size = chart->get_size();

  ParserCKYAllMaxRule::calculate_maxrule_probabilities();
}



void ParserCKYAllMaxRuleMultiple::calculate_best_edge()
{
  chart->opencells_apply_bottom_up( [](Cell&cell)
  {
    cell.apply_on_edges(
        &ProbaModel::init,
        &ProbaModel::pick_best_lexical,
        &ProbaModel::pick_best_binary );
    cell.apply_on_edges( &ProbaModel::pick_best_unary,
                         &ProbaModel::pick_best );
  }  );

  // std::cout << "init" << std::endl;
  // chart->opencells_apply_bottom_up([](Cell&cell) { cell.apply_on_edges(&ProbaModel::init);});
  // std::cout << "lexical" << std::endl;
  // chart->opencells_apply_bottom_up([](Cell&cell) { cell.apply_on_edges(&ProbaModel::pick_best_lexical);});
  // std::cout << "binary" << std::endl;
  // chart->opencells_apply_bottom_up([](Cell&cell) { cell.apply_on_edges(&ProbaModel::pick_best_binary);});
  // std::cout << "unary" << std::endl;
  // chart->opencells_apply_bottom_up([](Cell&cell) { cell.apply_on_edges(&ProbaModel::pick_best_unary);});
  // std::cout << "best" << std::endl;
  // chart->opencells_apply_bottom_up( [](Cell&cell){ cell.apply_on_edges(&ProbaModel::pick_best);});
}

void ParserCKYAllMaxRuleMultiple::backup_annotations() const
{
  chart->opencells_apply([](Cell&cell){
      cell.apply_on_edges(&MaxRuleProbabilityMultiple::backup_annotations);}
  );
}


void ParserCKYAllMaxRuleMultiple::extend_all_derivations()
{
  static int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  Cell& root = chart->get_root();

  if (!root.exists_edge(start_symbol))
    //   //   std::cout << "no axiom at root" << std::endl;
    return;

  for (unsigned i = 2; i <= k; ++i)
    {
      //      std::cout << "before extend" << std::endl;
      chart->get_root().get_edge(start_symbol).extend_derivation(i,true);
    }
}


void ParserCKYAllMaxRuleMultiple::initialise_candidates()
{
  MaxRuleProbabilityMultiple::set_size(k);
  MaxRuleProbabilityMultiple::set_nbgrammars(nb_grammars);


  calculate_best_edge();
}


#endif /* _PARSERCKYALLMAXVARMULTIPLE_CPP_ */
