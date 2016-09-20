#include "ParserCKYBest.h"
#include "SimpleChartCKY.hpp"

ParserCKYBest::~ParserCKYBest()
{
}

ParserCKYBest::ParserCKYBest(Grammar<Rule,Rule,Rule>& g):
  ParserCKY< Grammar<Rule,Rule,Rule> >(g)
{}

void ParserCKYBest::parse(Chart& chart) const
{
  bool isroot = chart.get_size() == 1;
  for(unsigned i = 0; i < chart.get_size(); ++i) {
    //    std::cout << "initialising position: " << i << std::endl;
    add_unary_init(chart.access(i,i),isroot);
    //    std::cout << chart.access(i,i) << std::endl;
  }

  process_internal_rules(chart);
}

inline
void ParserCKYBest::get_candidates(const Cell& left_cell,
				   const Cell& right_cell,
				   Cell& result_cell) const
{
  Edge current_candidate;

  // left and right edges can come from various
  // heights in unary chains
  // TODO: moe efficient implementation -> we do complete scans of NTS several times
  // even in the case the NT is absent at a lower height

  for(unsigned lh = 0; lh < Cell::unary_length; ++lh)
    for(unsigned rh = 0; rh < Cell::unary_length; ++rh)

  //iterating through all the rules P -> L R, indexed by L
  for(const auto& same_rhs0 : brules)
  {
    // is L present in  left_cell ?
    if(left_cell.exists_edge(lh, same_rhs0.rhs0)) {
      const Edge& left_edge = left_cell.at(lh, same_rhs0.rhs0);
      current_candidate.set_left_child(&left_edge);

      //iterating through all the rules P -> L R, indexed by R, L fixed
      for (const auto& same_rhs1 : same_rhs0)
      {

	// is R present in right_cell ?
	if(right_cell.exists_edge(rh, same_rhs1.rhs1)) {

	  const Edge& right_edge = right_cell.at(rh, same_rhs1.rhs1);
	  current_candidate.set_right_child(&right_edge);

	  double prob1 = left_edge.get_probability() + right_edge.get_probability();

	  //iterating through all the rules P -> L R, indexed by P, R and L fixed

	  for(const auto& b : same_rhs1)
          {
	    current_candidate.set_lhs(b->get_lhs());
	    current_candidate.set_probability(prob1 + b->get_probability());

	    //	    std::cout << *(*bitr) << std::endl;

            // always set at height 0
	    (void) result_cell.process_candidate(0, current_candidate);
	  }
	}
      }
    }
  }
}


void ParserCKYBest::process_internal_rules(Chart& chart) const
{
  unsigned sent_size=chart.get_size();
  for (unsigned span = 2; span <= sent_size; ++span) {
    unsigned end_of_begin=sent_size-span;
    for (unsigned begin=0; begin <= end_of_begin; ++begin) {
    	unsigned end = begin + span -1;

      //      std::cout << "begin: " << begin << ", end: " << end << std::endl;

    	Cell& result_cell = chart.access(begin,end);

    	if(!result_cell.is_closed()) {
    		// look for all possible new edges
    		for (unsigned m = begin; m < end; ++m) {
    			const Cell& left_cell = chart.access(begin,m);
    			if(!left_cell.is_closed()) {
    				const Cell& right_cell = chart.access(m+1,end);
    				if( !right_cell.is_closed())
    					get_candidates(left_cell,right_cell,result_cell);
    			}
    		}
    		// std::cout << result_cell << std::endl;

    		add_unary(result_cell, span == sent_size);

	//	result_cell.apply_beam();
      }
      // std::cout << result_cell << std::endl;
    }
  }
}


inline
void ParserCKYBest::add_unary_init(Cell& cell, bool isroot) const
{

  //for each unary rule set in the grammar [sets made up of all unary rules with a particular rhs]
  for (const auto&  unary_rhs : unary_rhs_from_pos)
  {
    if (cell.exists_edge(0, unary_rhs))
      follow_unary_chain(cell, &cell.at(0,unary_rhs), isroot);
  }
}


inline
void ParserCKYBest::add_unary(Cell& cell, bool isroot) const
{

  //for each unary rule set in the grammar [sets made up of all unary rules with a particular rhs being a lhs of a binary rule]
  for (const auto unary_rhs : unary_rhs_from_binary)
    if (cell.exists_edge(0, unary_rhs))
      follow_unary_chain(cell, &cell.at(0,unary_rhs),isroot);
}

//inline
void ParserCKYBest::follow_unary_chain(Cell& cell, const Edge * edge, bool isroot) const
{
  static int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  std::vector<std::pair<unsigned,const Edge*>> accumulator = {{0,edge}};

  Edge candidate;
  candidate.set_right_child(NULL);

  do {
    const auto& p = accumulator.back();
    unsigned current_height = p.first;
    const Edge * current_edge = p.second;
    accumulator.pop_back();

    if (current_height >= Cell::unary_length - 1)
      continue;
    candidate.set_left_child(current_edge);
    const auto& rules = unary_rhs_2_rules[current_edge->get_lhs()];

    for(const auto& rulep : rules)
    {
      if(isroot || rulep->get_lhs() != start_symbol)
      {
        candidate.set_lhs(rulep->get_lhs());
        candidate.set_probability(current_edge->get_probability() + rulep->get_probability());

        const Edge * new_edge = cell.process_candidate(current_height + 1, candidate);
        if(new_edge && rules_for_unary_exist(new_edge->get_lhs()))
        {
          accumulator.push_back(std::make_pair(current_height + 1, new_edge));
        }
      }
    }
  } while(!accumulator.empty());
}
