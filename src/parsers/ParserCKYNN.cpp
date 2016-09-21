#include "ParserCKYNN.h"
#include "SimpleChartCKY.hpp"

#define THRESHOLD -9

void ParserCKYNN::parse(Chart& chart, scorer& s) const
{
  bool isroot = chart.get_size() == 1;

  for(unsigned i = 0; i < chart.get_size(); ++i)
  {
    //std::cerr << "add unary init " << i << endl;
    add_unary_init(chart.access(i,i),isroot, s);
  }

  process_internal_rules(chart, s);
}



inline
void ParserCKYNN::add_unary_init(Cell& cell, bool isroot, scorer& s) const
{

  //for each unary rule set in the grammar [sets made up of all unary rules with a particular rhs]
  for (const auto& unary_rhs : unary_rhs_from_pos)
  {
     // std::cerr << "preterm is " << unary_rhs
     //          << " " << SymbolTable::instance_nt().translate(unary_rhs)
     //          << std::endl;

    if (cell.exists_edge(0, unary_rhs))
    {
      follow_unary_chain(cell,&cell.at(0,unary_rhs),isroot, s);
    }
  }
}


inline
void ParserCKYNN::add_unary(Cell& cell, bool isroot, scorer& s) const
{

  //for each unary rule set in the grammar [sets made up of all unary rules with a particular rhs being a lhs of a binary rule]
  for (const auto unary_rhs : unary_rhs_from_binary)
  {

     // std::cerr << "preterm is " << unary_rhs
     //          << " " << SymbolTable::instance_nt().translate(unary_rhs)
     //          << std::endl;

    if (cell.exists_edge(0, unary_rhs))
      follow_unary_chain(cell,&cell.at(0,unary_rhs),isroot, s);
  }
}





//inline
void ParserCKYNN::follow_unary_chain(Cell& cell, const Edge * edge, bool isroot, scorer& s) const
    {

      // std::cerr << "follow_unary_chain" << std::endl;

  static int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  std::vector<std::pair<unsigned,const Edge*>> accumulator = {{0,edge}};

  Edge candidate;
  candidate.set_right_child(NULL);

  size_t size = cell.get_top() ? Cell::unary_length +1 : Cell::unary_length;

  do {
    const auto& p = accumulator.back();
    unsigned current_height = p.first;
    const Edge * current_edge = p.second;
    accumulator.pop_back();



    if (current_height < size - 1)
    {

      //std::cerr << current_height << std::endl;

    candidate.set_left_child(current_edge);

    if (current_edge->get_lhs() >= (int) unary_rhs_2_rules.size()) continue;
    const auto& rules = unary_rhs_2_rules[current_edge->get_lhs()];

    // std::cerr << SymbolTable::instance_nt().translate(current_edge->get_lhs())
    //           << " " << current_edge->get_lhs() << "/" << unary_rhs_2_rules.size()
    //           << " L: " << rules.size() << std::endl;



    //if(rules.size() == 0) continue; //WTF!!
    for(const auto& rulep : rules)
    {
      //std::cerr << * static_cast<const Production*>(rulep) << std::endl;
      if(isroot || rulep->get_lhs() != start_symbol)
      {
        double pruning_probability = current_edge->get_pruning_probability() +
                                     static_cast<const Rule*>(rulep)->get_probability();
        if (pruning_probability < THRESHOLD)
          continue;

        candidate.set_lhs(rulep->get_lhs());
        candidate.set_pruning_probability(pruning_probability);
        double prob = current_edge->get_probability() +
                      s.compute_unary_score(cell.get_begin(), cell.get_end(), rulep);
        candidate.set_probability(prob);

        const Edge * new_edge = cell.process_candidate(current_height + 1, candidate);

        if (new_edge && new_edge->get_probability() == prob)
        {
          s.register_last_expression(new_edge);
        }


        if(new_edge && rules_for_unary_exist(new_edge->get_lhs()))
        {
          accumulator.push_back(std::make_pair(current_height + 1, new_edge));
        }
      }
    }
    }
  } while(!accumulator.empty());

  //std::cerr << "acc is empty" << std::endl;

}







void ParserCKYNN::process_internal_rules(Chart& chart, scorer& s) const
{
  unsigned sent_size=chart.get_size();
  for (unsigned span = 2; span <= sent_size; ++span)
  {
    //std::cerr << "span: " << span << std::endl;
    unsigned end_of_begin=sent_size-span;
    for (unsigned begin=0; begin <= end_of_begin; ++begin)
    {
      unsigned end = begin + span -1;
      //      std::cerr << "begin: " << begin << ", end: " << end << std::endl;

      Cell& result_cell = chart.access(begin,end);

      if(!result_cell.is_closed())
      {
        // look for all possible new edges
        for (unsigned m = begin; m < end; ++m)
        {
          const Cell& left_cell = chart.access(begin,m);
          if(!left_cell.is_closed())
          {
            const Cell& right_cell = chart.access(m+1,end);
            if( !right_cell.is_closed())
              get_candidates(left_cell,right_cell,result_cell, s);
          }
        }
        // std::cout << result_cell << std::endl;

        add_unary(result_cell, span == sent_size, s);

	//	result_cell.apply_beam();
      }
      // std::cout << result_cell << std::endl;
    }
  }
}


inline
void ParserCKYNN::get_candidates(const Cell& left_cell,
				   const Cell& right_cell,
				   Cell& result_cell,
                                   scorer& s) const
{
  Edge current_candidate;


  //iterating through all the rules P -> L R, indexed by L
  for(const auto& same_rhs0 : brules)
  {
    // is L present in  left_cell ?
    if(left_cell.exists_edge_all_height(same_rhs0.rhs0)) {
      const Edge& left_edge = left_cell.get_best_edge(same_rhs0.rhs0);
      current_candidate.set_left_child(&left_edge);

      //iterating through all the rules P -> L R, indexed by R, L fixed
      for (const auto& same_rhs1 : same_rhs0)
      {

	// is R present in right_cell ?
	if(right_cell.exists_edge_all_height(same_rhs1.rhs1)) {

	  const Edge& right_edge = right_cell.get_best_edge(same_rhs1.rhs1);
	  current_candidate.set_right_child(&right_edge);

	  double pruprob1 = left_edge.get_pruning_probability() + right_edge.get_pruning_probability();
          if (pruprob1 < THRESHOLD) continue;

	  double prob1 = left_edge.get_probability() + right_edge.get_probability();

	  //iterating through all the rules P -> L R, indexed by P, R and L fixed

	  for(const auto& b : same_rhs1)
          {
            double pruprob = pruprob1 + static_cast<const Rule*>(b)->get_probability();
            if (pruprob < THRESHOLD)
              continue;

            current_candidate.set_pruning_probability(pruprob);
	    current_candidate.set_lhs(b->get_lhs());
            double prob = prob1 + s.compute_binary_score(result_cell.get_begin(),
                                                         result_cell.get_end(),
                                                         right_cell.get_begin(),
                                                         b);
	    current_candidate.set_probability(prob);

	    //	    std::cout << *(*bitr) << std::endl;

	    auto ep = result_cell.process_candidate(0, current_candidate);
            if (ep and ep->get_probability() == prob)
              s.register_last_expression(ep);
	  }
	}
      }
    }
  }
}
