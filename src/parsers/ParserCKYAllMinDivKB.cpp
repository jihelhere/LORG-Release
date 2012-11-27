// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMINDIVKB_CPP_
#define _PARSERCKYALLMINDIVKB_CPP_

#include "ParserCKYAllMinDivKB.h"
#include "ParserCKYAll.hpp"

double MinDivProbabilityKB::log_normalisation_factor = 0;
unsigned MinDivProbabilityKB::size = 0;

void
MinDivProbabilityKB::finalize ()
{

}


ostream & MinDivProbabilityKB::operator>> (ostream & out) const 
{return out;}

void
MinDivProbabilityKB::extend_derivation (MinDivProbabilityKB::Edge * edge,
                                        unsigned int i, bool licence_unaries)
{
  if(derivations.size() == i) {
    return;
  }

  if(derivations.size() > 0) {

    packed_edge_probability_with_index& last = derivations[derivations.size() -1];

    //    std::cout << "last.probability " << last.probability << std::endl;

    assert(last.probability <= 0);

    find_succ(edge,last,licence_unaries);
    //    std::cout << "after find_succ" << std::endl;
  }

  if (!candidates.empty()) {

    //get next element from the candidatesidates and append it to derivations
    pop_heap(candidates.begin(),candidates.end());
    derivations.push_back(candidates.back());
    candidates.pop_back();

    //    std::cout << "in edge " << edge << " there are " << derivations.size() << " derivations." << std::endl;

  }



#ifndef NDEBUG
   for(unsigned j = 0; j < derivations.size(); ++j) {
  //   std::cout << "derivations " << j << ": " << derivations[j].probability << " " ;
  //   if(!(derivations[j].dtrs)) { std::cout << "NULL"; }
  //   else {
  //     if(derivations[j].dtrs->is_lexical())
  //       std::cout << *(static_cast<const LexicalPackedEdgeDaughters*>(derivations[j].dtrs)->get_rule());
  //     if(derivations[j].dtrs->is_binary())
  //       std::cout << *(static_cast<const BinaryPackedEdgeDaughters*>(derivations[j].dtrs)->get_rule());
  //     if(!derivations[j].dtrs->is_binary() && !derivations[j].dtrs->is_lexical())
  //       std::cout << *(static_cast<const UnaryPackedEdgeDaughters*>(derivations[j].dtrs)->get_rule());
  //   }
  //   std::cout << std::endl;

    assert(derivations[j].probability <= 0);
  }
#endif

}

void MinDivProbabilityKB::find_succ(PackedEdge<MinDivProbabilityKB>* edge, packed_edge_probability_with_index& pep, bool licence_unaries)
{
  typedef PackedEdge<MinDivProbabilityKB> P;

  if(pep.dtrs->is_lexical())  { return;}
  // binary -> extend left and right daughters
  if(pep.dtrs->is_binary()) {
    const P::BinaryDaughters* d = static_cast<const P::BinaryDaughters*>(pep.dtrs);

    //extend to the left
    unsigned left_pos = d->get_rule()->get_rhs0();
    P* left  = d->left_daughter()->get_edge_ptr(left_pos);
    unsigned nextleft = pep.left_index + 1;
    left->extend_derivation(nextleft+1,true);

    // we haven't reached the expected number of solutions
    if(nextleft < left->get_prob_model().n_deriv()) {

      packed_edge_probability_with_index p(pep);
      p.left_index = nextleft;
      p.probability = Updater::update_maxrule_probability(edge->get_annotations(), *d, log_normalisation_factor, p.left_index, p.right_index);

      assert(p.probability <= 0);

      //      std::cout << p.probability << std::endl;

      // TODO : Find a proper way to remove duplicates !
      if (std::find_if(candidates.begin(), candidates.end(), test_helper(p)) == candidates.end()) {
        candidates.push_back(p);
        push_heap(candidates.begin(), candidates.end());
      }

    }

    //extend to the right
    unsigned right_pos = d->get_rule()->get_rhs1();
    P* right = d->right_daughter()->get_edge_ptr(right_pos);
    unsigned nextright = pep.right_index + 1;

    right->extend_derivation(nextright+1,true);

    if(nextright < right->get_prob_model().n_deriv()) {
      //        std::cout << "bin extending on the right" << std::endl;


      packed_edge_probability_with_index p(pep);
      p.right_index = nextright;
      p.probability = Updater::update_maxrule_probability(edge->get_annotations(), *d, log_normalisation_factor, p.left_index, p.right_index);

      assert(p.probability <= 0);

      //      std::cout << p.probability << std::endl;

      if(std::find_if(candidates.begin(), candidates.end(), test_helper(p)) == candidates.end()){
        candidates.push_back(p);
        push_heap(candidates.begin(), candidates.end());
      }
    }
  }

  //unary
  else {
    if(!licence_unaries) return;

    //      std::cout << "unary case" << std::endl;

    const P::UnaryDaughters* d = static_cast<const P::UnaryDaughters*>(pep.dtrs);

    //        std::cout << * d->get_rule() << std::endl;


    //extend to the left
    unsigned left_pos = d->get_rule()->get_rhs0();
    P* left  = d->left_daughter()->get_edge_ptr(left_pos);
    unsigned nextleft = pep.left_index + 1;

    left->extend_derivation(nextleft+1, false);

    if(nextleft < left->get_prob_model().n_deriv() ) {
      //        std::cout << "un extending" << std::endl;
      packed_edge_probability_with_index p(pep);
      p.left_index = nextleft;
      p.probability = Updater::update_maxrule_probability(edge->get_annotations(), *d, log_normalisation_factor, p.left_index);

      assert(p.probability <= 0);

      //      std::cout << p.probability << std::endl;


      if(std::find_if(candidates.begin(), candidates.end(), test_helper(p)) == candidates.end()){
        candidates.push_back(p);
        push_heap(candidates.begin(), candidates.end());
      }
    }
  }
}



ParserCKYAllMinDivKB::ParserCKYAllMinDivKB(std::vector<AGrammar*>& cgs,
                                             const std::vector<double>& p, double b_t,
                                             const annot_descendants_type& annot_descendants_,
                                             bool accurate_, unsigned min_beam, int stubborn, unsigned k_)
: ParserCKYAll_Impl<MinDivProbabilityCell>(cgs, p, b_t, annot_descendants_, accurate_, min_beam, stubborn) , k(k_)
{
  // this is not in the super class because maxn parsing uses a
  //different mapping
  //create the coarse-to-fine map
  this->create_coarse_to_fine_mapping(this->grammars);

  Edge::set_unary_chains(this->grammars[this->grammars.size() - 1]->get_unary_decoding_paths());
}


void ParserCKYAllMinDivKB::extend_all_derivations()
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


void ParserCKYAllMinDivKB::extract_solution()
{
  //  std::cout << "in extract" << std::endl;


  compute_inside_outside_probabilities();

  initialise_candidates();

  //  std::cout << "after init cand" << std::endl;

  extend_all_derivations();
}


void ParserCKYAllMinDivKB::initialise_candidates()
{

  double sentence_probability = std::log(get_sentence_probability());
  //  unsigned sent_size = chart->get_size();

  MinDivProbabilityKB::set_log_normalisation_factor(sentence_probability);
  MinDivProbabilityKB::set_size(k);

//   copied from ParserCKYAllMaxRuleKB
//   calculate_maxrule_probabilities();
}





#endif /* _PARSERCKYALLMINDIVKB_CPP_ */
