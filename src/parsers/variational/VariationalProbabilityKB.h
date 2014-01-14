// -*- mode: c++ -*-
#ifndef _VARIATIONALPROBABILITYKB_H_
#define _VARIATIONALPROBABILITYKB_H_

#include "edges/PackedEdgeProbability.h"
#include "edges/PackedEdge.h"
#include "parsers/maxrule/MaxRuleTreeLogProbaComputer.h"
#include "emptystruct.h"
#include "ChartCKY.h"

#include "edges/PackedEdgeDaughters.h"
#include "VariationalTypes.h"

#include "parsers/ParserCKYAllFactory.h"

class VariationalProbabilityKB
{

public:

  typedef std::vector<packed_edge_probability> heap_type;
  typedef typename VariationalTypes::Edge Edge;
  typedef typename VariationalTypes::Cell Cell;
  typedef MaxRuleTreeLogProbaComputer<VariationalProbabilityKB> QInsideComputer;

  typedef typename VariationalTypes::LexicalDaughter LexicalDaughter;
  typedef typename VariationalTypes::UnaryDaughter UnaryDaughter;
  typedef typename VariationalTypes::BinaryDaughter BinaryDaughter;

  typedef typename VariationalTypes::LRule LRule;

  static void set_calculation(ParserCKYAllFactory::MaxParsing_Calculation c) {QInsideComputer::set_calculation(c);}


private:

  //TODO refactor
  //k-best attributes
  heap_type candidates;
  heap_type derivations;
  static unsigned size;

  //Variational normalisation factor
  // actually log-marginals
  double marginals;


public:

  VariationalProbabilityKB() :  candidates(), derivations(), marginals() {candidates.reserve(50);};
  ~VariationalProbabilityKB() {};

  inline static void set_size(unsigned k) {size = k;}

  inline const heap_type & get_candidates() const { return candidates; }
  inline const heap_type & get_derivations() const { return derivations; }

  inline const packed_edge_probability& get(unsigned idx) const {return derivations[idx];}
  inline packed_edge_probability& get(unsigned idx) { return derivations[idx]; }

  inline void init() {}
  inline void store_marginals(const Edge& edge);
  inline void update_lexical(Edge& e, LexicalDaughter& dtr, double /*unused*/);
  inline void update_unary(Edge& e, UnaryDaughter& dtr, double /*unused*/);
  inline void update_binary(Edge& e, BinaryDaughter& dtr, double /*unused*/);
  inline void finalize();

  inline void find_succ(Edge*,packed_edge_probability& pep, bool licence_unaries);
  inline void extend_derivation(Edge*, unsigned, bool, const std::vector<double>& /*unused*/) ;

  inline unsigned n_deriv() const {return derivations.size();};

  inline bool has_solution(unsigned i) const {return i <derivations.size();}

private:

  struct test_helper
  {
    const packed_edge_probability& pep;
    test_helper(const packed_edge_probability& p) : pep(p) {};

    inline bool operator()(const packed_edge_probability& p)
    {
      return (p.probability == pep.probability) //|| (p.dtrs == pep.dtrs)
      ;
    }
  };

  public:
    inline std::ostream& operator>>(std::ostream& out) const;
};


inline std::ostream& operator<<(std::ostream& out, const VariationalProbabilityKB & prob)
{
  return out << "((VariationalKBProb: " << &prob
             << "): nb_deriv." << prob.get_derivations().size()
             << " nb_candid." << prob.get_candidates().size() << ")";
}


void VariationalProbabilityKB::update_lexical(Edge& e, LexicalDaughter& dtr, double /*unused*/)
{
   //BLOCKTIMING("VariationalProbabilityKB::update_lexical");
  const AnnotationInfo & a = e.get_annotations();
  const LRule* rule = dtr.get_rule();
  assert(rule != NULL);

  packed_edge_probability pep;
  pep.probability = QInsideComputer::compute(a, rule, marginals);

//    std::cout << "lexical " << pep.probability << std::endl;
  assert(pep.probability <=0);

  pep.dtrs = &dtr;

  candidates.push_back(pep);

  if(derivations.empty())
     derivations.push_back(pep);
  else if(pep.probability > derivations[0].probability)
      derivations[0] = pep;

//   std::cout << *this << std::endl;
}

void VariationalProbabilityKB::update_unary(Edge& e, UnaryDaughter& dtr, double /*unused*/)
{
  //BLOCKTIMING("VariationalProbabilityKB::update_unary");
  const AnnotationInfo & a = e.get_annotations();
  packed_edge_probability pep;
  pep.dtrs = &dtr;
  //  std::cout << "before ump" << std::endl;
  pep.probability= QInsideComputer::compute(a, dtr, marginals);

//    std::cout << "unary "<< pep.probability << std::endl;
  assert(pep.probability <=0);

  candidates.push_back(pep);

  if(derivations.empty())
       derivations.push_back(pep);
  else if(pep.probability > derivations[0].probability)
       derivations[0] = pep;

//   std::cout << *this << std::endl;
}

void VariationalProbabilityKB::update_binary(Edge& e, BinaryDaughter& dtr, double /*unused*/)
{
  //BLOCKTIMING("VariationalProbabilityKB::update_binary");
  const AnnotationInfo & a = e.get_annotations();
  packed_edge_probability pep;
  pep.dtrs = &dtr;


  pep.probability= QInsideComputer::compute(a, dtr, marginals);

  //  std::cout << candidates.size() << std::endl;
//   std::cout << "binary " << pep.probability << std::endl;
  //  std::cout << pep.dtrs << std::endl;


  //  std::cout << pep.probability << std::endl;
  assert(pep.probability <=0);

  candidates.push_back(pep);


  if(derivations.empty())
    derivations.push_back(pep);
  else if(pep.probability > derivations[0].probability)
    derivations[0] = pep;

//     std::cout << *this << std::endl;
}

struct gt_pep
{
  bool operator()(const packed_edge_probability& p1, const packed_edge_probability& p2) const
  {
    return p1 > p2;
  }
};



void VariationalProbabilityKB:: finalize()
{
  //  std::cout << "size candidates: " << candidates.size() << std::endl;

  if(!candidates.empty()) {
    if(candidates.size() > size) {

      // std::cout << "BEFORE NTH " << size << std::endl;
      // for (unsigned i = 0; i < candidates.size(); ++i)
      //   {
      //     std::cout << candidates[i].probability << " ";

      //     if(candidates[i].dtrs->is_lexical())
      //       std::cout << *(static_cast<const LexicalPackedEdgeDaughters*>(candidates[i].dtrs)->get_rule());
      //     if(candidates[i].dtrs->is_binary())
      //       std::cout << *(static_cast<const BinaryPackedEdgeDaughters*>(candidates[i].dtrs)->get_rule());
      //     if(!candidates[i].dtrs->is_binary() && !candidates[i].dtrs->is_lexical())
      //       std::cout << *(static_cast<const UnaryPackedEdgeDaughters*>(candidates[i].dtrs)->get_rule());
      //     std::cout << std::endl;
      //   }

      std::nth_element(candidates.begin(),candidates.begin()+size,candidates.end(), gt_pep());
      candidates.resize(size);

      //       std::cout << "AFTER NTH " << size << std::endl;
      // for (int i = 0; i < candidates.size(); ++i)
      //   {
      //     std::cout << candidates[i].probability << " ";

      //     if(candidates[i].dtrs->is_lexical())
      //       std::cout << *(static_cast<const LexicalPackedEdgeDaughters*>(candidates[i].dtrs)->get_rule());
      //     if(candidates[i].dtrs->is_binary())
      //       std::cout << *(static_cast<const BinaryPackedEdgeDaughters*>(candidates[i].dtrs)->get_rule());
      //     if(!candidates[i].dtrs->is_binary() && !candidates[i].dtrs->is_lexical())
      //       std::cout << *(static_cast<const UnaryPackedEdgeDaughters*>(candidates[i].dtrs)->get_rule());
      //     std::cout << std::endl;
      //   }



    }
    std::make_heap(candidates.begin(),candidates.end());

    std::pop_heap(candidates.begin(),candidates.end());
    candidates.pop_back();
  }
}

void VariationalProbabilityKB::extend_derivation(Edge* edge, unsigned i, bool licence_unaries, const std::vector<double>& /*  unused */)
{
  if(derivations.size() == i) {
    return;
  }

  if(derivations.size() > 0) {

    packed_edge_probability& last = derivations[derivations.size() -1];

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

void VariationalProbabilityKB::find_succ(Edge* edge, packed_edge_probability& pep, bool licence_unaries)
{

  if(pep.dtrs->is_lexical())  { return;}
  // binary -> extend left and right daughters
  if(pep.dtrs->is_binary()) {
    const BinaryDaughter* d = static_cast<const BinaryDaughter*>(pep.dtrs);

    //extend to the left
    Edge& left  = d->left_daughter();
    unsigned nextleft = pep.left_index + 1;
    left.extend_derivation(nextleft+1,true, std::vector<double>());

    // we haven't reached the expected number of solutions
    if(nextleft < left.get_prob_model().n_deriv()) {

      packed_edge_probability p(pep);
      p.left_index = nextleft;
      p.probability = QInsideComputer::compute(edge->get_annotations(), *d, marginals, p.left_index, p.right_index);

      assert(p.probability <= 0);

      //      std::cout << p.probability << std::endl;

      // TODO : Find a proper way to remove duplicates !
      if (std::find_if(candidates.begin(), candidates.end(), test_helper(p)) == candidates.end()) {
        candidates.push_back(p);
        push_heap(candidates.begin(), candidates.end());
      }

    }

    //extend to the right
    Edge& right = d->right_daughter();
    unsigned nextright = pep.right_index + 1;

    right.extend_derivation(nextright+1,true, std::vector<double>());

    if(nextright < right.get_prob_model().n_deriv()) {
      //        std::cout << "bin extending on the right" << std::endl;


      packed_edge_probability p(pep);
      p.right_index = nextright;
      p.probability = QInsideComputer::compute(edge->get_annotations(), *d, marginals, p.left_index, p.right_index);

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

    const UnaryDaughter* d = static_cast<const UnaryDaughter*>(pep.dtrs);

    //        std::cout << * d->get_rule() << std::endl;


    //extend to the left
    Edge& left  = d->left_daughter();
    unsigned nextleft = pep.left_index + 1;

    left.extend_derivation(nextleft+1, false, std::vector<double>());

    if(nextleft < left.get_prob_model().n_deriv() ) {
      //        std::cout << "un extending" << std::endl;
      packed_edge_probability p(pep);
      p.left_index = nextleft;
      p.probability = QInsideComputer::compute(edge->get_annotations(), *d, marginals, p.left_index);

      assert(p.probability <= 0);

      //      std::cout << p.probability << std::endl;


      if(std::find_if(candidates.begin(), candidates.end(), test_helper(p)) == candidates.end()){
        candidates.push_back(p);
        push_heap(candidates.begin(), candidates.end());
      }
    }
  }
}


inline void VariationalProbabilityKB::store_marginals(const Edge& edge)
{
  marginals = edge.marginalise();
}

#endif /* _VARIATIONALPROBABILITYKB_H_ */
