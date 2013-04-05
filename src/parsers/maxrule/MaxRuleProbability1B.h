// // -*- mode: c++ -*-
#ifndef _MAXRULEPROBABILITY1B_H_
#define _MAXRULEPROBABILITY1B_H_
#include <numeric>

#include "edges/PackedEdgeProbability.h"
#include "edges/PackedEdge.h"
#include "MaxRuleTreeLogProbaComputer.h"
#include "emptystruct.h"
#include "ChartCKY.h"

#include "parsers/ParserCKYAllFactory.h"


class MaxRuleProbability1B;

struct MaxRule1BTypes {
  typedef MaxRuleProbability1B EdgeProbability ;
  typedef emptystruct EdgeDaughterProbability ;
  typedef Word ChartWord ;

  typedef BRuleC2f BRule;
  typedef URuleC2f URule;
  typedef LexicalRuleC2f LRule;
  typedef PackedEdge< MaxRule1BTypes > Edge ;
  typedef PCKYAllCell< MaxRule1BTypes > Cell ;
  typedef ChartCKY< MaxRule1BTypes > Chart ;
  typedef BinaryPackedEdgeDaughters<MaxRule1BTypes> BinaryDaughter;
  typedef UnaryPackedEdgeDaughters<MaxRule1BTypes>  UnaryDaughter;
  typedef LexicalPackedEdgeDaughters<MaxRule1BTypes> LexicalDaughter;
};


class MaxRuleProbability1B
{
public:
  typedef typename MaxRule1BTypes::Edge Edge;
  typedef typename MaxRule1BTypes::Cell Cell;
  typedef typename MaxRule1BTypes::UnaryDaughter UnaryDaughter;
  typedef typename MaxRule1BTypes::BinaryDaughter BinaryDaughter;
  typedef typename MaxRule1BTypes::LexicalDaughter LexicalDaughter;
  typedef MaxRuleTreeLogProbaComputer<MaxRuleProbability1B> QInsideComputer;

private:
  packed_edge_probability best;
  static double log_normalisation_factor;
public:
  MaxRuleProbability1B() : best() {};
  ~MaxRuleProbability1B() {};

  static void set_log_normalisation_factor(double lnf);
  static void set_calculation(ParserCKYAllFactory::MaxParsing_Calculation c) {QInsideComputer::set_calculation(c);}

  inline void init() {best.init();}
  inline const packed_edge_probability& get(unsigned/*ignored*/) const {return best;}
  inline packed_edge_probability& get(unsigned /*ignored*/) {return best;}

  inline void update_lexical(Edge& e,  const LexicalDaughter& dtr);
  inline void update_unary(Edge& e, const UnaryDaughter& dtr);
  inline void update_binary(Edge& e, const BinaryDaughter& dtr);
  inline void finalize();

  inline unsigned n_deriv() const {return 1;}

  inline bool has_solution(unsigned i) const {return i == 0;} ;
};

inline std::ostream& operator<<(std::ostream& out, const MaxRuleProbability1B & prob)
{
  return out << "((MaxRule1BProb: " << &prob << ")";
}





inline void MaxRuleProbability1B::update_lexical(Edge & edge,  const LexicalDaughter& dtr)
{
    const LexicalRuleC2f* rule = dtr.get_rule();

    double probability = QInsideComputer::compute(edge.get_annotations(), rule, log_normalisation_factor) + dtr.get_relaxation();

    // if (dtr.get_relaxation() != 0.0)
    //   std::cout << "lexical relax "<< dtr.get_relaxation() << std::endl;


    if (probability > best.probability) {
        best.probability = probability;
        best.dtrs = &dtr;
    }
}


inline void MaxRuleProbability1B::update_unary (Edge & e, const UnaryDaughter & dtr)
{
  //  std::cout << "update with " << *(dtr.get_rule()) << std::endl;

  double probability = -std::numeric_limits<double>::infinity();

  const Edge& left  = dtr.left_daughter();
  // std::cout << left << std::endl;
  // std::cout << left.get_prob_model().get(0).dtrs << std::endl;
  if(left.get_prob_model().get(0).dtrs && (left.get_prob_model().get(0).dtrs->is_lexical() || left.get_prob_model().get(0).dtrs->is_binary())) {
    //    std::cout << "unary case" << std::endl;

    // if (dtr.get_relaxation() != 0.0)
    //   std::cout << "unary relax "<< dtr.get_relaxation() << std::endl;


    probability =  QInsideComputer::compute(e.get_annotations(), dtr, log_normalisation_factor) + dtr.get_relaxation();
  }

  if (probability > best.probability) {
    best.probability = probability;
    best.dtrs = &dtr;
  }

  //  std::cout << probability << std::endl;
}

inline void
MaxRuleProbability1B::update_binary (Edge & e, const BinaryDaughter & dtr)
{
  //  std::cout << "update with " << *(dtr.get_rule()) << std::endl;

  double probability = QInsideComputer::compute(e.get_annotations(), dtr, log_normalisation_factor) + dtr.get_relaxation();

  // if (dtr.get_relaxation() != 0.0)
  //   std::cout << "binary relax "<< dtr.get_relaxation() << std::endl;

  if (probability > best.probability)
    {
      best.probability = probability;
      best.dtrs = &dtr;
    }

  //  std::cout << probability << std::endl;
}

// //uncomment the function to get the best indexes
// // only useful if you want to print them out
inline void MaxRuleProbability1B::finalize() {}
// // {
// //   if(up) {
// //     //   std::cout << "here" << std::endl;
// //     if(best.dtrs !=  NULL) {
// //       if (best.dtrs->is_binary())
// // 	{
// // 	  //	  std::cout << "there" << std::endl;
// // 	  compute_best_indexes(up->get_annotations(),
// // 			       *static_cast<const BinaryPackedEdgeDaughters*>(best.dtrs),
// // 			       log_normalisation_factor,best.left_index,best.right_index);
// // 	}
// //       else { //unary
// // 	compute_best_indexes(up->get_annotations(),
// // 			     *static_cast<const UnaryPackedEdgeDaughters*>(best.dtrs),
// // 			     log_normalisation_factor,best.left_index);
// //       }
// //     }
// //   }
// // }







////////////



#endif /* _MAXRULEPROBABILITY_H_ */