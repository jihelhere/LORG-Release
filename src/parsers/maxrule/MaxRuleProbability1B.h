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
public:
  MaxRuleProbability1B() : best() {};
  ~MaxRuleProbability1B() {};


  static void set_calculation(ParserCKYAllFactory::MaxParsing_Calculation c) {QInsideComputer::set_calculation(c);}

  inline void init() {best.init();}
  inline const packed_edge_probability& get(unsigned/*ignored*/) const {return best;}
  inline packed_edge_probability& get(unsigned /*ignored*/) {return best;}

  inline void update_lexical(Edge& e, LexicalDaughter& dtr, double log_normalisation_factor);
  inline void update_unary(Edge& e,   UnaryDaughter& dtr, double log_normalisation_factor);
  inline void update_binary(Edge& e,  BinaryDaughter& dtr, double log_normalisation_factor);
  inline void finalize();

  inline unsigned n_deriv() const {return 1;}

  inline bool has_solution(unsigned i) const {return i == 0;} ;
};

inline std::ostream& operator<<(std::ostream& out, const MaxRuleProbability1B & prob)
{
  return out << "((MaxRule1BProb: " << &prob << ")";
}


inline void MaxRuleProbability1B::update_lexical(Edge & edge, LexicalDaughter& dtr, double log_normalisation_factor)
{
    if (not dtr.is_calculated())
    {
      const LexicalRuleC2f* rule = dtr.get_rule();
      dtr.set_q_score(QInsideComputer::compute(edge.get_annotations(), rule, log_normalisation_factor));
      dtr.set_calculated(true);
    }

    double probability = dtr.get_q_score() + dtr.get_relaxation();

    // if (dtr.get_relaxation() != 0.0)
    //   std::cout << "lexical relax "<< dtr.get_relaxation() << std::endl;

    if (probability > best.probability)
    {
        best.probability = probability;
        best.dtrs = &dtr;
    }
}


inline void MaxRuleProbability1B::update_unary (Edge & e, UnaryDaughter & dtr, double log_normalisation_factor)
{
  //  std::cout << "update with " << *(dtr.get_rule()) << std::endl;

  const Edge& left  = dtr.left_daughter();

  double probability = - std::numeric_limits<double>::infinity();

  if(left.get_prob_model().get(0).dtrs &&
     (left.get_prob_model().get(0).dtrs->is_lexical() || left.get_prob_model().get(0).dtrs->is_binary()))
  {
    if (not dtr.is_calculated())
    {
      const AnnotationInfo& left_annot = dtr.left_daughter().get_annotations();
      const auto& probs = dtr.get_rule()->get_probability();

      dtr.set_q_score(QInsideComputer::compute_simple(e.get_annotations(), log_normalisation_factor, left_annot, probs));
      dtr.set_calculated(true);
    }

    probability = dtr.get_q_score() + dtr.get_relaxation()
                  + dtr.left_daughter().get_prob_model().get(0).probability;

    // if (dtr.get_relaxation() != 0.0)
    //   std::cout << "unary relax "<< dtr.get_relaxation() << std::endl;
  }

  if (probability > best.probability)
  {
    best.probability = probability;
    best.dtrs = &dtr;
  }

  //  std::cout << probability << std::endl;
}

inline void
MaxRuleProbability1B::update_binary (Edge & e, BinaryDaughter & dtr, double log_normalisation_factor)
{
  //std::cout << "update with " << *(dtr.get_rule()) << std::endl;

    if (not dtr.is_calculated())
    {
      const AnnotationInfo& left_annot = dtr.left_daughter().get_annotations();
      const AnnotationInfo& right_annot = dtr.right_daughter().get_annotations();
      const auto& probs = dtr.get_rule()->get_probability();

      dtr.set_q_score(QInsideComputer::compute_simple(e.get_annotations(), log_normalisation_factor,
                                                      left_annot, right_annot, probs));
      dtr.set_calculated(true);
    }

    double probability = dtr.get_q_score();

    // if score == -inf then dtr might not be initialized correctly
    // TODO: FIX THIS
    if(probability != -std::numeric_limits<double>::infinity())
    {
      probability += dtr.get_relaxation()
                     + dtr.left_daughter().get_prob_model().get(0).probability
                     + dtr.right_daughter().get_prob_model().get(0).probability;

    // if (dtr.get_relaxation() != 0.0)
    //   std::cout << "binary relax "<< dtr.get_relaxation() << std::endl;
    }
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
