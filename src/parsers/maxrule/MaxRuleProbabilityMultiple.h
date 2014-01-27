// -*- mode: c++ -*-

#ifndef _MAXRULEMULTIPLEPROBABILITY_H_
#define _MAXRULEMULTIPLEPROBABILITY_H_

#include "edges/PackedEdgeProbability.h"
#include "edges/PackedEdge.h"
#include "MaxRuleTreeLogProbaComputer.h"
#include "emptystruct.h"
#include "ChartCKY.h"

#include "parsers/ParserCKYAllFactory.h"

#include <vector>
#include <unordered_map>

class MaxRuleProbabilityMultiple;

struct MaxRuleMultipleTypes {
  typedef MaxRuleProbabilityMultiple EdgeProbability ;
  typedef emptystruct EdgeDaughterProbability ;
  typedef Word ChartWord ;

  typedef BRuleC2f BRule;
  typedef URuleC2f URule;
  typedef LexicalRuleC2f LRule;
  typedef PackedEdge< MaxRuleMultipleTypes > Edge ;
  typedef PCKYAllCell< MaxRuleMultipleTypes > Cell ;
  typedef ChartCKY< MaxRuleMultipleTypes > Chart ;
  typedef BinaryPackedEdgeDaughters<MaxRuleMultipleTypes> BinaryDaughter;
  typedef UnaryPackedEdgeDaughters<MaxRuleMultipleTypes>  UnaryDaughter;
  typedef LexicalPackedEdgeDaughters<MaxRuleMultipleTypes> LexicalDaughter;
};



class MaxRuleProbabilityMultiple
{
private:
  typedef std::vector<packed_edge_probability> heap_type;
  typedef MaxRuleTreeLogProbaComputer<MaxRuleProbabilityMultiple> QInsideComputer;

  std::vector<AnnotationInfo> annotations_backup;

  heap_type candidates;
  heap_type derivations;

  static unsigned size;
  static unsigned nb_grammars;


public:
  typedef typename MaxRuleMultipleTypes::Edge Edge;
  typedef typename MaxRuleMultipleTypes::Cell Cell;
  typedef typename MaxRuleMultipleTypes::UnaryDaughter UnaryDaughter;
  typedef typename MaxRuleMultipleTypes::BinaryDaughter BinaryDaughter;
  typedef typename MaxRuleMultipleTypes::LexicalDaughter LexicalDaughter;

  MaxRuleProbabilityMultiple() : candidates(),
                                 derivations(heap_type(1))
  {candidates.reserve(50);};
  ~MaxRuleProbabilityMultiple() {};

  inline static void set_size(unsigned k)       {size = k;}

  inline static void set_nbgrammars(unsigned n) {nb_grammars = n;}

  static void set_calculation(ParserCKYAllFactory::MaxParsing_Calculation c) {QInsideComputer::set_calculation(c);}


  inline void init()
  {
    //candidates = heap_type(0);
    derivations =  heap_type(1);
  }

  inline const packed_edge_probability& get(unsigned idx) const {return derivations[idx];}
  inline       packed_edge_probability& get(unsigned idx)       {return derivations[idx];}

  inline void update_lexical(Edge&, const LexicalDaughter&, double /*unused*/)
  { throw std::runtime_error("I shall not be called");}
  inline void update_unary(Edge&, const UnaryDaughter&, double /*unused*/)
  {throw std::runtime_error("I shall not be called");};
  inline void update_binary(Edge&, const BinaryDaughter&, double /*unused*/)
  {throw std::runtime_error("I shall not be called");};
  inline void finalize();

  inline void pick_best_lexical(const Edge& e, LexicalDaughter& dtr, const std::vector<double>& log_norms);
  inline void pick_best_binary(const Edge& e, BinaryDaughter& dtr, const std::vector<double>& log_norms);
  inline void pick_best_unary(const Edge& e, UnaryDaughter& dtr, const std::vector<double>& log_norms);
  inline void pick_best();

  inline void find_succ(Edge*,packed_edge_probability& pep, bool licence_unaries, const std::vector<double>& log_norms);
  inline void extend_derivation(Edge*, unsigned, bool, const std::vector<double>& log_norms);


  inline unsigned n_deriv() const {return derivations.size();}
  inline bool has_solution(unsigned i) const
  {
    //    std::cout << "i " << i << std::endl;
    //    std::cout << "derivations.size() " << derivations.size() << std::endl;

    // if(i < derivations.size())
    //   std::cout << derivations[i].probability << std::endl;;

    return
      i < derivations.size() ;
    //&& derivations[i].probability != -std::numeric_limits<double>::infinity();
  }

  inline       std::vector<AnnotationInfo>& get_annotations_backup();
  inline const std::vector<AnnotationInfo>& get_annotations_backup() const;

  inline void backup_annotations(const AnnotationInfo& annotations);

private:

  // struct test_helper
  // {
  //   const packed_edge_probability& pep;
  //   test_helper(const packed_edge_probability& p) : pep(p) {};

  //   bool operator()(const packed_edge_probability& p)
  //   {
  //     //      return false;
  //     return p.probability == pep.probability
  //       //  || p.dtrs == pep.dtrs
  //       ;}
  // };

};
inline std::ostream& operator<<(std::ostream& out, const MaxRuleProbabilityMultiple & prob)
{
  return out << "(MaxRuleMulltipleProb: " << &prob << ")";
}



void MaxRuleProbabilityMultiple::finalize()
{
}


void MaxRuleProbabilityMultiple::pick_best_lexical(const Edge& e,
                                                   LexicalDaughter & dtr, const std::vector<double>& lnfs)
{
  packed_edge_probability p;
  p.dtrs = &dtr;

  LexicalDaughter* d = static_cast<LexicalDaughter*>(p.dtrs);

  if (p.dtrs->is_calculated())
  {
    p.probability = p.dtrs->get_q_score();
  }
  else
  {
    const std::vector<AnnotationInfo>& upannots = get_annotations_backup();

    p.probability = 0;
    for (unsigned i = 0; i < upannots.size(); ++i)
    {
      if(lnfs[i] == -std::numeric_limits<double>::infinity())
        continue;

      const std::vector<double>& rule_probs = d->get_rule()->get_coarser(upannots.size() - i - 1)->get_probability();

      double probability = 0;
      for(unsigned j = 0; j < rule_probs.size(); ++j)
      {
        if(!upannots[i].valid_prob_at(j, LorgConstants::NullProba)) continue;
        //std::cout << "out: " << upannots[i].outside_probabilities.array[j] << std::endl;
        //std::cout << "prb: " << rule_probs[j] << std::endl;
        probability += rule_probs[j] * upannots[i].outside_probabilities.array[j];
      }


      double logprob =
          // (get_log_normalisation_factor(i) == -std::numeric_limits<double>::infinity()) ?
          // -std::numeric_limits<double>::infinity() :
          std::log(probability) - lnfs[i];

      // ???
      //    p.probability += logprob + std::exp(logprob);
      p.probability += logprob;

      if(p.probability ==-std::numeric_limits<double>::infinity())
        //std::cout << "it's happening!" << std::endl;
        break;
    }

    p.dtrs->set_calculated(true);
    p.dtrs->set_q_score(p.probability);
  }

  p.probability += e.get_relaxation();

  //p.probability += std::exp(p.probability);

  if (candidates.empty() || p.probability > derivations[0].probability)
    derivations[0] = p;
  candidates.emplace_back(p);
}

void MaxRuleProbabilityMultiple::pick_best_binary(const Edge& e,
                                                  BinaryDaughter& dtr, const std::vector<double>& log_norms)
{
  //    std::cout << "binary case" << std::endl;

  packed_edge_probability p;
  p.dtrs = &dtr;

  BinaryDaughter * d = static_cast<BinaryDaughter*>(p.dtrs);
  Edge& left  = d->left_daughter();
  Edge& right = d->right_daughter();

  if (p.dtrs->is_calculated())
  {
    p.probability = p.dtrs->get_q_score();
  }
  else
  {
    const std::vector<AnnotationInfo>& upannots = get_annotations_backup();
    const std::vector<AnnotationInfo>& leftannots = left.get_prob_model().get_annotations_backup();
    const std::vector<AnnotationInfo>& rightannots = right.get_prob_model().get_annotations_backup();

    p.probability = 0;

    for (unsigned i = 0; i < upannots.size(); ++i)
    {
      if(log_norms[i] == -std::numeric_limits<double>::infinity())
        continue;
      //          std::cout << *( d->get_rule()->get_coarser(upannots.size() - i + 1)) << std::endl;

      const std::vector<std::vector<std::vector<double> > >& rule_probs =
          d->get_rule()->get_coarser(upannots.size() - i - 1)->get_probability();

      p.probability += QInsideComputer::compute_simple(upannots[i], log_norms[i],
                                                       leftannots[i], rightannots[i],
                                                       rule_probs);

      if(p.probability == -std::numeric_limits<double>::infinity())
        //std::cout << "it's happening! (b)" << std::endl;
        break;
    }

    p.dtrs->set_calculated(true);
    p.dtrs->set_q_score(p.probability);
  }


  //  p.probability += std::exp(p.probability);
  p.probability += left.get_prob_model().get(0).probability + right.get_prob_model().get(0).probability;
  p.probability += e.get_relaxation();

  if (candidates.empty() || p.probability > derivations[0].probability)
    derivations[0] = p;
  candidates.emplace_back(p);
}


void MaxRuleProbabilityMultiple::pick_best_unary(const Edge& e,
                                                 UnaryDaughter & dtr, const std::vector<double>& log_norms)
{
  Edge& left  = dtr.left_daughter();

  if(left.get_prob_model().get(0).dtrs && // should be an assertion
       // only branch on a binary or a lexical : prevent unary chains
       (left.get_prob_model().get(0).dtrs->is_lexical() || left.get_prob_model().get(0).dtrs->is_binary())
       )
  {
    packed_edge_probability p;
    p.dtrs = &dtr;

    UnaryDaughter* d = static_cast<UnaryDaughter*>(p.dtrs);

    if (p.dtrs->is_calculated())
    {
      p.probability = p.dtrs->get_q_score();
    }
    else
    {
      const std::vector<AnnotationInfo>& upannots = get_annotations_backup();
      const std::vector<AnnotationInfo>& leftannots = left.get_prob_model().get_annotations_backup();

      p.probability = 0;

      for (unsigned i = 0; i < upannots.size(); ++i)
      {
        if(log_norms[i] == -std::numeric_limits<double>::infinity())
          continue;
        const auto& rule_probs = d->get_rule()->get_coarser(upannots.size() - i - 1)->get_probability();

        p.probability += QInsideComputer::compute_simple(upannots[i], log_norms[i],
                                                         leftannots[i], rule_probs);

        if(p.probability ==-std::numeric_limits<double>::infinity())
          //std::cout << "it's happening! (u)" << std::endl;
          break;
      }

      p.dtrs->set_calculated(true);
      p.dtrs->set_q_score(p.probability);

    }
    //    p.probability += std::exp(p.probability);
    p.probability +=  left.get_prob_model().get(0).probability;
    p.probability += e.get_relaxation();

    if (candidates.empty() || p.probability > derivations[0].probability)
      derivations[0] = p;

    candidates.emplace_back(p);
  }
}

//read scores and pick best
void MaxRuleProbabilityMultiple::pick_best()
{
  if(!candidates.empty()) {
    if(candidates.size() > size) {
      //      std::cout << candidates.size() << std::endl;
      std::nth_element(candidates.begin(),candidates.begin()+size,candidates.end(), std::greater<packed_edge_probability>());
      candidates.resize(size);

    }
    std::make_heap(candidates.begin(),candidates.end());

    std::pop_heap(candidates.begin(),candidates.end());
    derivations[0] = candidates.back();
    candidates.pop_back();
  }
  else
    derivations.resize(0);

}


void MaxRuleProbabilityMultiple::backup_annotations(const AnnotationInfo& annotations)
{
  // std::cout << "outside " ;
  //   for (int i = 0; i < annotations.outside_probabilities.array.size(); ++i)
  //     {
  //       std::cout << annotations.outside_probabilities.array[i] << " " ;
  //     }
  //   std::cout << std::endl;

  //   std::cout << "inside " ;
  //   for (int i = 0; i < annotations.inside_probabilities.array.size(); ++i)
  //     {
  //       std::cout << annotations.inside_probabilities.array[i] << " " ;
  //     }
  //   std::cout << std::endl;
  annotations_backup.push_back(annotations);
}


inline
std::vector<AnnotationInfo>&
MaxRuleProbabilityMultiple::get_annotations_backup() {return annotations_backup;}

inline
const std::vector<AnnotationInfo>& MaxRuleProbabilityMultiple::get_annotations_backup() const {return annotations_backup;}

///////////////////////////


void MaxRuleProbabilityMultiple::extend_derivation(Edge* edge, unsigned i, bool licence_unaries,
                                                   const std::vector<double>& log_norms)
{
  if(derivations.size() == i) {
    return;
  }

  if(derivations.size() > 0) {

    packed_edge_probability& last = derivations[derivations.size() -1];

    //    std::cout << "last.probability " << last.probability << std::endl;

    assert(last.probability <= 0);
    // if (!(last.probability <=0)) {
    //   std::cout << derivations.size() << std::endl;
    //   std::cout << last.probability << std::endl << std::endl;
    //   last.probability = 0;
    // }

    find_succ(edge,last,licence_unaries, log_norms);
    //    std::cout << "after find_succ" << std::endl;
  }

  if (!candidates.empty()) {

    //get next element from the candidates and append it to derivations
    pop_heap(candidates.begin(),candidates.end());
    derivations.emplace_back(candidates.back());
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


void MaxRuleProbabilityMultiple::find_succ(Edge* edge, packed_edge_probability& pep, bool licence_unaries, const std::vector<double>& log_norms)
{
  if(pep.dtrs->is_lexical())  {
    //std::cout << "find_suc lex" << std::endl;
    return;}
  // binary -> extend left and right daughters


  const std::vector<AnnotationInfo>& upannots = edge->get_prob_model().get_annotations_backup();

  if(pep.dtrs->is_binary()) {

    //    std::cout << "binary case" << std::endl;

    const BinaryDaughter* d = static_cast<const BinaryDaughter*>(pep.dtrs);

    Edge& left  = d->left_daughter();
    const std::vector<AnnotationInfo>& leftannots = left.get_prob_model().get_annotations_backup();

    Edge& right = d->right_daughter();
    const std::vector<AnnotationInfo>& rightannots = right.get_prob_model().get_annotations_backup();

    //extend to the left
    //    std::cout << "bin extending on the left" << std::endl;
    unsigned nextleft = pep.left_index + 1;
    left.extend_derivation(nextleft+1,true, log_norms);

    // we haven't reached the expected number of solutions
    if(nextleft < left.get_prob_model().n_deriv()) {

      packed_edge_probability p(pep);
      p.left_index = nextleft;
      p.probability = 0;


      for (unsigned i = 0; i < upannots.size(); ++i)
        {
          if(log_norms[i] == -std::numeric_limits<double>::infinity())
            continue;
          //          std::cout << *( d->get_rule()->get_coarser(upannots.size() - i + 1)) << std::endl;

          const std::vector<std::vector<std::vector<double> > >& rule_probs =
            d->get_rule()->get_coarser(upannots.size() - i - 1)->get_probability();

          p.probability += QInsideComputer::compute_simple(upannots[i],
                                                           log_norms[i],
                                                           leftannots[i],
                                                           rightannots[i],
                                                           rule_probs);
        }

      p.probability += left.get_prob_model().get(p.left_index).probability + right.get_prob_model().get(p.right_index).probability;


      assert(p.probability <= 0);

      //      std::cout << p.probability << std::endl;

      // TODO : Find a proper way to remove duplicates !
      //      if (std::find_if(candidates.begin(), candidates.end(), test_helper(p)) == candidates.end()) {
        candidates.emplace_back(p);
        push_heap(candidates.begin(), candidates.end());
      // }
      // else
      // {
      //   std::cerr << "FOUND BIN" << std::endl;
      // }

    }

    //extend to the right
    unsigned nextright = pep.right_index + 1;

    right.extend_derivation(nextright+1,true, log_norms);

    if(nextright < right.get_prob_model().n_deriv()) {
      //      std::cout << "bin extending on the right" << std::endl;


      packed_edge_probability p(pep);
      p.right_index = nextright;
      p.probability = 0;


      for (unsigned i = 0; i < upannots.size(); ++i)
        {
          if(log_norms[i] == -std::numeric_limits<double>::infinity())
            continue;
          const std::vector<std::vector<std::vector<double> > >& rule_probs =
            d->get_rule()->get_coarser(upannots.size() - i - 1)->get_probability();

          p.probability += QInsideComputer::compute_simple(upannots[i],
                                                           log_norms[i],
                                                           leftannots[i],
                                                           rightannots[i],
                                                           rule_probs);
        }

      p.probability += left.get_prob_model().get(p.left_index).probability + right.get_prob_model().get(p.right_index).probability;



      assert(p.probability <= 0);

      //      std::cout << p.probability << std::endl;

      //      if(std::find_if(candidates.begin(), candidates.end(), test_helper(p)) == candidates.end()){
        candidates.emplace_back(p);
        push_heap(candidates.begin(), candidates.end());
      // }
      // else
      // {
      //   std::cerr << "FOUND BIN 2" << std::endl;
      // }
    }
  }

  //unary
  else {
    if(!licence_unaries) return;

    //    std::cout << "unary case" << std::endl;

    const UnaryDaughter* d = static_cast<const UnaryDaughter*>(pep.dtrs);

    Edge& left  = d->left_daughter();
    const std::vector<AnnotationInfo>& leftannots = left.get_prob_model().get_annotations_backup();

    //        std::cout << * d->get_rule() << std::endl;


    //extend to the left
    unsigned nextleft = pep.left_index + 1;

    left.extend_derivation(nextleft+1, false, log_norms);

    if(nextleft < left.get_prob_model().n_deriv() ) {
      //        std::cout << "un extending" << std::endl;
      packed_edge_probability p(pep);
      p.left_index = nextleft;
      p.probability = 0;


      for (unsigned i = 0; i < upannots.size(); ++i)
        {
          if(log_norms[i] == -std::numeric_limits<double>::infinity())
            continue;
          const std::vector<std::vector<double> >& rule_probs =
            d->get_rule()->get_coarser(upannots.size() - i - 1)->get_probability();

          p.probability += QInsideComputer::compute_simple(upannots[i],
                                                           log_norms[i],
                                                           leftannots[i],
                                                           rule_probs);
        }

      p.probability +=  left.get_prob_model().get(p.left_index).probability;



      assert(p.probability <= 0);

      //            std::cout << p.probability << std::endl;


      //      if(std::find_if(candidates.begin(), candidates.end(), test_helper(p)) == candidates.end()){
        candidates.push_back(p);
        push_heap(candidates.begin(), candidates.end());
      // }
      // else
      // {
      //   std::cerr << "FOUND UN" << std::endl;
      // }
    }
  }
}



#endif /* _MAXRULEMULTIPLEPROBABILITY_H_ */
