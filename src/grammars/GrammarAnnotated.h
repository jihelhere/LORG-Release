// -*- mode: c++ -*-
#pragma once


#include "Grammar.h"


#include <unordered_map>
#include "AnnotatedLabelsInfo.h"
#include "utils/SymbolTable.h"

#include <numeric>


#include "rules/LexicalRuleC2f.h"
#include "rules/URuleC2f.h"
#include "rules/BRuleC2f.h"


#include "utils/Tree.h"
#include "utils/hash_impl.h"


#define uomap std::unordered_map

#include <boost/serialization/map.hpp>
#include "utils/unordered_map_serialization.hpp"


typedef std::pair< int, unsigned> asymb;
typedef uomap<asymb, uomap< asymb, asymb> > PathMatrix;

template <typename Bin, typename Un, typename Lex>
class GrammarAnnotated : public Grammar<Bin,Un,Lex>, public AnnotatedContents
{
public:
  GrammarAnnotated();
  virtual ~GrammarAnnotated() {};

  void init();
  void set_logmode();
  void remove_unlikely_annotations_all_rules(double threshold);

  inline const PathMatrix&
  get_unary_decoding_paths() const {return viterbi_decoding_paths;}

  inline PathMatrix&
  get_unary_decoding_paths() {return viterbi_decoding_paths;}

  GrammarAnnotated(const std::string& filename);

  GrammarAnnotated * create_projection(const std::vector<std::vector<double> >& expected_counts,
                                       const std::vector<std::vector<std::vector<unsigned> > >& annotation_mapping) const;

  void
  compute_transition_probabilities(uomap<int, uomap<unsigned, uomap<int, uomap<unsigned, double > > > >&) const;

  std::vector<double> compute_priors() const;

  void linear_smooth(double internal_threshold, double lexical_threshold);


  inline const std::map< short, Tree<unsigned> >&  get_history_trees() const {return history_trees;}
  inline std::map< short, Tree<unsigned> >&  get_history_trees() {return history_trees;}


  // TODO: private + accessors
  std::map<std::string, std::string>  gram_conf;


private:
    /**
     \brief create new unary rules encoding unary chains
     \param decoding_paths the matrix of best paths for decoding after parsing
     \param max_path_length the maximum length of chains to encode
     \param unaries the unary_rules to enrich with unary chains (modified)
   */
  void compute_unary_chains(PathMatrix& decoding_paths,
                            unsigned max_path_length,
                            std::vector<Un>& unaries);

  GrammarAnnotated (const std::vector<std::vector<double> >& conditional_probabilities,
                    const std::vector<std::vector<std::vector<unsigned> > >& mapping,
                    const std::vector<Bin>& old_binary_rules,
                    const std::vector<Un>& old_unary_rules,
                    const std::vector<Lex>& old_lexical_rules);


  // for unary chains encoding/decoding
  PathMatrix viterbi_decoding_paths;

  // history of annotations (trees of numbers)
  std::map<short, Tree<unsigned> > history_trees;


 private:
  friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/)
    {
      ar & boost::serialization::base_object<Grammar<Bin,Un,Lex>>(*this);
      ar & boost::serialization::base_object<AnnotatedContents>(*this);
      ar & gram_conf;
      ar & viterbi_decoding_paths;
      ar & history_trees;
    }



};





/**************************************************************************************
 * utils
 *
 * functions. Defined in GrammarAnnotated.cpp
 **************************************************************************************/
void calculate_expected_counts(const uomap<int, uomap<unsigned, uomap<int, uomap<unsigned, double > > > >& trans,
                               const AnnotatedLabelsInfo& ali,
                               std::vector<std::vector<double> >& result);


std::vector<std::vector<std::vector<unsigned> > > compute_mapping(unsigned from,
                                                                  unsigned to,
                                                                  const std::vector< std::vector<std::vector< std::vector<unsigned>>>> & annot_descendants);

std::vector<std::vector<double> >
calculate_conditional_probs(const std::vector<std::vector<double> >& expected_counts,
                            const std::vector<std::vector<std::vector<unsigned> > >& mapping);


/**************************************************************************************
 * struct project_rule used in GrammarAnnotated.hpp. Defined in GrammarAnnotated.cpp
 **************************************************************************************/
class project_rule
{
public:
  project_rule(const std::vector<std::vector<double> >& conditional_probabilities_,
               const std::vector<std::vector<std::vector<unsigned> > >& mapping_)
  :
  conditional_probabilities(conditional_probabilities_),
  mapping(mapping_),
  inverted(invert_mapping(mapping_))
  {}

  LexicalRuleC2f operator()(const LexicalRuleC2f& old_rule) const;
  URuleC2f operator()(const URuleC2f& old_rule) const;
  BRuleC2f operator()(const BRuleC2f& old_rule) const;

private:
  const std::vector<std::vector<double> >& conditional_probabilities;
  const std::vector<std::vector<std::vector<unsigned> > >& mapping;
  std::vector<uomap<unsigned,unsigned> > inverted;


  //TODO: iterate on new_probs instead of old_probs
  // would be more efficient (especially for binary rules)

  std::vector<uomap<unsigned,unsigned> > invert_mapping(std::vector<std::vector<std::vector<unsigned> > > mapping);
};
