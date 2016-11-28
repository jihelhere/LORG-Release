// -*- mode: c++ -*-
#pragma once

#include "rules/MetaProduction.h"
#include "rules/Production.h"
#include "rules/Rule.h"
#include <unordered_map>
#include "utils/hash_impl.h"
#include "Word.h"


#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winconsistent-missing-override"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weffc++"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-local-typedef"
#else
// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif


#include "dynet/dynet.h"
#include "dynet/training.h"
#include "dynet/expr.h"


#include "dynet/rnn.h"
#include "dynet/gru.h"
#include "dynet/lstm.h"

#if defined(__clang__)
#pragma clang diagnostic pop
#pragma clang diagnostic pop
#pragma clang diagnostic pop
#pragma clang diagnostic pop
#else
//#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#endif

#include "computation_attachment.h"
#include "cfg_rule_representation.h"
#include "word_representation.h"
#include "lexicon_representation.h"
#include "lexical_rule_representation.h"
#include "span_representation.h"

typedef std::tuple<int,int,int,Production> anchored_binrule_type;
typedef std::tuple<int,int,Production> anchored_unirule_type;
typedef std::tuple<int,Production> anchored_lexrule_type;




struct nn_scorer : public computation_attachment
{
  static bool train_mode;

  static cfg_rule_representation cfg;

  word_representation * word_repr;
  lexicon_representation * lex_repr;
  lexical_rule_representation * lexrule_repr;
  span_representation * span_repr;




  const std::unordered_set<anchored_binrule_type>* anchored_binaries;
  const std::unordered_set<anchored_unirule_type>* anchored_unaries;
  const std::unordered_set<anchored_lexrule_type>* anchored_lexicals;

  bool gold;
  unsigned lexical_level;
  unsigned span_level;

  const std::vector<Word>* words;

  bool use_char_emb;
  bool use_span_midpoints;

  nn_scorer(bool init_global,
            dynet::Model& m, unsigned lex_level, unsigned span_level,
            unsigned word_embedding_size,
            unsigned nt_embedding_size,
            unsigned hidden_size,
            unsigned lstm_hidden_size,
            bool char_emb,
            bool span_midp);

  ~nn_scorer() {};

  static void set_cg(dynet::ComputationGraph& g) {cg = &g;};


  void set_gold(std::unordered_set<anchored_binrule_type>& ancbin,
                std::unordered_set<anchored_unirule_type>& ancuni,
                std::unordered_set<anchored_lexrule_type>& anclex);

  void unset_gold();


  void set_words(const std::vector<Word>& w);

  double compute_lexical_score(int position, const MetaProduction* mp);
  double compute_unary_score(int begin, int end, const MetaProduction* mp);

  double compute_binary_score(int begin, int end, int mid, const MetaProduction* mp);

  double
  compute_internal_span_score(int begin, int end, int medium, int lhs, int rhs0);

  void clear();

  void precompute_rule_expressions(const std::vector<Rule>& brules,
                                   const std::vector<Rule>& urules);

  void precompute_span_expressions(const std::vector<int>& lhs_int,
                                   const std::vector<int>& rhs0_int,
                                   const std::vector<int>& rhs1_int);

  void precompute_embeddings();


  dynet::expr::Expression& span_expression(int lhs, int word_position_start, int word_position_end, int word_medium);
  dynet::expr::Expression& span_init(int lhs, int begin);
  dynet::expr::Expression& span_end(int lhs, int end);
  dynet::expr::Expression& span_split(int lhs, int split);

  // dynet::expr::Expression& span_init_un(int lhs, int begin);
  // dynet::expr::Expression& span_end_un(int lhs, int end);

  dynet::expr::Expression& span_rhs0_init(int rhs0, int begin);
  dynet::expr::Expression& span_rhs0_end(int rhs0, int end);
  dynet::expr::Expression& span_rhs0_split(int rhs0, int split);


  void set_dropout(float d);
  void unset_dropout();

 private:
  double
  compute_internal_rule_score(const Production* r);
};
