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
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
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
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#endif


class Edge;


typedef std::tuple<int,int,int,Production> anchored_binrule_type;
typedef std::tuple<int,int,Production> anchored_unirule_type;
typedef std::tuple<int,Production> anchored_lexrule_type;


struct nn_scorer
{
  dynet::ComputationGraph * cg;

  static bool model_initialized;

  //internal rules FF
  static dynet::Parameter _p_W_int;
  static dynet::Parameter _p_b_int;
  static dynet::Parameter _p_o_int;

  //lexical rules FF
  static dynet::Parameter _p_W_lex;
  static dynet::Parameter _p_b_lex;
  static dynet::Parameter _p_o_lex;

  //span FF
  static dynet::Parameter _p_W_span;
  static dynet::Parameter _p_b_span;
  static dynet::Parameter _p_o_span;

  // embeddings for words and symbols
  static dynet::LookupParameter _p_word;
  static dynet::LookupParameter _p_nts;

  static std::vector<dynet::LSTMBuilder> l2r_builders;
  static std::vector<dynet::LSTMBuilder> r2l_builders;

  std::vector<dynet::expr::Expression> lstm_embeddings;
  std::unordered_map<const Production*, double> rule_scores;
  std::unordered_map<std::tuple<int,int,int>, double> span_scores;

  const std::unordered_set<anchored_binrule_type>* anchored_binaries;
  const std::unordered_set<anchored_unirule_type>* anchored_unaries;
  const std::unordered_set<anchored_lexrule_type>* anchored_lexicals;
  bool gold;
  unsigned lexical_level;
  unsigned span_level;



  const std::vector<Word>* words;

  nn_scorer(dynet::Model& m, unsigned lex_level, unsigned span_level,
            unsigned word_embedding_size,
            unsigned nt_embedding_size,
            unsigned hidden_size,
            unsigned lstm_hidden_size);

  ~nn_scorer() {};

  void set_cg(dynet::ComputationGraph& g) {cg = &g;};


  void set_gold(std::unordered_set<anchored_binrule_type>& ancbin,
                std::unordered_set<anchored_unirule_type>& ancuni,
                std::unordered_set<anchored_lexrule_type>& anclex);

  void unset_gold();


  void set_words(const std::vector<Word>& w);



  double compute_lexical_score(int position, const MetaProduction* mp);
  double compute_unary_score(int begin, int end, const MetaProduction* mp);

  double compute_binary_score(int begin, int end, int mid, const MetaProduction* mp);

  double
  compute_internal_span_score(int begin, int end, int medium, int lhs);

  void clear();
  void precompute_rule_expressions(const std::vector<Rule>& brules,
                                   const std::vector<Rule>& urules);
  void precompute_span_expressions(const std::vector<int>& lhs_int);

  void precompute_embeddings();


  dynet::expr::Expression rule_expression(int lhs, int rhs0, int rhs1);
  dynet::expr::Expression lexical_rule_expression(int lhs, unsigned word_position);
  dynet::expr::Expression span_expression(int lhs, int word_position_start, int word_position_end);


 private:
  double
  compute_internal_rule_score(const Production* r);

};
