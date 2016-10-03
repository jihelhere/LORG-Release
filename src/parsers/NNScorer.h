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


#include "cnn/cnn.h"
#include "cnn/training.h"
#include "cnn/expr.h"


#include "cnn/rnn.h"
#include "cnn/gru.h"
#include "cnn/lstm.h"

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
  cnn::ComputationGraph * cg;


  //internal rules FF
  cnn::Parameters* _p_W_int;
  cnn::Parameters* _p_b_int;
  cnn::Parameters* _p_o_int;

  //lexical rules FF
  cnn::Parameters* _p_W_lex;
  cnn::Parameters* _p_b_lex;
  cnn::Parameters* _p_o_lex;


  //span FF
  cnn::Parameters* _p_W_span;
  cnn::Parameters* _p_b_span;
  cnn::Parameters* _p_o_span;


  // embeddings for words and symbols
  cnn::LookupParameters* _p_word;
  cnn::LookupParameters* _p_nts;


  cnn::LSTMBuilder l2r_builder;
  cnn::LSTMBuilder r2l_builder;
  std::vector<cnn::expr::Expression> lstm_embeddings;

  std::unordered_map<const Production*, cnn::expr::Expression> rules_expressions;
  //std::unordered_map<std::tuple<int,int,int,int>, cnn::expr::Expression> spans_expressions;
  std::unordered_map<std::tuple<int,int>, cnn::expr::Expression> spans_expressions;
  std::unordered_map<const Edge*, std::vector<cnn::expr::Expression>> edges_expressions;


  std::unordered_set<anchored_binrule_type> anchored_binaries;
  std::unordered_set<anchored_unirule_type> anchored_unaries;
  std::unordered_set<anchored_lexrule_type> anchored_lexicals;
  bool gold;


  std::vector<Word> words;


  nn_scorer(cnn::Model& m);

  void set_cg(cnn::ComputationGraph& g) {cg = &g;};


  void set_gold(std::vector<anchored_binrule_type>& ancbin,
                std::vector<anchored_unirule_type>& ancuni,
                std::vector<anchored_lexrule_type>& anclex);

  void unset_gold();


  void set_words(const std::vector<Word>& w); //todo : a ref/pointer to an existing vector ?



  double compute_lexical_score(int position, const MetaProduction* mp,
                               std::vector<cnn::expr::Expression>& expv);
  double compute_unary_score(int begin, int end, const MetaProduction* mp,
                             std::vector<cnn::expr::Expression>& expv);

  double compute_binary_score(int begin, int end, int mid, const MetaProduction* mp,
                              std::vector<cnn::expr::Expression>& expv);
  double
  compute_internal_span_score(int begin,
                              int end,
                              int medium,
                              int lhs,
                              std::vector<cnn::expr::Expression>& e);

  void register_expression(const Edge*ep, std::vector<cnn::expr::Expression>& expv);

  void clear();
  void precompute_rule_expressions(const std::vector<Rule>& brules,
                                   const std::vector<Rule>& urules);
  void precompute_span_expressions(const std::unordered_set<int>& lhs);

  void precompute_embeddings();


  cnn::expr::Expression rule_expression(int lhs, int rhs0, int rhs1);
  cnn::expr::Expression lexical_rule_expression(int lhs, int word_position);
  cnn::expr::Expression span_expression(int lhs, int word_position);


 private:
  double
  compute_internal_rule_score(const Production* r, std::vector<cnn::expr::Expression>& epp);

};
