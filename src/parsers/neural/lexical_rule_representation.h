// -*- mode: c++ -*-
#pragma once



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

#include "lexicon_representation.h"

#include "utils/hash_impl.h"

class lexical_rule_representation
{
 public:
  lexical_rule_representation(lexicon_representation * l,
                              cfg_rule_representation * c)
      : lr(l), cfg(c) {};
  virtual ~lexical_rule_representation() {};
  virtual
  dynet::expr::Expression lexical_rule_expression(int lhs, unsigned word_idx)=0;

  dynet::expr::Expression& retrieve_lexical_rule_expression(int lhs, unsigned word_idx);

  void clear();

 protected:
  lexicon_representation * lr;
  cfg_rule_representation * cfg;

  // (lhs,word_position) -> expression
  std::unordered_map<std::tuple<int,int>,dynet::expr::Expression> expression_map;

  static dynet::Parameter _p_W_lex;
  static dynet::Parameter _p_b_lex;
  static dynet::Parameter _p_o_lex;
};

class simple_lexical_rule : public lexical_rule_representation,
                            public computation_attachment
{
 public:
  simple_lexical_rule(bool init_global,
                      dynet::Model& m,
                      unsigned nt_embedding_size,
                      unsigned word_embedding_size,
                      unsigned hidden_size,
                      lexicon_representation * l,
                      cfg_rule_representation * c);

  dynet::expr::Expression lexical_rule_expression(int lhs, unsigned word_idx);
};


class bilstm_lexical_rule : public lexical_rule_representation,
                            public computation_attachment
{
 public:

  bilstm_lexical_rule(bool init_global,
                      dynet::Model& m,
                      unsigned nt_embedding_size,
                      unsigned hidden_size,
                      unsigned lstm_hidden_size,
                      lexicon_representation * l,
                      cfg_rule_representation * c);

  dynet::expr::Expression lexical_rule_expression(int lhs, unsigned word_idx);

  void set_dropout(float d);
  void unset_dropout();
};
