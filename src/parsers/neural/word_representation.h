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
#include "dynet/expr.h"

#include "dynet/rnn.h"
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


#include "Word.h"
#include "computation_attachment.h"


class word_representation
{
 public:
  virtual ~word_representation() {};
  virtual dynet::expr::Expression word_reprentation(const Word& w)=0;
  virtual dynet::expr::Expression get_pad_expression()=0;
};

class simple_word_representation : public word_representation,
                                   public computation_attachment
{
 public:
  simple_word_representation(bool init_global, dynet::Model& m,
                             unsigned word_embedding_size);
  virtual ~simple_word_representation() {};

  dynet::expr::Expression word_reprentation(const Word& w);

  dynet::expr::Expression get_pad_expression();

 private:
  static dynet::LookupParameter _p_word;
};


class char_word_representation : public word_representation,
                                 public computation_attachment

{
 public:
  char_word_representation(bool init_global,
                           dynet::Model& m,
                           unsigned word_embedding_size);

  virtual ~char_word_representation() {};

  dynet::expr::Expression word_reprentation(const Word& w);
  dynet::expr::Expression get_pad_expression();

 private:
  static dynet::LookupParameter _p_word;

  static dynet::LSTMBuilder letter_l2r_builder;
  static dynet::LSTMBuilder letter_r2l_builder;
};
