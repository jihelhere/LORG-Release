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

#include "word_representation.h"

class lexicon_representation
{
 public:
  lexicon_representation(word_representation* wordrep): embeddings(), wr(wordrep) {};
  virtual ~lexicon_representation() {};
  std::vector<dynet::expr::Expression>& get_embeddings() {return embeddings;};
  virtual void precompute_embeddings(const std::vector<Word>& words)=0;

  word_representation * get_word_representation() {return wr;};

  virtual void set_dropout(float) {};
  virtual void unset_dropout() {};

 protected:
  std::vector<dynet::expr::Expression> embeddings;
  word_representation* wr;
};


class simple_lexicon_representation : public lexicon_representation,
                                      public computation_attachment
{
 public:
  simple_lexicon_representation(word_representation* word_repr)
      : lexicon_representation(word_repr), computation_attachment()
  {};

  virtual ~simple_lexicon_representation() {};

  void precompute_embeddings(const std::vector<Word>& words);
};


class bilstm_lexicon_representation : public lexicon_representation,
                                      public computation_attachment
{
 public:
  bilstm_lexicon_representation(word_representation* word_repr,
                                bool init_global,
                                dynet::Model& m,
                                unsigned word_embedding_size,
                                unsigned lstm_hidden_size,
                                unsigned l);

  void precompute_embeddings(const std::vector<Word>& words);

  virtual ~bilstm_lexicon_representation() {};

  void set_dropout(float);
  void unset_dropout();

private:
  unsigned layers;

  static std::vector<dynet::LSTMBuilder> word_l2r_builders;
  static std::vector<dynet::LSTMBuilder> word_r2l_builders;

};
