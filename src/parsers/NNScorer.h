// -*- mode: c++ -*-
#pragma once

#include "rules/MetaProduction.h"
#include "rules/Production.h"
#include <unordered_map>
#include "utils/hash_impl.h"

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


#include <cnn/cnn.h>
#include <cnn/training.h>
#include "cnn/expr.h"


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

  cnn::Parameter _p_W_int;
  cnn::Parameter _p_b_int;
  cnn::Parameter _p_o_int;


  cnn::Parameter _p_W_lex;
  cnn::Parameter _p_b_lex;
  cnn::Parameter _p_o_lex;


  cnn::LookupParameter _p_word;
  cnn::LookupParameter _p_nts;


  std::unordered_map<const Production*, cnn::expr::Expression> rules_expressions;
  std::unordered_map<const Edge*, cnn::expr::Expression> expressions;

  std::unordered_set<anchored_binrule_type> anchored_binaries;
  std::unordered_set<anchored_unirule_type> anchored_unaries;
  std::unordered_set<anchored_lexrule_type> anchored_lexicals;
  bool gold;


  nn_scorer(cnn::Model& m);

  void set_cg(cnn::ComputationGraph& g) {cg = &g;};


  void set_gold(std::vector<anchored_binrule_type>& ancbin,
                std::vector<anchored_unirule_type>& ancuni,
                std::vector<anchored_lexrule_type>& anclex);

  void unset_gold();




  double compute_lexical_score(int position, const MetaProduction* mp,
                                cnn::expr::Expression* expp);
  double compute_unary_score(int begin, int end, const MetaProduction* mp,
                              cnn::expr::Expression* expp);
  double compute_binary_score(int begin, int end, int mid, const MetaProduction* mp,
                               cnn::expr::Expression* expp);


  void register_expression(const Edge*ep, const cnn::expr::Expression& expp);

  void clear();
};
