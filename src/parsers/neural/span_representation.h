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


#include <tuple>
#include <unordered_map>
#include "utils/hash_impl.h"

#include "lexicon_representation.h"
#include "cfg_rule_representation.h"

inline bool
is_artificial(int lhs)
{
  return SymbolTable::instance_nt().translate(lhs)[0] == '[' ;
}

class span_representation
{
 public:
  span_representation() {};
  virtual ~span_representation() {};

  virtual
  void precompute_span_expressions(const std::vector<int>&,
                                   const std::vector<int>&,
                                   const std::vector<int>&,
                                   const std::vector<Word>&,
                                   bool) {};
  virtual void clear() {};

  virtual double get_span_score_lhs_begin(int/*lhs*/,int /*begin*/)=0;
  virtual double get_span_score_lhs_end(int/*lhs*/, int/*end*/)=0;
  virtual double get_span_score_lhs_begin_unary(int/*lhs*/,int /*begin*/)=0;
  virtual double get_span_score_lhs_end_unary(int/*lhs*/, int/*end*/)=0;

  virtual double get_span_score_lhs_split(int/*lhs*/, int/*split*/)=0;

  virtual double get_span_score_rhs0_begin(int/*rhs*/,int /*begin*/)=0;
  virtual double get_span_score_rhs0_end(int/*rhs*/, int/*end*/)=0;
  virtual double get_span_score_rhs0_split(int/*rhs*/, int/*split*/)=0;

  virtual double get_span_score_bin_info(int, // begin
                                         int, // end
                                         int, // split
                                         int // root_info
                                         )=0;

  virtual double get_span_score_una_info(int, // begin
                                         int, // end
                                         int // root_info
                                         )=0;

  virtual dynet::expr::Expression& get_span_expr_lhs_init(int, int)=0;
  virtual dynet::expr::Expression& get_span_expr_lhs_init_unary(int, int)=0;

  virtual  dynet::expr::Expression& get_span_expr_lhs_end (int, int)=0;
  virtual  dynet::expr::Expression& get_span_expr_lhs_end_unary(int, int)=0;

  virtual dynet::expr::Expression& get_span_expr_lhs_split(int, int)=0;

  virtual dynet::expr::Expression& get_span_expr_rhs0_init(int, int)=0;
  // virtual dynet::expr::Expression& get_span_expr_rhs0_init_unary(int, int)=0;
  virtual dynet::expr::Expression& get_span_expr_rhs0_end(int, int)=0;
  // virtual dynet::expr::Expression& get_span_expr_rhs0_end_unary(int, int)=0;

  virtual dynet::expr::Expression& get_span_expr_rhs0_split(int, int)=0;

  virtual dynet::expr::Expression& get_span_expr_lhs_info(int, int, int, int)=0;

};


class empty_span_representation : public span_representation
{
 public:
  empty_span_representation() {};
  virtual ~empty_span_representation() {};


  virtual double get_span_score_lhs_begin(int/*lhs*/,int /*begin*/) override {return 0.0;}
  virtual double get_span_score_lhs_end(int/*lhs*/, int/*end*/) override    {return 0.0;}
  virtual double get_span_score_lhs_begin_unary(int,int) override {return 0.0;}
  virtual double get_span_score_lhs_end_unary(int, int)  override {return 0.0;}


  virtual double get_span_score_lhs_split(int/*lhs*/, int/*split*/) override {return 0.0;}

  virtual double get_span_score_rhs0_begin(int/*rhs*/,int /*begin*/) override {return 0.0;}
  virtual double get_span_score_rhs0_end(int/*rhs*/, int/*end*/)  override {return 0.0;}
  virtual double get_span_score_rhs0_split(int/*rhs*/, int/*split*/) override {return 0.0;}

  virtual double get_span_score_bin_info(int, // begin
                                         int, // end
                                         int, // split
                                         int // root_info
                                         ) override {return 0.0;};

  virtual double get_span_score_una_info(int, // begin
                                         int, // end
                                         int // root_info
                                         ) override {return 0.0;};


  virtual dynet::expr::Expression& get_span_expr_lhs_init(int, int) override
  {
    throw std::logic_error("should not be called");
  }
  virtual dynet::expr::Expression& get_span_expr_lhs_init_unary(int, int) override
  {
    throw std::logic_error("should not be called");
  }


  virtual  dynet::expr::Expression& get_span_expr_lhs_end (int, int) override
  {
    throw std::logic_error("should not be called");
  }
  virtual  dynet::expr::Expression& get_span_expr_lhs_end_unary(int, int) override
  {
    throw std::logic_error("should not be called");
  }

  virtual dynet::expr::Expression& get_span_expr_lhs_split(int, int) override
  {
    throw std::logic_error("should not be called");
  }


  virtual dynet::expr::Expression& get_span_expr_rhs0_init(int, int) override
  {
    throw std::logic_error("should not be called");
  }


  virtual dynet::expr::Expression& get_span_expr_rhs0_end(int, int) override
  {
    throw std::logic_error("should not be called");
  }


  virtual dynet::expr::Expression& get_span_expr_rhs0_split(int, int) override
  {
    throw std::logic_error("should not be called");
  }


  virtual dynet::expr::Expression& get_span_expr_lhs_info(int, int, int, int) override
  {
    throw std::logic_error("should not be called");
  }

};


class all_span_representation : public span_representation,
                                public computation_attachment
{
 public:
  all_span_representation(bool init_global,
                          dynet::Model& m,
                          unsigned span_level, // should be an enum?
                          unsigned input_size,
                          unsigned nt_embedding_size,
                          unsigned hidden_size,
                          bool span_midpoints,
                          lexicon_representation * l,
                          cfg_rule_representation * c);

  virtual ~all_span_representation() {};

  void precompute_span_expressions(const std::vector<int>& lhs_int,
                                   const std::vector<int>& rhs0_int,
                                   const std::vector<int>& rhs1_int,
                                   const std::vector<Word>& words,
                                   bool train_mode);

  void clear();

  double get_span_score_lhs_begin(int lhs, int begin);
  double get_span_score_lhs_begin_unary(int lhs, int begin);
  double get_span_score_lhs_end(int lhs, int end);
  double get_span_score_lhs_end_unary(int lhs, int end);
  double get_span_score_lhs_split(int lhs, int split);

  double get_span_score_rhs0_begin(int rhs,int begin);
  double get_span_score_rhs0_end(int rhs,int end);
  double get_span_score_rhs0_split(int rhs,int split);


  double get_span_score_bin_info(int begin, int end, int split, int root_info);

  double get_span_score_una_info(int begin, int end, int root_info);

  dynet::expr::Expression& get_span_expr_lhs_init(int lhs, int begin);
  dynet::expr::Expression& get_span_expr_lhs_end (int lhs, int end);
  dynet::expr::Expression& get_span_expr_lhs_init_unary(int lhs, int begin);
  dynet::expr::Expression& get_span_expr_lhs_end_unary (int lhs, int end);

  dynet::expr::Expression& get_span_expr_lhs_split(int lhs, int split);

  dynet::expr::Expression& get_span_expr_rhs0_init(int rhs0, int begin);
  dynet::expr::Expression& get_span_expr_rhs0_end(int rhs0, int end);
  dynet::expr::Expression& get_span_expr_rhs0_split(int rhs0, int split);

  dynet::expr::Expression& get_span_expr_lhs_info(int lhs, int begin, int end, int medium);

 private:
  lexicon_representation * lr;
  cfg_rule_representation * cfg;

  unsigned span_level;
  bool use_span_midpoints;


  //span FF
  static dynet::Parameter _p_W_span_init;
  // static dynet::Parameter _p_W_span_end;
  // static dynet::Parameter _p_W_span_split;

  static dynet::Parameter _p_b_span_init;
  static dynet::Parameter _p_o_span_init;
  static dynet::Parameter _p_b_span_init_un;
  static dynet::Parameter _p_o_span_init_un;

  static dynet::Parameter _p_b_span_end;
  static dynet::Parameter _p_o_span_end;
  static dynet::Parameter _p_b_span_end_un;
  static dynet::Parameter _p_o_span_end_un;

  static dynet::Parameter _p_b_span_split;
  static dynet::Parameter _p_o_span_split;


  static dynet::Parameter _p_W_span0_init;
  static dynet::Parameter _p_b_span0_init;
  static dynet::Parameter _p_o_span0_init;
  static dynet::Parameter _p_b_span0_end;
  static dynet::Parameter _p_o_span0_end;
  static dynet::Parameter _p_b_span0_split;
  static dynet::Parameter _p_o_span0_split;


  static dynet::Parameter _p_W_span_left;
  static dynet::Parameter _p_W_span_right;
  static dynet::Parameter _p_W_span_mid;

  static dynet::Parameter _p_W_span_distance;
  static dynet::Parameter _p_W_span_extra;

  static dynet::Parameter _p_b_span_bin;
  static dynet::Parameter _p_o_span_bin;

  static dynet::Parameter _p_b_span_un;
  static dynet::Parameter _p_o_span_un;

  std::unordered_map<std::tuple<int,int,int,int>, dynet::expr::Expression> span_expressions_bin;
  std::unordered_map<std::tuple<int,int,int>, dynet::expr::Expression> span_expressions_un;
  std::unordered_map<std::tuple<int,int,int,int>, double> span_scores_bin;
  std::unordered_map<std::tuple<int,int,int>, double> span_scores_un;

  std::vector<std::vector<dynet::expr::Expression>> span_expressions_init;
  std::vector<std::vector<double>> span_scores_init;

  std::vector<std::vector<dynet::expr::Expression>> span_expressions_init_un;
  std::vector<std::vector<double>> span_scores_init_un;

  std::vector<std::vector<dynet::expr::Expression>> span_expressions_end;
  std::vector<std::vector<double>> span_scores_end;

  std::vector<std::vector<dynet::expr::Expression>> span_expressions_end_un;
  std::vector<std::vector<double>> span_scores_end_un;

  std::vector<std::vector<dynet::expr::Expression>> span_expressions_split;
  std::vector<std::vector<double>> span_scores_split;



  std::vector<std::vector<dynet::expr::Expression>> span_expressions_rhs0_init;
  std::vector<std::vector<double>> span_scores_rhs0_init;

  std::vector<std::vector<dynet::expr::Expression>> span_expressions_rhs0_end;
  std::vector<std::vector<double>> span_scores_rhs0_end;

  std::vector<std::vector<dynet::expr::Expression>> span_expressions_rhs0_split;
  std::vector<std::vector<double>> span_scores_rhs0_split;



};
