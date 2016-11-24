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

#include "rules/Production.h"
#include "rules/Rule.h"

#include "utils/SymbolTable.h"
#include "computation_attachment.h"



inline int nt_triple_to_index(int lhs, int rhs0, int rhs1 = -1)
{
  // should not change after 1st call.
  static int b = SymbolTable::instance_nt().get_symbol_count();
  return (b+1) * (lhs * b + rhs0) + (rhs1 >= 0 ? rhs1 : b);

}


class cfg_rule_representation : public computation_attachment
{
 public:
  cfg_rule_representation(){}; // for static creation

  cfg_rule_representation(bool /*init_global*/,
                          dynet::Model& m,
                          unsigned nt_embedding_size,
                          unsigned hidden_size);

  dynet::expr::Expression get_nt_expr(int lhs);

  double compute_internal_rule_score(const Production* r);

  void precompute_rule_expressions(const std::vector<Rule>& brules,
                                   const std::vector<Rule>& urules,
                                   bool train_mode);

  dynet::expr::Expression& brule_expression(const Production& p);
  dynet::expr::Expression& urule_expression(const Production& p);

private:
  dynet::expr::Expression rule_expression(int lhs, int rhs0, int rhs1);

  dynet::LookupParameter _p_nts;
  std::vector<dynet::expr::Expression> rule_expressions;
  std::vector<double> rule_scores;

  //internal rules FF
  dynet::Parameter _p_W_int;
  dynet::Parameter _p_b_int;
  dynet::Parameter _p_o_int;
};
