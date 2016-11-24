#include "cfg_rule_representation.h"

namespace dy = dynet;
namespace de = dy::expr;



cfg_rule_representation::cfg_rule_representation(bool /*init_global*/,
                                                 dynet::Model& m,
                                                 unsigned nt_embedding_size,
                                                 unsigned hidden_size)
    : computation_attachment()
{
  unsigned s = SymbolTable::instance_nt().get_symbol_count();

  rule_expressions.resize(s * s * (s+1));
  rule_scores.resize(s * s * (s+1));

  _p_nts = m.add_lookup_parameters(s+1, {nt_embedding_size});

  // grammar rules
  _p_W_int = m.add_parameters({hidden_size, nt_embedding_size*3});
  _p_b_int = m.add_parameters({hidden_size});
  _p_o_int = m.add_parameters({1,hidden_size});
}

dynet::expr::Expression cfg_rule_representation::get_nt_expr(int lhs)
{
  return dynet::expr::lookup(*cg, _p_nts, lhs);
}

void
cfg_rule_representation::precompute_rule_expressions(const std::vector<Rule>& brules,
                                                     const std::vector<Rule>& urules,
                                                     bool train_mode)
{
  std::cout << "before rule precomputation" << std::endl;

  static int pad = SymbolTable::instance_nt().get_symbol_count();

  for (const auto& r : brules)
  {
    auto&& e = rule_expression(r.get_lhs(), r.get_rhs0(), r.get_rhs1());
    auto&& t = nt_triple_to_index(r.get_lhs(),r.get_rhs0(),r.get_rhs1());
    if (train_mode) rule_expressions[t] = e;
    rule_scores[t] = as_scalar(cg->get_value(e.i));
  }

  std::cout << " after rule precomputation for binaries " << std::endl;

  for (const auto& r : urules)
  {
    auto&& e =  rule_expression(r.get_lhs(), r.get_rhs0(), pad);
    auto&& t = nt_triple_to_index(r.get_lhs(),r.get_rhs0());
    if (train_mode) rule_expressions[t] = e;
    rule_scores[t] = as_scalar(cg->get_value(e.i));
  }

  std::cout << "after rule precomputation" << std::endl;
}


de::Expression cfg_rule_representation::rule_expression(int lhs, int rhs0, int rhs1)
{
  std::lock_guard<std::mutex> guard(computation_attachment::cg_mutex);

  auto&& i = de::concatenate({de::lookup(*cg,_p_nts,rhs0),
                              de::lookup(*cg,_p_nts,rhs1),
                              de::lookup(*cg,_p_nts,lhs)});

  auto&& W = de::parameter(*cg, _p_W_int);
  auto&& b = de::parameter(*cg, _p_b_int);
  auto&& o = de::parameter(*cg, _p_o_int);
  return o * de::rectify(W*i + b);
}


double cfg_rule_representation::compute_internal_rule_score(const Production* r)
  {
    static int pad = SymbolTable::instance_nt().get_symbol_count();
    return rule_scores[nt_triple_to_index(r->get_lhs(), r->get_rhs0(), (r->get_rhs().size() > 1 ? r->get_rhs1() : pad))];
}



dynet::expr::Expression&
cfg_rule_representation::brule_expression(const Production& p)
{
  return
      rule_expressions[nt_triple_to_index(p.get_lhs(), p.get_rhs0(), p.get_rhs1())];
}

dynet::expr::Expression&
cfg_rule_representation::urule_expression(const Production& p)
{
  return
      rule_expressions[nt_triple_to_index(p.get_lhs(), p.get_rhs0(), -1)];
}
