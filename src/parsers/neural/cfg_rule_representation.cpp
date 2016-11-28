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
                                                     const std::vector<Rule>& urules)
{
  // first time we initialize structures
  // TODO: do it clean
  static bool init = true;
  if (init)
  {
    auto s = brules.size() + urules.size();
    rule_expressions.resize(s);
    rule_scores.resize(s);

    unsigned i = 0;
    for (; i < brules.size(); ++i)
    {
      auto&& t = nt_triple_to_index(brules[i].get_lhs(),
                                    brules[i].get_rhs0(),
                                    brules[i].get_rhs1());
      rule_to_idx_map[t] = i;
    }

    unsigned j = 0;
    for (; j < urules.size(); ++j)
    {
      auto&& t = nt_triple_to_index(urules[j].get_lhs(),
                                    urules[j].get_rhs0());
      rule_to_idx_map[t] = i + j;
    }

    init = false;
  }

  ////////////////////////////
  static int pad = SymbolTable::instance_nt().get_symbol_count();
  unsigned i = 0;
  for (; i < brules.size(); ++i)
  {
    auto&& r = brules[i];
    auto&& e = rule_expressions[i] = rule_expression(r.get_lhs(), r.get_rhs0(), r.get_rhs1());
    rule_scores[i] = as_scalar(e.value());
  }

  unsigned j = 0;
  for (; j < urules.size(); ++j)
  {
    auto&& r = urules[j];
    auto&& e =rule_expressions[i+j] = rule_expression(r.get_lhs(), r.get_rhs0(), pad);
    rule_scores[i+j] = as_scalar(e.value());
  }
  //  std::cout << "after rule precomputation" << std::endl;
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
    //static int pad = SymbolTable::instance_nt().get_symbol_count();
    return rule_scores[rule_to_idx_map[nt_triple_to_index(r->get_lhs(), r->get_rhs0(), (r->get_rhs().size() > 1 ? r->get_rhs1() : -1))]];
}

dynet::expr::Expression&
cfg_rule_representation::brule_expression(const Production& p)
{
  return
      rule_expressions[rule_to_idx_map[nt_triple_to_index(p.get_lhs(), p.get_rhs0(), p.get_rhs1())]];
}

dynet::expr::Expression&
cfg_rule_representation::urule_expression(const Production& p)
{
  return
      rule_expressions[rule_to_idx_map[nt_triple_to_index(p.get_lhs(), p.get_rhs0())]];
}
