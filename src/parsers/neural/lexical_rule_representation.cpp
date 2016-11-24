#include "lexical_rule_representation.h"

#include <unordered_map>
#include "utils/hash_impl.h"


dynet::Parameter lexical_rule_representation::_p_W_lex;
dynet::Parameter lexical_rule_representation::_p_b_lex;
dynet::Parameter lexical_rule_representation::_p_o_lex;



void lexical_rule_representation::clear()
{
  expression_map.clear();
}

dynet::expr::Expression&
lexical_rule_representation::retrieve_lexical_rule_expression(int lhs, unsigned word_idx)
{
  return expression_map[std::make_tuple(lhs,word_idx)];
}


simple_lexical_rule::simple_lexical_rule(bool init_global,
                                         dynet::Model& m,
                                         unsigned nt_embedding_size,
                                         unsigned word_embedding_size,
                                         unsigned hidden_size,
                                         lexicon_representation * l,
                                         cfg_rule_representation * c)
    : lexical_rule_representation(l,c), computation_attachment()
{
  if (init_global)
  {
    _p_W_lex = m.add_parameters({hidden_size, nt_embedding_size + word_embedding_size * 3});
    _p_b_lex = m.add_parameters({hidden_size});
    _p_o_lex = m.add_parameters({1,hidden_size});
  }
}

dynet::expr::Expression simple_lexical_rule::lexical_rule_expression(
    int lhs, unsigned word_idx)
{
  auto&& t = std::make_tuple(lhs,word_idx);

  std::lock_guard<std::mutex> guard(computation_attachment::cg_mutex);

  auto&& W = dynet::expr::parameter(*cg, _p_W_lex);
  auto&& b = dynet::expr::parameter(*cg, _p_b_lex);
  auto&& o = dynet::expr::parameter(*cg, _p_o_lex);

  auto&& emb = lr->get_embeddings();


  auto&& i =
      dynet::expr::concatenate({emb[word_idx],
                                cfg->get_nt_expr(lhs),
                                word_idx == 0 ? lr->get_word_representation()->get_pad_expression() : emb[word_idx-1],
                                word_idx == emb.size() - 1 ? lr->get_word_representation()->get_pad_expression() : emb[word_idx+1]
        })
      ;

  return expression_map[t] = o * dynet::expr::rectify(W*i + b);
}



bilstm_lexical_rule::bilstm_lexical_rule(bool init_global,
                                         dynet::Model& m,
                                         unsigned nt_embedding_size,
                                         unsigned hidden_size,
                                         unsigned lstm_hidden_size,
                                         lexicon_representation * l,
                                         cfg_rule_representation * c)
  : lexical_rule_representation(l,c), computation_attachment()
  {
    if (init_global)
  {
    _p_W_lex = m.add_parameters({hidden_size, nt_embedding_size + lstm_hidden_size});
    _p_b_lex = m.add_parameters({hidden_size});
    _p_o_lex = m.add_parameters({1,hidden_size});
  }
}




dynet::expr::Expression
bilstm_lexical_rule::lexical_rule_expression(int lhs, unsigned word_idx)
{
  auto&& t = std::make_tuple(lhs,word_idx);
  auto&& emb = lr->get_embeddings();

  std::lock_guard<std::mutex> guard(computation_attachment::cg_mutex);

  auto&& W = dynet::expr::parameter(*cg, _p_W_lex);
  auto&& b = dynet::expr::parameter(*cg, _p_b_lex);
  auto&& o = dynet::expr::parameter(*cg, _p_o_lex);

  auto&& i = dynet::expr::concatenate({emb[word_idx],
                                       cfg->get_nt_expr(lhs)
    });

  return expression_map[t] = o * dynet::expr::rectify(W*i + b);

}
