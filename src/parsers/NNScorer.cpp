#include "NNScorer.h"

#include "utils/SymbolTable.h"
#include "rules/Production.h"

#include "utils/hash_impl.h"


#define WORD_EMBEDDING_SIZE 50
#define NT_EMBEDDING_SIZE 20
#define HIDDEN_SIZE 80

#define LSTM_HIDDEN_SIZE 100

nn_scorer::nn_scorer(dynet::Model& m) :
    cg(nullptr),

    _p_W_int(m.add_parameters({HIDDEN_SIZE, NT_EMBEDDING_SIZE*3})),
    _p_b_int(m.add_parameters({HIDDEN_SIZE})),
    _p_o_int(m.add_parameters({1,HIDDEN_SIZE})),

    _p_W_lex(m.add_parameters({HIDDEN_SIZE, NT_EMBEDDING_SIZE+2*LSTM_HIDDEN_SIZE})),
    _p_b_lex(m.add_parameters({HIDDEN_SIZE})),
    _p_o_lex(m.add_parameters({1,HIDDEN_SIZE})),

    // linear classifier with 2 simple features (POS and WORD, POS)
    // _p_W_lex(m.add_parameters({1, SymbolTable::instance_nt().get_symbol_count() * SymbolTable::instance_word().get_symbol_count()
    //                               + SymbolTable::instance_nt().get_symbol_count()})),


    _p_Wleft_span(m.add_parameters({HIDDEN_SIZE,  2*LSTM_HIDDEN_SIZE + 1})),
    _p_Wright_span(m.add_parameters({HIDDEN_SIZE, 2*LSTM_HIDDEN_SIZE + 1})),
    _p_b_span(m.add_parameters({HIDDEN_SIZE})),
    _p_o_span(m.add_parameters({1,HIDDEN_SIZE})),


    _p_word(m.add_lookup_parameters(SymbolTable::instance_word().get_symbol_count()+1,
                                    {WORD_EMBEDDING_SIZE})),
    _p_nts(m.add_lookup_parameters(SymbolTable::instance_nt().get_symbol_count()+1,
                                   {NT_EMBEDDING_SIZE})),


    l2r_builder(2, WORD_EMBEDDING_SIZE, LSTM_HIDDEN_SIZE, &m),
    r2l_builder(2, WORD_EMBEDDING_SIZE, LSTM_HIDDEN_SIZE, &m),



    rules_expressions(),
    spans_expressions(),

    words(nullptr)

{}


void nn_scorer::clear()
{
  rules_expressions.clear();
  spans_expressions.clear();

  // other members are managed somewhere else!
}

void nn_scorer::set_gold(std::unordered_set<anchored_binrule_type>& ancbin,
                         std::unordered_set<anchored_unirule_type>& ancuni,
                         std::unordered_set<anchored_lexrule_type>& anclex)
{
  anchored_binaries = &ancbin;
  anchored_unaries =  &ancuni;
  anchored_lexicals= &anclex;

  gold = true;
}

void nn_scorer::unset_gold()
{
  anchored_binaries = nullptr;
  anchored_unaries = nullptr;
  anchored_lexicals = nullptr;
  gold = false;
}



double nn_scorer::compute_lexical_score(int position, const MetaProduction* mp)
{
    auto r = static_cast<const Production*>(mp);

    //std::cerr << *r << std::endl;

    //double v = rules_expressions.at(r).second;

    double v = as_scalar(cg->get_value(lexical_rule_expression(mp->get_lhs(), position).i));

    //std::cerr << "after" << std::endl;

    if (gold and not anchored_lexicals->count(std::make_tuple(position,*r)))
    {
      v += 1.0;
    }
    return v;
}


double
nn_scorer::compute_internal_rule_score(const Production* r)
{

  double v = rules_expressions.at(r);
  return v;
}

double
nn_scorer::compute_internal_span_score(int begin,
                                       int end,
                                       int /*medium*/,
                                       int lhs)
{
  //  return 0.0;



  //  auto t = std::make_tuple(begin,end,medium,lhs);
  auto t = std::make_tuple(begin,end,
                           SymbolTable::instance_nt().translate(lhs)[0] == '[' ? 0 : 1
                           );

  double v = spans_expressions.at(t);

  return v;
}



void nn_scorer::set_words(const std::vector<Word>& w)
{
  words = &w;
}

double nn_scorer::compute_unary_score(int begin, int end, const MetaProduction* mp)
{

    auto r = static_cast<const Production*>(mp);

    //std::cerr << *r << std::endl;

    double v = compute_internal_rule_score(r);

    //std::cerr << "b4 un" << std::endl;

    //std::cerr << "un: " << begin << " " << words[begin]  << (end -1) << " " << words[end-1] << std::endl;

    if (end - begin > 2)
      v += compute_internal_span_score(begin, end -1, -1, r->get_lhs());

    //       std::cerr << "after un" << std::endl;


    if (gold and not anchored_unaries->count(std::make_tuple(begin,end,*r)))
      v += 1.0;
    // else
    // {
    //   std::cerr << "gold unary: " << *r << " " << begin << " " << end << std::endl;
    // }


    return v;

  }

double nn_scorer::compute_binary_score(int s, int e, int m, const MetaProduction* mp)
{
    auto r = static_cast<const Production*>(mp);

    double v = compute_internal_rule_score(r);

    if (e - s > 2) v+= compute_internal_span_score(s, e - 1, m, r->get_lhs());

    if (gold and not anchored_binaries->count(std::make_tuple(s,e,m,*r))) v += 1.0;

    return v;
  }

void nn_scorer::precompute_rule_expressions(const std::vector<Rule>& brules,
                                            const std::vector<Rule>& urules)
{
  for (const auto& r : brules)
  {
    auto e = rule_expression(r.get_lhs(), r.get_rhs0(), r.get_rhs1());
    auto v = as_scalar(cg->get_value(e.i));
    rules_expressions[&r] = v;
  }

  for (const auto& r : urules)
  {
    auto e =  rule_expression(r.get_lhs(),
                              r.get_rhs0(),
                              SymbolTable::instance_nt().get_symbol_count());
    double v = as_scalar(cg->get_value(e.i));
    rules_expressions[&r] = v;
  }
}


void nn_scorer::precompute_span_expressions(const std::vector<int>& /*lhs_set*/)
{
  dynet::expr::Expression Wleft = dynet::expr::parameter(*cg, _p_Wleft_span);
  dynet::expr::Expression Wright = dynet::expr::parameter(*cg, _p_Wright_span);

  std::vector<std::vector<dynet::expr::Expression>> lefts, rights;
  for (unsigned l = 0; l < 2; ++l)
  {
    lefts.push_back(std::vector<dynet::expr::Expression>());
    rights.push_back(std::vector<dynet::expr::Expression>());
    for (unsigned i = 0; i < words->size(); ++i)
    {
      auto e = dynet::expr::concatenate({lstm_embeddings[i],
                                         dynet::expr::input(*cg, l)});
      //if (i < words->size() - 2)
      lefts[l].push_back(Wleft * e);
      //if (i >= 2)
      rights[l].push_back(Wright * e);
    }
  }

  dynet::expr::Expression b = dynet::expr::parameter(*cg, _p_b_span);
  dynet::expr::Expression o = dynet::expr::parameter(*cg, _p_o_span);

  for (unsigned l = 0; l < 2; ++l)
    for (unsigned i = 0; i < words->size(); ++i)
      for (unsigned j = i; j < words->size(); ++j)
      {
        auto e = o * dynet::expr::rectify(lefts[l][i] + rights[l][j] + b);
        spans_expressions[std::make_tuple(i,j,l)] = as_scalar(cg->get_value(e.i));
      }
}

dynet::expr::Expression nn_scorer::span_expression(int lhs, int word_position_begin, int word_position_end)
{
  std::cerr << lstm_embeddings.size() << std::endl;
  // std::cerr << word_position_end << std::endl;

  //std::cerr << lhs << SymbolTable::instance_nt().translate(lhs) << std::endl;

  dynet::expr::Expression Wleft = dynet::expr::parameter(*cg, _p_Wleft_span);
  dynet::expr::Expression Wright = dynet::expr::parameter(*cg, _p_Wright_span);


  // seems that we cannot access ptbpstree namespace...
  auto art = SymbolTable::instance_nt().translate(lhs)[0] == '[' ? 0.0 : 1.0;


  // todo do we need this information on both ends ?
  dynet::expr::Expression i_left = Wleft * dynet::expr::concatenate({lstm_embeddings[word_position_begin],
                                                                     dynet::expr::input(*cg, art)});

  dynet::expr::Expression i_right = Wright * dynet::expr::concatenate({lstm_embeddings[word_position_end],
                                                                       dynet::expr::input(*cg, art)});

  dynet::expr::Expression b = dynet::expr::parameter(*cg, _p_b_span);
  dynet::expr::Expression o = dynet::expr::parameter(*cg, _p_o_span);

  return o * dynet::expr::rectify(i_left + i_right + b);
}

dynet::expr::Expression nn_scorer::rule_expression(int lhs, int rhs0, int rhs1)
{
  dynet::expr::Expression i = dynet::expr::concatenate({dynet::expr::lookup(*cg,_p_nts,rhs0),
                                                        dynet::expr::lookup(*cg,_p_nts,rhs1),
                                                        dynet::expr::lookup(*cg,_p_nts,lhs)});

  dynet::expr::Expression W = dynet::expr::parameter(*cg, _p_W_int);
  dynet::expr::Expression b = dynet::expr::parameter(*cg, _p_b_int);
  dynet::expr::Expression o = dynet::expr::parameter(*cg, _p_o_int);

  return o * dynet::expr::rectify(W*i + b);
}


dynet::expr::Expression nn_scorer::lexical_rule_expression(int lhs, int word_idx)
{
  //  static int pad = SymbolTable::instance_word().get_symbol_count();

  dynet::expr::Expression i = dynet::expr::concatenate({lstm_embeddings[word_idx],
                                                        dynet::expr::lookup(*cg, _p_nts, lhs)
                                                        // ,
                                                        // word_idx == 0 ? dynet::expr::lookup(*cg, _p_word, pad) : lstm_embeddings[word_idx-1],
                                                        // word_idx == lstm_embeddings.size() - 1 ? dynet::expr::lookup(*cg, _p_word, pad) : lstm_embeddings[word_idx+1]
    });

  dynet::expr::Expression W = dynet::expr::parameter(*cg, _p_W_lex);
  dynet::expr::Expression b = dynet::expr::parameter(*cg, _p_b_lex);
  dynet::expr::Expression o = dynet::expr::parameter(*cg, _p_o_lex);

  return o * dynet::expr::rectify(W*i + b);


  // linear classifier with 2 simple features (POS and WORD, POS)
  // auto v = std::vector<float>(SymbolTable::instance_nt().get_symbol_count() * SymbolTable::instance_word().get_symbol_count()
  //                             + SymbolTable::instance_nt().get_symbol_count(), 0.0);
  // v[lhs * SymbolTable::instance_word().get_symbol_count() + words[word_position].get_id()] = 1.0;
  // v[SymbolTable::instance_nt().get_symbol_count() * SymbolTable::instance_word().get_symbol_count() + lhs] = 1.0;

  // dynet::expr::Expression i = dynet::expr::input(*cg,{SymbolTable::instance_nt().get_symbol_count() * SymbolTable::instance_word().get_symbol_count()
  //                           + SymbolTable::instance_nt().get_symbol_count()}, v);

  // return W * i;


}

// adapted from caio
void nn_scorer::precompute_embeddings()
{
  lstm_embeddings.clear();

  std::vector<dynet::expr::Expression> embeddings;
  for (const auto& w : (*words))
  {
    embeddings.push_back(dynet::expr::lookup(*cg,_p_word, w.get_id()));
  }


  std::vector<dynet::expr::Expression> lstm_forward, lstm_backward;

  // Build forward LSTM
  l2r_builder.new_graph(*cg);
  l2r_builder.start_new_sequence();
  for (const Expression& input : embeddings)
    lstm_forward.push_back(l2r_builder.add_input(input));

  // Build backward LSTM
  r2l_builder.new_graph(*cg);
  r2l_builder.start_new_sequence();
  for (int i = embeddings.size() -1; i >= 0; --i)
  {
    lstm_backward.push_back(r2l_builder.add_input(embeddings[i]));
  }

  for (unsigned int i = 0 ; i < lstm_forward.size() ; ++i)
  {
    Expression e = dynet::expr::concatenate({lstm_backward[lstm_backward.size() - i - 1], lstm_forward[i]});
    lstm_embeddings.push_back(e);
  }



  // lstm_embeddings = embeddings;
}
