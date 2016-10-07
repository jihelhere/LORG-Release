#include "NNScorer.h"

#include "utils/SymbolTable.h"
#include "rules/Production.h"

#define WORD_EMBEDDING_SIZE 100
#define NT_EMBEDDING_SIZE 50
#define HIDDEN_SIZE 100

#define LSTM_HIDDEN_SIZE 150

nn_scorer::nn_scorer(cnn::Model& m) :
    cg(nullptr),

    _p_W_int(m.add_parameters({HIDDEN_SIZE, NT_EMBEDDING_SIZE*3})),
    _p_b_int(m.add_parameters({HIDDEN_SIZE})),
    _p_o_int(m.add_parameters({1,HIDDEN_SIZE})),

    _p_W_lex(m.add_parameters({HIDDEN_SIZE, NT_EMBEDDING_SIZE+3*WORD_EMBEDDING_SIZE})),
    _p_b_lex(m.add_parameters({HIDDEN_SIZE})),
    _p_o_lex(m.add_parameters({1,HIDDEN_SIZE})),

    // linear classifier with 2 simple features (POS and WORD, POS)
    // _p_W_lex(m.add_parameters({1, SymbolTable::instance_nt().get_symbol_count() * SymbolTable::instance_word().get_symbol_count()
    //                               + SymbolTable::instance_nt().get_symbol_count()})),

    // _p_W_span(m.add_parameters({HIDDEN_SIZE, 2*LSTM_HIDDEN_SIZE+NT_EMBEDDING_SIZE})),
    // _p_b_span(m.add_parameters({HIDDEN_SIZE})),
    // _p_o_span(m.add_parameters({1,HIDDEN_SIZE})),


    _p_word(m.add_lookup_parameters(SymbolTable::instance_word().get_symbol_count()+1,
                                    {WORD_EMBEDDING_SIZE})),
    _p_nts(m.add_lookup_parameters(SymbolTable::instance_nt().get_symbol_count()+1,
                                   {NT_EMBEDDING_SIZE})),


    // l2r_builder(2, WORD_EMBEDDING_SIZE, LSTM_HIDDEN_SIZE, &m),
    // r2l_builder(2, WORD_EMBEDDING_SIZE, LSTM_HIDDEN_SIZE, &m),



    rules_expressions(),
    spans_expressions()

{}


void nn_scorer::clear()
{
  rules_expressions.clear();
  spans_expressions.clear();

  // other members are managed somewhere else!
}

void nn_scorer::set_gold(std::vector<anchored_binrule_type>& ancbin,
                         std::vector<anchored_unirule_type>& ancuni,
                         std::vector<anchored_lexrule_type>& anclex)
{
  anchored_binaries.clear();
  anchored_binaries.insert(ancbin.begin(), ancbin.end());

  anchored_unaries.clear();
  anchored_unaries.insert(ancuni.begin(), ancuni.end());

  anchored_lexicals.clear();
  anchored_lexicals.insert(anclex.begin(), anclex.end());
  gold = true;
}

void nn_scorer::unset_gold()
{
  anchored_binaries.clear();
  anchored_unaries.clear();
  anchored_lexicals.clear();
  gold = false;
}



double nn_scorer::compute_lexical_score(int position, const MetaProduction* mp)
{
    auto r = static_cast<const Production*>(mp);

    //std::cerr << *r << std::endl;

    //double v = rules_expressions.at(r).second;

    double v = as_scalar(cg->get_value(lexical_rule_expression(mp->get_lhs(), position).i));

    //std::cerr << "after" << std::endl;

    if (gold and not anchored_lexicals.count(std::make_tuple(position,*r)))
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
                                       int medium,
                                       int lhs)
{
  //  return 0.0;



  //  auto t = std::make_tuple(begin,end,medium,lhs);
  auto t = std::make_tuple(begin,lhs);

  double v = as_scalar(cg->get_value(spans_expressions[t].i));

  return v;
}



void nn_scorer::set_words(const std::vector<Word>& w)
{
  words = w;
}

double nn_scorer::compute_unary_score(int begin, int end, const MetaProduction* mp)
{

    auto r = static_cast<const Production*>(mp);

    //std::cerr << *r << std::endl;

    double v = compute_internal_rule_score(r);

    //std::cerr << "b4 un" << std::endl;



    //std::cerr << "un: " << begin << " " << words[begin]  << (end -1) << " " << words[end-1] << std::endl;

    // v+= compute_internal_span_score(begin,
    //                                 end - 1,
    //                                 -1,
    //                                 r->get_lhs(),
    //                                 expv);

    //       std::cerr << "after un" << std::endl;


    if (gold and not anchored_unaries.count(std::make_tuple(begin,end,*r)))
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

    //    std::cerr << "b4 bin" << std::endl;

    //std::cerr << "bin: " << s << " " << words[s]  << (e-1) << " " << words[e-1]  << m << " " << words[m] << std::endl;


    // v+= compute_internal_span_score(s,
    //                                 e - 1,
    //                                 m,
    //                                 r->get_lhs(),
    //                                 expp);
    //          std::cerr << "after bin" << std::endl;


    if (gold and not anchored_binaries.count(std::make_tuple(s,e,m,*r)))
      v += 1.0;
    // else
    // {
    //   std::cerr << "gold binary: " << *r << " " << s << " " << e << " " << m << std::endl;
    // }
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

  cg->forward();
}


void nn_scorer::precompute_span_expressions(const std::unordered_set<int>& lhs_set)
{
  for (auto lhs : lhs_set)
  for (unsigned i = 0; i < words.size(); ++i)
  {
    spans_expressions[std::make_tuple(i,lhs)] = span_expression(lhs, i);
  }
  cg->incremental_forward();
}

cnn::expr::Expression nn_scorer::span_expression(int lhs, int word_position)
{
   cnn::expr::Expression i = cnn::expr::concatenate({lstm_embeddings[word_position],
                                                     // cnn::expr::lookup(*cg,_p_word,end_id),
                                                     // cnn::expr::lookup(*cg,_p_word,medium_id),
                                                     cnn::expr::lookup(*cg,_p_nts,lhs)});

      cnn::expr::Expression W = cnn::expr::parameter(*cg, _p_W_span);
      cnn::expr::Expression b = cnn::expr::parameter(*cg, _p_b_span);
      cnn::expr::Expression o = cnn::expr::parameter(*cg, _p_o_span);

      return o * cnn::expr::rectify(W*i + b);
}


cnn::expr::Expression nn_scorer::rule_expression(int lhs, int rhs0, int rhs1)
{
  cnn::expr::Expression i = cnn::expr::concatenate({cnn::expr::lookup(*cg,_p_nts,rhs0),
                                                    cnn::expr::lookup(*cg,_p_nts,rhs1),
                                                    cnn::expr::lookup(*cg,_p_nts,lhs)});

  cnn::expr::Expression W = cnn::expr::parameter(*cg, _p_W_int);
  cnn::expr::Expression b = cnn::expr::parameter(*cg, _p_b_int);
  cnn::expr::Expression o = cnn::expr::parameter(*cg, _p_o_int);

  return o * cnn::expr::rectify(W*i + b);
}


cnn::expr::Expression nn_scorer::lexical_rule_expression(int lhs, int word_position)
{
  cnn::expr::Expression i = cnn::expr::concatenate({lstm_embeddings[word_position],
                                                    cnn::expr::lookup(*cg, _p_nts, lhs)
                                                    ,
                                                    word_position == 0 ? cnn::expr::lookup(*cg, _p_word, SymbolTable::instance_word().get_symbol_count()) : lstm_embeddings[word_position-1],
                                                    word_position == lstm_embeddings.size() - 1 ? cnn::expr::lookup(*cg, _p_word, SymbolTable::instance_word().get_symbol_count()) : lstm_embeddings[word_position+1]
    });

  cnn::expr::Expression W = cnn::expr::parameter(*cg, _p_W_lex);
  cnn::expr::Expression b = cnn::expr::parameter(*cg, _p_b_lex);
  cnn::expr::Expression o = cnn::expr::parameter(*cg, _p_o_lex);

  return o * cnn::expr::rectify(W*i + b);


  // linear classifier with 2 simple features (POS and WORD, POS)
  // auto v = std::vector<float>(SymbolTable::instance_nt().get_symbol_count() * SymbolTable::instance_word().get_symbol_count()
  //                             + SymbolTable::instance_nt().get_symbol_count(), 0.0);
  // v[lhs * SymbolTable::instance_word().get_symbol_count() + words[word_position].get_id()] = 1.0;
  // v[SymbolTable::instance_nt().get_symbol_count() * SymbolTable::instance_word().get_symbol_count() + lhs] = 1.0;

  // cnn::expr::Expression i = cnn::expr::input(*cg,{SymbolTable::instance_nt().get_symbol_count() * SymbolTable::instance_word().get_symbol_count()
  //                           + SymbolTable::instance_nt().get_symbol_count()}, v);

  // return W * i;


}

// adapted from caio
void nn_scorer::precompute_embeddings()
{
  lstm_embeddings.clear();


  std::vector<cnn::expr::Expression> embeddings;

  // base embeddings (could be avoided)
  for (const auto& w : words)
  {
    //std::cerr << "set 1 embedding" << std::endl;
    embeddings.push_back(cnn::expr::lookup(*cg,_p_word, w.get_id()));
  }


  // std::vector<cnn::expr::Expression> lstm_forward, lstm_backward;

  // // Build forward LSTM
  // l2r_builder.new_graph(*cg);
  // l2r_builder.start_new_sequence();
  // for (const Expression& input : embeddings)
  //   lstm_forward.push_back(l2r_builder.add_input(input));

  // // Build backward LSTM
  // r2l_builder.new_graph(*cg);
  // r2l_builder.start_new_sequence();
  // for (int i = embeddings.size() -1; i >= 0; --i)
  // {
  //   lstm_backward.push_back(r2l_builder.add_input(embeddings[i]));
  // }

  // for (unsigned int i = 0 ; i < lstm_forward.size() ; ++i)
  // {
  //   Expression e = cnn::expr::concatenate({lstm_backward[lstm_backward.size() - i - 1], lstm_forward[i]});
  //   lstm_embeddings.push_back(e);
  // }


 lstm_embeddings = embeddings;
}
