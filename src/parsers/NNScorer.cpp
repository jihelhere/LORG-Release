#include "NNScorer.h"

#include "utils/SymbolTable.h"
#include "rules/Production.h"

#include "utils/hash_impl.h"

#include <mutex>

#define WORD_EMBEDDING_SIZE 100
#define NT_EMBEDDING_SIZE 50
#define HIDDEN_SIZE 200
#define LSTM_HIDDEN_SIZE 200

namespace d = dynet;
namespace de = d::expr;

// for concurrent accesses to the cg and the model
std::mutex cg_mutex;

bool nn_scorer::model_initialized = false;
d::Parameter nn_scorer::_p_W_int;
d::Parameter nn_scorer::_p_b_int;
d::Parameter nn_scorer::_p_o_int;

d::Parameter nn_scorer::_p_W_lex;
d::Parameter nn_scorer::_p_b_lex;
d::Parameter nn_scorer::_p_o_lex;

d::Parameter nn_scorer::_p_Wleft_span;
d::Parameter nn_scorer::_p_Wright_span;
d::Parameter nn_scorer::_p_b_span;
d::Parameter nn_scorer::_p_o_span;

d::LookupParameter nn_scorer::_p_word;
d::LookupParameter nn_scorer::_p_nts;


std::vector<d::LSTMBuilder> nn_scorer::l2r_builders;
std::vector<d::LSTMBuilder> nn_scorer::r2l_builders;




inline bool
is_artificial(int lhs)
{
  return SymbolTable::instance_nt().translate(lhs)[0] == '[' ;
}


nn_scorer::nn_scorer(d::Model& m, unsigned lex_level, unsigned sp_level) :
    cg(nullptr),
    rule_scores(),
    span_scores(),
    lexical_level(lex_level),
    span_level(sp_level),
    words(nullptr)


{
  std::lock_guard<std::mutex> guard(cg_mutex);

  if (not model_initialized)
  {
    _p_word = m.add_lookup_parameters(SymbolTable::instance_word().get_symbol_count()+1,
                                      {WORD_EMBEDDING_SIZE});
    _p_nts = m.add_lookup_parameters(SymbolTable::instance_nt().get_symbol_count()+1,
                                     {NT_EMBEDDING_SIZE});



    unsigned lex_input_dim = 0;
    unsigned span_input_dim_base = 0;
    switch (lex_level) {
      case 0: {
        lex_input_dim = 3*WORD_EMBEDDING_SIZE;
        span_input_dim_base = WORD_EMBEDDING_SIZE;

        break;
      }
      default:
        lex_input_dim = 2*LSTM_HIDDEN_SIZE;
        span_input_dim_base = 2*LSTM_HIDDEN_SIZE;
        break;
    }



    _p_W_int = m.add_parameters({HIDDEN_SIZE, NT_EMBEDDING_SIZE*3});
    _p_b_int = m.add_parameters({HIDDEN_SIZE});
    _p_o_int = m.add_parameters({1,HIDDEN_SIZE});

    _p_W_lex = m.add_parameters({HIDDEN_SIZE, NT_EMBEDDING_SIZE + lex_input_dim});
    _p_b_lex = m.add_parameters({HIDDEN_SIZE});
    _p_o_lex = m.add_parameters({1,HIDDEN_SIZE});

    // linear classifier with 2 simple features (POS and WORD, POS)
    // _p_W_lex(m.add_parameters({1, SymbolTable::instance_nt().get_symbol_count() * SymbolTable::instance_word().get_symbol_count()
    //                               + SymbolTable::instance_nt().get_symbol_count()})),

    if (span_level > 0)
    {
      unsigned span_input_dim = 0;
      switch (span_level) {
        case 1: { // only the span
          span_input_dim = span_input_dim_base;
          break;
        }
        case 2: { // the span + natural/artificial nt
          span_input_dim = span_input_dim_base + 1;
          break;
        }
        default: { // the span + nt embeddding
          span_input_dim = span_input_dim_base + NT_EMBEDDING_SIZE;
          break;
        }
      }

      _p_Wleft_span = m.add_parameters({HIDDEN_SIZE, span_input_dim });
      _p_Wright_span = m.add_parameters({HIDDEN_SIZE, span_input_dim });
      _p_b_span = m.add_parameters({HIDDEN_SIZE});
      _p_o_span = m.add_parameters({1,HIDDEN_SIZE});
    }



    switch (lex_level) {
      case 0: {

        break;
      }
      case 1: {
        l2r_builders = std::vector<d::LSTMBuilder>(1,d::LSTMBuilder(2, WORD_EMBEDDING_SIZE, LSTM_HIDDEN_SIZE, &m));
        r2l_builders = std::vector<d::LSTMBuilder>(1,d::LSTMBuilder(2, WORD_EMBEDDING_SIZE, LSTM_HIDDEN_SIZE, &m));
        break;
      }
      default:
        l2r_builders = std::vector<d::LSTMBuilder>(1,d::LSTMBuilder(1, WORD_EMBEDDING_SIZE, LSTM_HIDDEN_SIZE, &m));
        r2l_builders = std::vector<d::LSTMBuilder>(1,d::LSTMBuilder(1, WORD_EMBEDDING_SIZE, LSTM_HIDDEN_SIZE, &m));

        for (unsigned s = 1; s < lex_level; ++s)
        {
          l2r_builders.push_back(d::LSTMBuilder(1, 2*LSTM_HIDDEN_SIZE, LSTM_HIDDEN_SIZE, &m));
          r2l_builders.push_back(d::LSTMBuilder(1, 2*LSTM_HIDDEN_SIZE, LSTM_HIDDEN_SIZE, &m));
        }
        break;
    }

    model_initialized = true;
  }
}


void nn_scorer::clear()
{
  rule_scores.clear();
  span_scores.clear();

  // other members are managed somewhere else!
}

void nn_scorer::set_gold(std::unordered_set<anchored_binrule_type>& ancbin,
                         std::unordered_set<anchored_unirule_type>& ancuni,
                         std::unordered_set<anchored_lexrule_type>& anclex)
{
  anchored_binaries = &ancbin;
  anchored_unaries  = &ancuni;
  anchored_lexicals = &anclex;

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

  // next function call is protected by non-recursive mutex
  auto e = lexical_rule_expression(mp->get_lhs(), position);

  std::lock_guard<std::mutex> guard(cg_mutex);
  double v = as_scalar(cg->get_value(e.i));


  if (gold && !anchored_lexicals->count(std::make_tuple(position,*r)))
  {
    v += 1.0;
  }
  return v;
}


double
nn_scorer::compute_internal_rule_score(const Production* r)
{

  double v = rule_scores.at(r);
  return v;
}

double
nn_scorer::compute_internal_span_score(int begin,
                                       int end,
                                       int /*medium*/,
                                       int lhs)
{

  //  auto t = std::make_tuple(begin,end,medium,lhs);
  double v;
  int lhs_info = 0;
  switch (span_level) {
    case 1: {
      lhs_info = 0;
      break;
    }
    case 2: {
      lhs_info = is_artificial(lhs) ? 0 : 1;
      break;
    }
    case 3: {
      lhs_info = lhs;
      break;
    }
    default:
      break;
  }

  v = span_scores.at(std::make_tuple(begin,end, lhs_info));
  return v;
}



void nn_scorer::set_words(const std::vector<Word>& w)
{
  words = &w;
}

double nn_scorer::compute_unary_score(int begin, int end, const MetaProduction* mp)
{

  auto r = static_cast<const Production*>(mp);

  double v = compute_internal_rule_score(r);

  if (gold && !anchored_unaries->count(std::make_tuple(begin,end,*r)))
    v += 1.0;

  return v;

}

double nn_scorer::compute_binary_score(int s, int e, int m, const MetaProduction* mp)
{
  auto r = static_cast<const Production*>(mp);

  double v = compute_internal_rule_score(r);

  if ((span_level > 0) and (e - s > 2)) v+= compute_internal_span_score(s, e - 1, m, r->get_lhs());

  if (gold and not anchored_binaries->count(std::make_tuple(s,e,m,*r))) v += 1.0;

  return v;
}

void nn_scorer::precompute_rule_expressions(const std::vector<Rule>& brules,
                                            const std::vector<Rule>& urules)
{
  int pad = SymbolTable::instance_nt().get_symbol_count();

  for (const auto& r : brules)
  {
    auto e = rule_expression(r.get_lhs(), r.get_rhs0(), r.get_rhs1());
    auto v = as_scalar(cg->get_value(e.i));
    rule_scores[&r] = v;
  }

  for (const auto& r : urules)
  {
    auto e =  rule_expression(r.get_lhs(), r.get_rhs0(), pad);
    double v = as_scalar(cg->get_value(e.i));
    rule_scores[&r] = v;
  }
}


void nn_scorer::precompute_span_expressions(const std::vector<int>& lhs_int)
{
  std::lock_guard<std::mutex> guard(cg_mutex);

  auto Wleft  = de::parameter(*cg, _p_Wleft_span);
  auto Wright = de::parameter(*cg, _p_Wright_span);
  auto b = de::parameter(*cg, _p_b_span);
  auto o = de::parameter(*cg, _p_o_span);



  switch (span_level) {
    case 1: {
      std::vector<de::Expression> lefts, rights;
      for (unsigned i = 0; i < words->size(); ++i)
      {
        auto e = lstm_embeddings[i];
        //if (i < words->size() - 2)
        lefts.push_back(Wleft * e);
        //if (i >= 2)
        rights.push_back(Wright * e);
      }

      for (unsigned i = 0; i < words->size(); ++i)
        for (unsigned j = i; j < words->size(); ++j)
        {
          auto e = o * de::rectify(lefts[i] + rights[j] + b);
          span_scores[std::make_tuple(i,j,0)] = as_scalar(cg->get_value(e.i));
        }

      break;
    }
    case 2: {
      std::vector<std::vector<de::Expression>> lefts, rights;

      for (unsigned l = 0; l < 2; ++l)
      {
        lefts.push_back(std::vector<de::Expression>());
        rights.push_back(std::vector<de::Expression>());
        for (unsigned i = 0; i < words->size(); ++i)
        {
          auto e = de::concatenate({lstm_embeddings[i],
                                    de::input(*cg, l)});

          //if (i < words->size() - 2)
          lefts[l].push_back(Wleft * e);
          //if (i >= 2)
          rights[l].push_back(Wright * e);
        }
      }

      for (unsigned l = 0; l < 2; ++l)
        for (unsigned i = 0; i < words->size(); ++i)
          for (unsigned j = i; j < words->size(); ++j)
          {
            auto e = o * de::rectify(lefts[l][i] + rights[l][j] + b);
            span_scores[std::make_tuple(i,j,l)] = as_scalar(cg->get_value(e.i));
          }


      break;
    }
    default:
      std::vector<std::vector<de::Expression>> lefts, rights;

      for (unsigned l = 0; l < lhs_int.size(); ++l)
      {
        lefts.push_back(std::vector<de::Expression>());
        rights.push_back(std::vector<de::Expression>());
        for (unsigned i = 0; i < words->size(); ++i)
        {
          auto e = de::concatenate({lstm_embeddings[i],
                                    de::lookup(*cg, _p_nts, lhs_int[l])});
          //if (i < words->size() - 2)
          lefts[l].push_back(Wleft * e);
          //if (i >= 2)
          rights[l].push_back(Wright * e);
        }
      }

      for (unsigned l =0;l < lhs_int.size(); ++l)
        for (unsigned i = 0; i < words->size(); ++i)
          for (unsigned j = i; j < words->size(); ++j)
          {
            auto e = o * de::rectify(lefts[l][i] + rights[l][j] + b);
            span_scores[std::make_tuple(i,j,lhs_int[l])] = as_scalar(cg->get_value(e.i));
          }
      break;
  }




}

de::Expression nn_scorer::span_expression(int lhs, int begin, int end)
{
  std::lock_guard<std::mutex> guard(cg_mutex);

  auto Wleft  = de::parameter(*cg, _p_Wleft_span);
  auto Wright = de::parameter(*cg, _p_Wright_span);
  auto b = de::parameter(*cg, _p_b_span);
  auto o = de::parameter(*cg, _p_o_span);


  switch (span_level) {
    case 1: {
      auto i_left = lstm_embeddings[begin];
      auto i_right = lstm_embeddings[end];
      return o * de::rectify(Wleft * i_left + Wright * i_right + b);
      break;
    }
    case 2:
      {
        auto type_symbol = is_artificial(lhs) ? 0.0 : 1.0;
        // todo do we need this information on both ends ?
        auto i_left = de::concatenate({lstm_embeddings[begin],
                                       de::input(*cg, type_symbol)});

        auto i_right = de::concatenate({lstm_embeddings[end],
                                        de::input(*cg, type_symbol)});
        return o * de::rectify(Wleft * i_left + Wright * i_right + b);
        break;
      }
    default:
        // todo do we need this information on both ends ?
        auto i_left = de::concatenate({lstm_embeddings[begin],
                                       de::lookup(*cg, _p_nts, lhs)});

        auto i_right = de::concatenate({lstm_embeddings[end],
                                        de::lookup(*cg, _p_nts, lhs)});
        return o * de::rectify(Wleft * i_left + Wright * i_right + b);
        break;
  }

}

de::Expression nn_scorer::rule_expression(int lhs, int rhs0, int rhs1)
{
  std::lock_guard<std::mutex> guard(cg_mutex);

  auto i = de::concatenate({de::lookup(*cg,_p_nts,rhs0),
                            de::lookup(*cg,_p_nts,rhs1),
                            de::lookup(*cg,_p_nts,lhs)});

  auto W = de::parameter(*cg, _p_W_int);
  auto b = de::parameter(*cg, _p_b_int);
  auto o = de::parameter(*cg, _p_o_int);
  return o * de::rectify(W*i + b);
}


de::Expression nn_scorer::lexical_rule_expression(int lhs, unsigned word_idx)
{
  static int pad = SymbolTable::instance_word().get_symbol_count();

  std::lock_guard<std::mutex> guard(cg_mutex);

  auto i = lexical_level > 0 ? de::concatenate({lstm_embeddings[word_idx],
                                                de::lookup(*cg, _p_nts, lhs)})
           :
           de::concatenate({lstm_embeddings[word_idx],
                            de::lookup(*cg, _p_nts, lhs),
                            word_idx == 0 ? de::lookup(*cg, _p_word, pad) : lstm_embeddings[word_idx-1],
                            word_idx == lstm_embeddings.size() - 1 ? de::lookup(*cg, _p_word, pad) : lstm_embeddings[word_idx+1]
             })
           ;

  auto W = de::parameter(*cg, _p_W_lex);
  auto b = de::parameter(*cg, _p_b_lex);
  auto o = de::parameter(*cg, _p_o_lex);

  return o * de::rectify(W*i + b);


  // linear classifier with 2 simple features (POS and WORD, POS)
  // auto v = std::vector<float>(SymbolTable::instance_nt().get_symbol_count() * SymbolTable::instance_word().get_symbol_count()
  //                             + SymbolTable::instance_nt().get_symbol_count(), 0.0);
  // v[lhs * SymbolTable::instance_word().get_symbol_count() + words[word_position].get_id()] = 1.0;
  // v[SymbolTable::instance_nt().get_symbol_count() * SymbolTable::instance_word().get_symbol_count() + lhs] = 1.0;

  // de::Expression i = de::input(*cg,{SymbolTable::instance_nt().get_symbol_count() * SymbolTable::instance_word().get_symbol_count()
  //                           + SymbolTable::instance_nt().get_symbol_count()}, v);

  // return W * i;


}

// adapted from caio
void nn_scorer::precompute_embeddings()
{
  lstm_embeddings.clear();

  std::lock_guard<std::mutex> guard(cg_mutex);

  for (const auto& w : (*words))
  {
    lstm_embeddings.push_back(de::lookup(*cg,_p_word, w.get_id()));
  }

  // finish here without lstms
  if (lexical_level ==  0)
    return;
  else
  {
    for (unsigned s = 0; s < lexical_level; ++s)
    {
      std::vector<de::Expression> lstm_forward, lstm_backward;

      // Build forward LSTM
      l2r_builders[s].new_graph(*cg);
      l2r_builders[s].start_new_sequence();
      for (const auto& input : lstm_embeddings)
      {
        lstm_forward.push_back(l2r_builders[s].add_input(input));
      }

      // Build backward LSTM
      r2l_builders[s].new_graph(*cg);
      r2l_builders[s].start_new_sequence();
      for (int i = lstm_embeddings.size() -1; i >= 0; --i)
      {
        lstm_backward.push_back(r2l_builders[s].add_input(lstm_embeddings[i]));
      }

      for (unsigned int i = 0 ; i < lstm_forward.size() ; ++i)
      {
        auto e = de::concatenate({lstm_backward[lstm_backward.size() - i - 1],
                                  lstm_forward[i]});

        lstm_embeddings[i] = e;
      }
    }
  }
}
