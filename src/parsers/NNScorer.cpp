#include "NNScorer.h"

#include "utils/SymbolTable.h"
#include "rules/Production.h"

#include "utils/hash_impl.h"

#include <mutex>

namespace d = dynet;
namespace de = d::expr;

// for concurrent accesses to the cg and the model
std::mutex cg_mutex;

bool nn_scorer::model_initialized = false;
bool nn_scorer::train_mode = false;


d::ComputationGraph * nn_scorer::cg = nullptr;

d::Parameter nn_scorer::_p_W_int;
d::Parameter nn_scorer::_p_b_int;
d::Parameter nn_scorer::_p_o_int;

d::Parameter nn_scorer::_p_W_lex;
d::Parameter nn_scorer::_p_b_lex;
d::Parameter nn_scorer::_p_o_lex;

d::Parameter nn_scorer::_p_W_span_init;
// d::Parameter nn_scorer::_p_W_span_end;
// d::Parameter nn_scorer::_p_W_span_split;

d::Parameter nn_scorer::_p_b_span_init;
d::Parameter nn_scorer::_p_o_span_init;

d::Parameter nn_scorer::_p_W_span_left;
d::Parameter nn_scorer::_p_W_span_right;
d::Parameter nn_scorer::_p_W_span_mid;

d::Parameter nn_scorer::_p_W_span_distance;
d::Parameter nn_scorer::_p_W_span_extra;
d::Parameter nn_scorer::_p_b_span_bin;
d::Parameter nn_scorer::_p_o_span_bin;
d::Parameter nn_scorer::_p_b_span_un;
d::Parameter nn_scorer::_p_o_span_un;

d::LookupParameter nn_scorer::_p_word;
d::LookupParameter nn_scorer::_p_nts;


std::vector<d::LSTMBuilder> nn_scorer::word_l2r_builders;
std::vector<d::LSTMBuilder> nn_scorer::word_r2l_builders;

d::LSTMBuilder nn_scorer::letter_l2r_builder;
d::LSTMBuilder nn_scorer::letter_r2l_builder;

std::vector<de::Expression> nn_scorer::rule_expressions;
std::vector<double> nn_scorer::rule_scores;

inline bool
is_artificial(int lhs)
{
  return SymbolTable::instance_nt().translate(lhs)[0] == '[' ;
}


nn_scorer::nn_scorer(d::Model& m, unsigned lex_level, unsigned sp_level,
                     unsigned word_embedding_size,
                     unsigned nt_embedding_size,
                     unsigned hidden_size,
                     unsigned lstm_hidden_size,
                     bool char_emb,
                     bool span_midp
                     ) :
    span_scores_bin(),
    span_scores_un(),
    lexical_level(lex_level),
    span_level(sp_level),
    words(nullptr),
    use_char_emb(char_emb),
    use_span_midpoints(span_midp)
{
  std::lock_guard<std::mutex> guard(cg_mutex);

  if (not model_initialized)
  {
    rule_scores.resize(SymbolTable::instance_nt().get_symbol_count() *
                       SymbolTable::instance_nt().get_symbol_count() *
                       (SymbolTable::instance_nt().get_symbol_count() + 1) );

    rule_expressions.resize(SymbolTable::instance_nt().get_symbol_count() *
                            SymbolTable::instance_nt().get_symbol_count() *
                            (SymbolTable::instance_nt().get_symbol_count() +1) );

    if (use_char_emb)
    {
      // todo manage unicode ??? use a dict ???
      letter_l2r_builder = d::LSTMBuilder(1, 1, word_embedding_size / 2, &m);
      letter_r2l_builder = d::LSTMBuilder(1, 1, word_embedding_size - (word_embedding_size / 2), &m);
    }

    // needed for padding when not using word lstms
    // todo : use one vector for padding
    _p_word = m.add_lookup_parameters(SymbolTable::instance_word().get_symbol_count()+1,
                                      {word_embedding_size});

    _p_nts = m.add_lookup_parameters(SymbolTable::instance_nt().get_symbol_count()+1,
                                     {nt_embedding_size});



    unsigned lex_input_dim = 0;
    unsigned span_input_dim_base = 0;
    switch (lex_level) {
      case 0: {
        lex_input_dim = 3*word_embedding_size;
        span_input_dim_base = word_embedding_size;

        break;
      }
      default:
        lex_input_dim = lstm_hidden_size;
        span_input_dim_base = lstm_hidden_size;
        break;
    }

    // grammar rules
    _p_W_int = m.add_parameters({hidden_size, nt_embedding_size*3});
    _p_b_int = m.add_parameters({hidden_size});
    _p_o_int = m.add_parameters({1,hidden_size});

    // lexical rules
    _p_W_lex = m.add_parameters({hidden_size, nt_embedding_size + lex_input_dim});
    _p_b_lex = m.add_parameters({hidden_size});
    _p_o_lex = m.add_parameters({1,hidden_size});


    if (span_level > 0)
    {
      _p_W_span_init = m.add_parameters({hidden_size, span_input_dim_base + nt_embedding_size});
      _p_b_span_init = m.add_parameters({hidden_size});
      _p_o_span_init = m.add_parameters({1,hidden_size});
      // _p_W_span_end =  m.add_parameters({hidden_size, span_input_dim_base + nt_embedding_size});
      // _p_W_span_mid =  m.add_parameters({hidden_size, span_input_dim_base + nt_embedding_size})


      _p_W_span_left = m.add_parameters({hidden_size, span_input_dim_base });
      _p_W_span_right = m.add_parameters({hidden_size, span_input_dim_base });
      //if (use_span_midpoints) // TODO templatize code to be able to remove this
      _p_W_span_mid = m.add_parameters({hidden_size, span_input_dim_base });
      _p_W_span_distance = m.add_parameters({hidden_size,1});

      switch (span_level) {
        case 1: { // only the span with distance
          break;
        }
        case 2: { // the span + distance + natural/artificial nt
          _p_W_span_extra = m.add_parameters({hidden_size,1});
          break;
        }
        default: { // the span + distance + nt embeddding
          _p_W_span_extra = m.add_parameters({hidden_size,nt_embedding_size});
          break;
        }
      }
      _p_b_span_bin = m.add_parameters({hidden_size});
      _p_o_span_bin = m.add_parameters({1,hidden_size});

      _p_b_span_un = m.add_parameters({hidden_size});
      _p_o_span_un = m.add_parameters({1,hidden_size});
    }

    switch (lex_level) {
      case 0: {

        break;
      }
      case 1: {
        word_l2r_builders = std::vector<d::LSTMBuilder>(1,d::LSTMBuilder(2, word_embedding_size, lstm_hidden_size, &m));
        word_r2l_builders = std::vector<d::LSTMBuilder>(1,d::LSTMBuilder(2, word_embedding_size, lstm_hidden_size, &m));
        break;
      }
      default:
        word_l2r_builders = std::vector<d::LSTMBuilder>(1,d::LSTMBuilder(1, word_embedding_size, lstm_hidden_size, &m));
        word_r2l_builders = std::vector<d::LSTMBuilder>(1,d::LSTMBuilder(1, word_embedding_size, lstm_hidden_size, &m));

        for (unsigned s = 1; s < lex_level; ++s)
        {
          word_l2r_builders.push_back(d::LSTMBuilder(1, lstm_hidden_size, lstm_hidden_size, &m));
          word_r2l_builders.push_back(d::LSTMBuilder(1, lstm_hidden_size, lstm_hidden_size, &m));
        }
        break;
    }

    model_initialized = true;
  }
}


void nn_scorer::clear()
{
  span_scores_bin.clear();
  span_scores_un.clear();

  lexical_expressions.clear();
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
  double v;

  auto&& r = static_cast<const Production*>(mp);

  // next function call is protected by non-recursive mutex
  auto&& e = lexical_rule_expression(mp->get_lhs(), position);

  {
    std::lock_guard<std::mutex> guard(cg_mutex);
    v = as_scalar(cg->get_value(e.i));
  }

  if (gold && !anchored_lexicals->count(std::make_tuple(position,*r)))
  {
    v += 1.0;
  }
  return v;
}


double
nn_scorer::compute_internal_rule_score(const Production* r)
{
  static int pad = SymbolTable::instance_nt().get_symbol_count();
  return rule_scores[nt_triple_to_index(r->get_lhs(), r->get_rhs0(), (r->get_rhs().size() > 1 ? r->get_rhs1() : pad))];
}


// same code as span_expression ??
double
nn_scorer::compute_internal_span_score(int begin,
                                       int end,
                                       int medium,
                                       int lhs)
{
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


  double v = span_scores_init[begin][lhs];


  if (medium >=0 && not use_span_midpoints) medium = 0;

  return v + (medium >= 0 ?
              span_scores_bin.at(std::make_tuple(begin,end, medium, lhs_info))
              :
              span_scores_un.at(std::make_tuple(begin,end, lhs_info)));
}



void nn_scorer::set_words(const std::vector<Word>& w)
{
  words = &w;
}

double nn_scorer::compute_unary_score(int begin, int end, const MetaProduction* mp)
{

  auto r = static_cast<const Production*>(mp);

  double v = compute_internal_rule_score(r);

  if (span_level > 0) v+= compute_internal_span_score(begin, end - 1, -1, r->get_lhs());

  if (gold && !anchored_unaries->count(std::make_tuple(begin,end,*r)))
    v += 1.0;

  return v;

}

double nn_scorer::compute_binary_score(int s, int e, int m, const MetaProduction* mp)
{
  auto r = static_cast<const Production*>(mp);

  double v = compute_internal_rule_score(r);

  if (span_level > 0) v+= compute_internal_span_score(s, e - 1, m, r->get_lhs());

  if (gold and not anchored_binaries->count(std::make_tuple(s,e,m,*r))) v += 1.0;

  return v;
}

void nn_scorer::precompute_rule_expressions(const std::vector<Rule>& brules,
                                            const std::vector<Rule>& urules)
{
  int pad = SymbolTable::instance_nt().get_symbol_count();

  for (const auto& r : brules)
  {
    auto&& e = rule_expression(r.get_lhs(), r.get_rhs0(), r.get_rhs1());
    auto&& t = nt_triple_to_index(r.get_lhs(),r.get_rhs0(),r.get_rhs1());
    if (train_mode) rule_expressions[t] = e;
    rule_scores[t] = as_scalar(cg->get_value(e.i));
  }

  for (const auto& r : urules)
  {
    auto&& e =  rule_expression(r.get_lhs(), r.get_rhs0(), pad);
    auto&& t = nt_triple_to_index(r.get_lhs(),r.get_rhs0(),-1);
    if (train_mode) rule_expressions[t] = e;
    rule_scores[t] = as_scalar(cg->get_value(e.i));
  }
}


void nn_scorer::precompute_span_expressions(const std::vector<int>& lhs_int)
{
  std::lock_guard<std::mutex> guard(cg_mutex);

  auto&& Wl = de::parameter(*cg, _p_W_span_left);
  auto&& Wr = de::parameter(*cg, _p_W_span_right);
  auto&& Wm = de::parameter(*cg, _p_W_span_mid);

  auto&& Wd = de::parameter(*cg, _p_W_span_distance);
  auto&& bbin = de::parameter(*cg, _p_b_span_bin);
  auto&& obin = de::parameter(*cg, _p_o_span_bin);
  auto&& bun = de::parameter(*cg, _p_b_span_un);
  auto&& oun = de::parameter(*cg, _p_o_span_un);

  auto&& Wi = de::parameter(*cg,_p_W_span_init);
  auto&& bi = de::parameter(*cg,_p_b_span_init);
  auto&& oi = de::parameter(*cg,_p_o_span_init);

  std::vector<de::Expression> lefts,rights,mids,distances,extras;

  span_expressions_init.resize(words->size());
  span_scores_init.resize(words->size());

  for (unsigned i = 0; i < words->size(); ++i)
  {
    lefts.push_back(Wl * embeddings[i]);
    rights.push_back(Wr * embeddings[i]);
    if (use_span_midpoints) mids.push_back(Wm * embeddings[i]);
    distances.push_back(Wd * de::input(*cg, i));

    span_expressions_init[i].resize(SymbolTable::instance_nt().get_symbol_count());
    span_scores_init[i].resize(SymbolTable::instance_nt().get_symbol_count());
    for (auto l : lhs_int)
    {
      auto&& inp = de::concatenate({de::lookup(*cg,_p_nts,l), embeddings[i]});
      auto&& e = oi * de::rectify(Wi * inp + bi);
      if (train_mode) span_expressions_init[i][l] = e;
      span_scores_init[i][l] = as_scalar(cg->get_value(e.i));
    }
  }

  switch (span_level) {
    case 1: {
      for (unsigned i = 0; i < words->size(); ++i)
      {
        for (unsigned j = i; j < words->size(); ++j)
        {
          auto&& e = lefts[i] + rights[j] + distances[j-i];
          auto&& g = oun * de::rectify(e + bun);
          auto&& t = std::make_tuple(i,j,0);
          if (train_mode) span_expressions_un[t] = g;
          span_scores_un[t] = as_scalar(cg->get_value(g.i));


          if (use_span_midpoints)
            for (unsigned k = i+1; k <= j; ++k)
            {
              auto&& eb = e + mids[k];
              auto&& f = obin * de::rectify(eb + bbin);
              auto&& t = std::make_tuple(i,j,k,0);
              if (train_mode) span_expressions_bin[t] = f;
              span_scores_bin[t] = as_scalar(cg->get_value(f.i));
            }
          else
          {
            auto&& f = obin * de::rectify(e + bbin);
            auto&& t = std::make_tuple(i,j,0,0);
            if (train_mode) span_expressions_bin[t] = f;
            span_scores_bin[t] = as_scalar(cg->get_value(f.i));
          }
        }
      }
      break;
    }
    case 2: {
      auto&& We = de::parameter(*cg, _p_W_span_extra);
      std::vector<de::Expression> extras;
      for (unsigned l = 0; l < 2; ++l)
        extras.push_back(We * de::input(*cg, l));

      for (unsigned l = 0; l < 2; ++l)
      {
        auto&& e1 = extras[l];
        for (unsigned i = 0; i < words->size(); ++i)
        {
          auto&& e2 = e1 + lefts[i];
          for (unsigned j = i; j < words->size(); ++j)
          {
            auto&& e = e2 + rights[j] + distances[j-i];
            auto&& g = oun * de::rectify(e + bun);
            if (l==1) // artificial for unaries makes no sense
            {
              auto&& t = std::make_tuple(i,j,l);
              if (train_mode) span_expressions_un[t] = g;
              span_scores_un[t] = as_scalar(cg->get_value(g.i));
            }
            if (use_span_midpoints)
              for (unsigned k = i+1; k <= j; ++k)
              {
                auto&& eb = e + mids[k];
                auto&& f = obin * de::rectify(eb + bbin);
                auto&& t = std::make_tuple(i,j,k,l);
                if (train_mode) span_expressions_bin[t] = f;
                span_scores_bin[t] = as_scalar(cg->get_value(f.i));
              }
            else
            {
              auto&& f = obin * de::rectify(e + bbin);
              auto&& t = std::make_tuple(i,j,0,l);
              if (train_mode) span_expressions_bin[t] = f;
              span_scores_bin[t] = as_scalar(cg->get_value(f.i));
            }
          }
        }
      }
      break;
    }
    default:
      auto&& We = de::parameter(*cg, _p_W_span_extra);
      std::vector<de::Expression> extras;
      for (unsigned l =0;l < lhs_int.size(); ++l)
        extras.push_back(We *de::lookup(*cg, _p_nts, lhs_int[l]));

      for (unsigned l = 0;l < lhs_int.size(); ++l)
      {
        auto&& e1 = extras[l];
        for (unsigned i = 0; i < words->size(); ++i)
        {
          auto&& e2 = e1 + lefts[i];
          for (unsigned j = i; j < words->size(); ++j)
          {
            auto&& e = e2 + rights[j] + distances[j-i];
            auto&& g = oun * de::rectify(e + bun);

            auto&& t = std::make_tuple(i,j,lhs_int[l]);
            if (train_mode) span_expressions_un[t] = g;
            span_scores_un[t] = as_scalar(cg->get_value(g.i));

            if (use_span_midpoints)
              for (unsigned k = i+1; k <= j; ++k)
              {
                auto&& eb = e + mids[k];
                auto&& f = obin * de::rectify(eb + bbin);

                auto&& t = std::make_tuple(i,j,k,lhs_int[l]);
                if (train_mode) span_expressions_bin[t] = f;
                span_scores_bin[t] = as_scalar(cg->get_value(f.i));

              }
            else
            {
                auto&& f = obin * de::rectify(e + bbin);

                auto&& t = std::make_tuple(i,j,0,lhs_int[l]);
                if (train_mode) span_expressions_bin[t] = f;
                span_scores_bin[t] = as_scalar(cg->get_value(f.i));
            }
          }
        }
      }
      break;
  }
}


de::Expression& nn_scorer::span_init(int lhs, int begin)
{
  return span_expressions_init[begin][lhs];
}

de::Expression& nn_scorer::span_expression(int lhs, int begin, int end, int medium)
{
  int lhs_code = 0;
  switch (span_level) {
    case 1: {
      lhs_code = 0;
      break;
    }
    case 2:
      {
        lhs_code = is_artificial(lhs) ? 0 : 1;
        break;
      }
    default:
      lhs_code = lhs;
      break;
  }

  if (medium >=0 && not use_span_midpoints) medium = 0;

  return
      medium >= 0 ?
      span_expressions_bin.at(std::make_tuple(begin,end, medium, lhs_code))
      :
      span_expressions_un.at(std::make_tuple(begin,end, lhs_code));
}

de::Expression nn_scorer::rule_expression(int lhs, int rhs0, int rhs1)
{
  std::lock_guard<std::mutex> guard(cg_mutex);

  auto&& i = de::concatenate({de::lookup(*cg,_p_nts,rhs0),
                              de::lookup(*cg,_p_nts,rhs1),
                              de::lookup(*cg,_p_nts,lhs)});

  auto&& W = de::parameter(*cg, _p_W_int);
  auto&& b = de::parameter(*cg, _p_b_int);
  auto&& o = de::parameter(*cg, _p_o_int);
  return o * de::rectify(W*i + b);
}


de::Expression nn_scorer::lexical_rule_expression(int lhs, unsigned word_idx)
{
  static int pad = SymbolTable::instance_word().get_symbol_count();

  auto&& t = std::make_tuple(lhs,word_idx);

  if (lexical_expressions.count(t))
    return lexical_expressions[t];
  else
  {
    std::lock_guard<std::mutex> guard(cg_mutex);

    auto&& i = lexical_level > 0 ? de::concatenate({embeddings[word_idx],
                                                    de::lookup(*cg, _p_nts, lhs)})
               :
               de::concatenate({embeddings[word_idx],
                                de::lookup(*cg, _p_nts, lhs),
                                word_idx == 0 ? de::lookup(*cg, _p_word, pad) : embeddings[word_idx-1],
                                word_idx == embeddings.size() - 1 ? de::lookup(*cg, _p_word, pad) : embeddings[word_idx+1]
                 })
               ;

    auto&& W = de::parameter(*cg, _p_W_lex);
    auto&& b = de::parameter(*cg, _p_b_lex);
    auto&& o = de::parameter(*cg, _p_o_lex);

    return lexical_expressions[t] = o * de::rectify(W*i + b);
  }
}

// adapted from caio
void nn_scorer::precompute_embeddings()
{
  embeddings.clear();

  std::lock_guard<std::mutex> guard(cg_mutex);

  if (use_char_emb)
  {
    for (const auto& w : (*words))
    {
      // Build forward LSTM
      letter_l2r_builder.new_graph(*cg);
      letter_l2r_builder.start_new_sequence();
      // initialize with <w>
      letter_l2r_builder.add_input(de::input(*cg, 0.0));
      auto&& form = w.get_form();
      for (unsigned c = 0; c < form.size(); ++c)
      {
        letter_l2r_builder.add_input(de::input(*cg, float(form[c])));
      }

      // Build backward LSTM
      letter_r2l_builder.new_graph(*cg);
      letter_r2l_builder.start_new_sequence();
      // initialize with <w>
      letter_r2l_builder.add_input(de::input(*cg,0.0));
      for (int i = form.size() - 1; i >= 0; --i)
      {
        letter_r2l_builder.add_input(de::input(*cg,float(form[i])));
      }

      // finalize with </w> and concatenate
      embeddings.push_back(letter_l2r_builder.add_input(de::input(*cg,1.0)) +
                           letter_r2l_builder.add_input(de::input(*cg,1.0))
                           );
    }
  }

  else
  {
    for (const auto& w : (*words))
    {
      embeddings.push_back(de::lookup(*cg,_p_word, w.get_id()));
    }
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
      word_l2r_builders[s].new_graph(*cg);
      word_l2r_builders[s].start_new_sequence();
      for (const auto& input : embeddings)
      {
        lstm_forward.push_back(word_l2r_builders[s].add_input(input));
      }

      // Build backward LSTM
      word_r2l_builders[s].new_graph(*cg);
      word_r2l_builders[s].start_new_sequence();
      for (int i = embeddings.size() -1; i >= 0; --i)
      {
        lstm_backward.push_back(word_r2l_builders[s].add_input(embeddings[i]));
      }

      for (unsigned int i = 0 ; i < lstm_forward.size() ; ++i)
      {
        auto&& e = lstm_backward[lstm_backward.size() - i - 1] +
                   lstm_forward[i];
        embeddings[i] = e;
      }
    }
  }
}

void nn_scorer::set_dropout(float d)
{
  for (auto& b :  word_l2r_builders) {
    b.set_dropout(d);
  }
  for (auto& b :  word_r2l_builders) {
    b.set_dropout(d);
  }
}

void nn_scorer::unset_dropout()
{
  for (auto& b :  word_l2r_builders) {
    b.disable_dropout();
  }
  for (auto& b :  word_r2l_builders) {
    b.disable_dropout();
  }
}
