#include "NNScorer.h"

#include "utils/SymbolTable.h"
#include "rules/Production.h"

#include "utils/hash_impl.h"

#include <mutex>

namespace dy = dynet;
namespace de = dy::expr;

bool nn_scorer::train_mode = false;

cfg_rule_representation nn_scorer::cfg;

nn_scorer::nn_scorer(bool init_global,
                     dy::Model& m, unsigned lex_level, unsigned sp_level,
                     unsigned word_embedding_size,
                     unsigned nt_embedding_size,
                     unsigned hidden_size,
                     unsigned lstm_hidden_size,
                     bool char_emb,
                     bool span_midp
                     ) :
    lexical_level(lex_level),
    span_level(sp_level),
    words(nullptr),
    use_char_emb(char_emb),
    use_span_midpoints(span_midp)
{
  std::lock_guard<std::mutex> guard(computation_attachment::cg_mutex);


  if (init_global)
    cfg = cfg_rule_representation(init_global, m,
                                  nt_embedding_size, hidden_size);

  word_repr = use_char_emb ?
              static_cast<word_representation*>
              (new char_word_representation(init_global, m, word_embedding_size))
              :
              static_cast<word_representation*>
              (new simple_word_representation(init_global, m, word_embedding_size));

  lex_repr = lexical_level == 0 ?
             static_cast<lexicon_representation*>
             (new simple_lexicon_representation(word_repr))
             :
             static_cast<lexicon_representation*>
             (new bilstm_lexicon_representation(word_repr, init_global, m,
                                                word_embedding_size, lstm_hidden_size,
                                                lexical_level));

  lexrule_repr = lexical_level == 0 ?
                 static_cast<lexical_rule_representation*>
                 (new simple_lexical_rule(init_global, m,
                                          nt_embedding_size,
                                          word_embedding_size, hidden_size,
                                          lex_repr,
                                          &cfg))
                 :
                 static_cast<lexical_rule_representation*>
                 (new bilstm_lexical_rule(init_global,m,
                                          nt_embedding_size,
                                          hidden_size, lstm_hidden_size,
                                          lex_repr,
                                          &cfg));


  span_repr = span_level == 0 ?
              static_cast<span_representation*>
              (new empty_span_representation())
              :
              static_cast<span_representation*>
              (new all_span_representation(init_global, m,
                                           span_level,
                                           lexical_level == 0 ? word_embedding_size : lstm_hidden_size,
                                           nt_embedding_size,
                                           hidden_size,
                                           use_span_midpoints,
                                           lex_repr,&cfg));
}

void nn_scorer::clear()
{
  span_repr->clear();
  // span_scores_bin.clear();
  // span_scores_un.clear();
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
    std::lock_guard<std::mutex> guard(computation_attachment::cg_mutex);
    v = as_scalar(cg->get_value(e.i));
  }

  if (gold && !anchored_lexicals->count(std::make_tuple(position,*r)))
  {
    v += 1.0;
  }
  return v;
}


double nn_scorer::compute_internal_rule_score(const Production* r)
{
  return cfg.compute_internal_rule_score(r);
}


// same code as span_expression ??

double
nn_scorer::compute_internal_span_score(int begin,
                                       int end,
                                       int medium,
                                       int lhs, int rhs0)
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


  if (medium >= 0) // binary rule
  {
    double v = span_repr->get_span_score_lhs_begin(lhs, begin) +
               span_repr->get_span_score_lhs_end(lhs, end) +
               span_repr->get_span_score_lhs_split(lhs, medium);

    // v += span_repr->get_span_score_rhs0_begin(rhs0, begin) +
    //      span_repr->get_span_score_rhs0_end(rhs0, end) +
    //      span_repr->get_span_score_rhs0_split(rhs0, split);

    if (not use_span_midpoints) medium = 0;

    return v + span_repr->get_span_score_bin_info(begin,end, medium, lhs_info);

  }
  else
  {
    double v = span_repr->get_span_score_lhs_begin(lhs, begin) +
               span_repr->get_span_score_lhs_end(lhs, end);

    // v += span_repr->get_span_score_rhs0_begin(rhs0, begin) +
    //      span_repr->get_span_score_rhs0_end(rhs0, end);

    return v + span_repr->get_span_score_una_info(begin,end, lhs_info);
  }

}

void nn_scorer::set_words(const std::vector<Word>& w)
{
  words = &w;
}

double nn_scorer::compute_unary_score(int begin, int end, const MetaProduction* mp)
{

  auto r = static_cast<const Production*>(mp);

  double v = compute_internal_rule_score(r);

  if (span_level > 0) v+= compute_internal_span_score(begin, end - 1, -1, r->get_lhs(),
                                                      r->get_rhs0());

  if (gold && !anchored_unaries->count(std::make_tuple(begin,end,*r)))
    v += 1.0;

  return v;

}

double nn_scorer::compute_binary_score(int s, int e, int m, const MetaProduction* mp)
{
  auto r = static_cast<const Production*>(mp);

  double v = compute_internal_rule_score(r);

  if (span_level > 0) v+= compute_internal_span_score(s, e - 1, m, r->get_lhs(), r->get_rhs0());

  if (gold and not anchored_binaries->count(std::make_tuple(s,e,m,*r))) v += 1.0;

  return v;
}

void nn_scorer::precompute_rule_expressions(const std::vector<Rule>& brules,
                                            const std::vector<Rule>& urules)
{
  cfg.precompute_rule_expressions(brules, urules, train_mode);
}

void nn_scorer::precompute_span_expressions(const std::vector<int>& lhs_int,
                                            const std::vector<int>& rhs0_int,
                                            const std::vector<int>& rhs1_int)
{
  span_repr->precompute_span_expressions(lhs_int, rhs0_int, rhs1_int, *words, train_mode);
}

de::Expression nn_scorer::lexical_rule_expression(int lhs, unsigned word_idx)
{
  return lexrule_repr->lexical_rule_expression(lhs, word_idx);
}

void nn_scorer::precompute_embeddings()
{
  lex_repr->precompute_embeddings(*words);
}


void nn_scorer::set_dropout(float d)
{
  lex_repr->set_dropout(d);
}

void nn_scorer::unset_dropout()
{
  lex_repr->unset_dropout();
}
