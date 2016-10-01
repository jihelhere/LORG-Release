#include "NNScorer.h"

#include "utils/SymbolTable.h"
#include "rules/Production.h"

#define WORD_EMBEDDING_SIZE 60
#define NT_EMBEDDING_SIZE 20
#define HIDDEN_SIZE 40


#ifdef USE_THREADS
#include <tbb/mutex.h>
#include <tbb/spin_mutex.h>


tbb::spin_mutex rule_mutex;
tbb::spin_mutex expression_mutex;
#endif

nn_scorer::nn_scorer(cnn::Model& m) :
    cg(nullptr),

    _p_W_int(m.add_parameters({HIDDEN_SIZE, NT_EMBEDDING_SIZE*3})),
    _p_b_int(m.add_parameters({HIDDEN_SIZE})),
    _p_o_int(m.add_parameters({1,HIDDEN_SIZE})),

    _p_W_lex(m.add_parameters({HIDDEN_SIZE, NT_EMBEDDING_SIZE+WORD_EMBEDDING_SIZE})),
    _p_b_lex(m.add_parameters({HIDDEN_SIZE})),
    _p_o_lex(m.add_parameters({1,HIDDEN_SIZE})),

    _p_W_span(m.add_parameters({HIDDEN_SIZE, WORD_EMBEDDING_SIZE*1+NT_EMBEDDING_SIZE})),
    _p_b_span(m.add_parameters({HIDDEN_SIZE})),
    _p_o_span(m.add_parameters({1,HIDDEN_SIZE})),


    _p_word(m.add_lookup_parameters(SymbolTable::instance_word().get_symbol_count()+1,
                                    {WORD_EMBEDDING_SIZE})),
    _p_nts(m.add_lookup_parameters(SymbolTable::instance_nt().get_symbol_count()+1,
                                   {NT_EMBEDDING_SIZE})),
    rules_expressions(),
    spans_expressions(),
    edges_expressions()

{}

void nn_scorer::register_expression(const Edge * ep,
                                    std::vector<cnn::expr::Expression>& expv)
{
  {
    tbb::spin_mutex::scoped_lock lock(expression_mutex);
    edges_expressions[ep] = expv;
  }
}


void nn_scorer::clear()
{
  rules_expressions.clear();
  spans_expressions.clear();
  edges_expressions.clear();

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



double nn_scorer::compute_lexical_score(int position, const MetaProduction* mp,
                                        std::vector<cnn::expr::Expression>& expv)
{
    auto r = static_cast<const Production*>(mp);

    //std::cerr << *r << std::endl;

    double v = as_scalar(cg->get_value(rules_expressions.at(r).i));
    expv.push_back(rules_expressions[r]);

    //std::cerr << "after" << std::endl;


    if (gold and not anchored_lexicals.count(std::make_tuple(position,*r)))
    {
      v += 1.0;
    }
    return v;
}


double
nn_scorer::compute_internal_rule_score(const Production* r, std::vector<cnn::expr::Expression>& expp)
{

  double v = as_scalar(cg->get_value(rules_expressions.at(r).i));
  expp.push_back(rules_expressions[r]);

  return v;
}

double
nn_scorer::compute_internal_span_score(int begin,
                                       int end,
                                       int medium,
                                       int lhs,
                                       std::vector<cnn::expr::Expression>& e)
{
  //  return 0.0;



  //  auto t = std::make_tuple(begin,end,medium,lhs);
  auto t = std::make_tuple(begin,lhs);

  double v = as_scalar(cg->get_value(spans_expressions[t].i));
  e.push_back(spans_expressions[t]);

  return v;
}



void nn_scorer::set_words(const std::vector<Word>& w)
{
  words = w;
}

double nn_scorer::compute_unary_score(int begin, int end, const MetaProduction* mp,
                                      std::vector<cnn::expr::Expression>& expv)
  {

    auto r = static_cast<const Production*>(mp);

    //std::cerr << *r << std::endl;

    double v = compute_internal_rule_score(r, expv);

    //std::cerr << "b4 un" << std::endl;



    //std::cerr << "un: " << begin << " " << words[begin]  << (end -1) << " " << words[end-1] << std::endl;

    v+= compute_internal_span_score(begin,
                                    end - 1,
                                    -1,
                                    r->get_lhs(),
                                    expv);
    //       std::cerr << "after un" << std::endl;


    if (gold and not anchored_unaries.count(std::make_tuple(begin,end,*r)))
      v += 1.0;
    // else
    // {
    //   std::cerr << "gold unary: " << *r << " " << begin << " " << end << std::endl;
    // }


    return v;

  }

double nn_scorer::compute_binary_score(int s, int e, int m, const MetaProduction* mp,
                                       std::vector<cnn::expr::Expression>& expp)
{
    auto r = static_cast<const Production*>(mp);

    double v = compute_internal_rule_score(r, expp);

    //    std::cerr << "b4 bin" << std::endl;

    //std::cerr << "bin: " << s << " " << words[s]  << (e-1) << " " << words[e-1]  << m << " " << words[m] << std::endl;


    v+= compute_internal_span_score(s,
                                    e - 1,
                                    m,
                                    r->get_lhs(),
                                    expp);
    //       std::cerr << "after bin" << std::endl;


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
    rules_expressions[&r] = rule_expression(r.get_lhs(),
                                            r.get_rhs0(),
                                            r.get_rhs1());
  }

  for (const auto& r : urules)
  {
    rules_expressions[&r] = rule_expression(r.get_lhs(),
                                            r.get_rhs0(),
                                            SymbolTable::instance_nt().get_symbol_count());
  }

  for (const auto& w: words)
  {
    for (const auto mppt : w.get_rules())
    {
      auto ppt = static_cast<const Production*>(mppt);
      if (not rules_expressions.count(ppt))
        rules_expressions[ppt] = lexical_rule_expression(ppt->get_lhs(), ppt->get_rhs0());
    }
  }

  cg->forward();
}


void nn_scorer::precompute_span_expressions(const std::unordered_set<int>& lhs_set)
{
  for (auto lhs : lhs_set)
  for (unsigned i = 0; i < words.size(); ++i)
  {
    const auto& word = words[i];

    spans_expressions[std::make_tuple(i,lhs)] = span_expression(lhs, word.get_id());

  }
  cg->incremental_forward();
}

cnn::expr::Expression nn_scorer::span_expression(int lhs, int first_word_id)
{
   cnn::expr::Expression i = cnn::expr::concatenate({cnn::expr::lookup(*cg,_p_word, first_word_id),
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


cnn::expr::Expression nn_scorer::lexical_rule_expression(int lhs, int rhs0)
{
  cnn::expr::Expression i = cnn::expr::concatenate({cnn::expr::lookup(*cg, _p_word, rhs0),
                                                    cnn::expr::lookup(*cg, _p_nts, lhs)});

      cnn::expr::Expression W = cnn::expr::parameter(*cg, _p_W_lex);
      cnn::expr::Expression b = cnn::expr::parameter(*cg, _p_b_lex);
      cnn::expr::Expression o = cnn::expr::parameter(*cg, _p_o_lex);

      return o * cnn::expr::rectify(W*i + b);
}
