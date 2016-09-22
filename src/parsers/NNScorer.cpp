#include "NNScorer.h"

#include "utils/SymbolTable.h"
#include "rules/Production.h"

#define WORD_EMBEDDING_SIZE 50
#define NT_EMBEDDING_SIZE 20
#define HIDDEN_SIZE 80

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
    _p_word(m.add_lookup_parameters(SymbolTable::instance_word().get_symbol_count(),
                                    {WORD_EMBEDDING_SIZE})),
    _p_nts(m.add_lookup_parameters(SymbolTable::instance_nt().get_symbol_count()+1,
                                   {NT_EMBEDDING_SIZE})),
    rules_expressions(),
    expressions()

{}

void nn_scorer::register_expression(const Edge * ep,
                                    const cnn::expr::Expression& expp)
{
  {
    tbb::spin_mutex::scoped_lock lock(expression_mutex);
    expressions[ep] = expp;
  }
}


void nn_scorer::clear()
{
  rules_expressions.clear();
  expressions.clear();

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
                                        cnn::expr::Expression* expp)
{
    auto r = static_cast<const Production*>(mp);

    double v;

    {
      tbb::spin_mutex::scoped_lock lock(rule_mutex);
    if (rules_expressions.count(r))
    {
      v = as_scalar(cg->get_value(rules_expressions[r].i));
      *expp = rules_expressions[r];
    }
    else
    {

      cnn::expr::Expression i = cnn::expr::concatenate({cnn::expr::lookup(*cg,_p_word,r->get_rhs0()),
                                                        cnn::expr::lookup(*cg,_p_nts,r->get_lhs())});

      cnn::expr::Expression W = cnn::expr::parameter(*cg, _p_W_lex);
      cnn::expr::Expression b = cnn::expr::parameter(*cg, _p_b_lex);
      cnn::expr::Expression o = cnn::expr::parameter(*cg, _p_o_lex);

      cnn::expr::Expression out = o * cnn::expr::rectify(W*i + b);

      //std::cerr << "before as_scalar" << std::endl;
      cg->incremental_forward();
      v = as_scalar(cg->get_value(out.i));
      //std::cerr << v << std::endl;
      //std::cerr << "after as_scalar" << std::endl;

      rules_expressions[r] = out;
      *expp = rules_expressions[r];
    }
    }


    if (gold and not anchored_lexicals.count(std::make_tuple(position,*r)))
    {
      //std::cerr << "wrong rule at position: " << position << " " << *r << std::endl;
      v += 1.0;
    }
    // else
    // {
    //   std::cerr << "correct rule at position: " << position << " " << *r << std::endl;
    // }

    return v;
    //return 0.0;
  }

double nn_scorer::compute_unary_score(int begin, int end, const MetaProduction* mp,
                                      cnn::expr::Expression *expp)
  {

    auto r = static_cast<const Production*>(mp);

    double v;

    {
      tbb::spin_mutex::scoped_lock lock(rule_mutex);
    if (rules_expressions.count(r))
    {
      v = as_scalar(cg->get_value(rules_expressions[r].i));
      *expp = rules_expressions[r];
    }
    else
    {

      cnn::expr::Expression i = cnn::expr::concatenate({cnn::expr::lookup(*cg,_p_nts,r->get_rhs0()),
                                                        cnn::expr::lookup(*cg,_p_nts,SymbolTable::instance_nt().get_symbol_count()),
                                                        cnn::expr::lookup(*cg,_p_nts,r->get_lhs())});


    cnn::expr::Expression W = cnn::expr::parameter(*cg, _p_W_int);
    cnn::expr::Expression b = cnn::expr::parameter(*cg, _p_b_int);
    cnn::expr::Expression o = cnn::expr::parameter(*cg, _p_o_int);

    cnn::expr::Expression out = o * cnn::expr::rectify(W*i + b);

    cg->incremental_forward();
    // return 0.0;
    v = as_scalar(cg->get_value(out.i));

    rules_expressions[r] = out;
    *expp = rules_expressions[r];

    }
    }

    if (gold and not anchored_unaries.count(std::make_tuple(begin,end,*r)))
      v += 1.0;

   return v;


   //return 0.0;
  }

double nn_scorer::compute_binary_score(int s, int e, int m, const MetaProduction* mp,
                                       cnn::expr::Expression* expp)
{
    auto r = static_cast<const Production*>(mp);

    double v;

    {
      tbb::spin_mutex::scoped_lock lock(rule_mutex);
      if (rules_expressions.count(r))
      {
        v = as_scalar(cg->get_value(rules_expressions[r].i));
        *expp = rules_expressions[r];
      }
      else
      {
        cnn::expr::Expression i = cnn::expr::concatenate({cnn::expr::lookup(*cg,_p_nts,r->get_rhs0()),
                                                          cnn::expr::lookup(*cg,_p_nts,r->get_rhs1()),
                                                          cnn::expr::lookup(*cg,_p_nts,r->get_lhs())});

        cnn::expr::Expression W = cnn::expr::parameter(*cg, _p_W_int);
        cnn::expr::Expression b = cnn::expr::parameter(*cg, _p_b_int);
        cnn::expr::Expression o = cnn::expr::parameter(*cg, _p_o_int);

        cnn::expr::Expression out = o * cnn::expr::rectify(W*i + b);

        cg->incremental_forward();
        v = as_scalar(cg->get_value(out.i));

        rules_expressions[r] = out;
        *expp = rules_expressions[r];
      }
    }

    if (gold and not anchored_binaries.count(std::make_tuple(s,e,m,*r)))
      v += 1.0;

    return v;
    //  return 0.0;
  }
