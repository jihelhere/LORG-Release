#include "NNScorer.h"

#include "utils/SymbolTable.h"
#include "rules/Production.h"


#define WORD_EMBEDDING_SIZE 50
#define NT_EMBEDDING_SIZE 20
#define HIDDEN_SIZE 80

nn_scorer::nn_scorer(cnn::Model& m, cnn::Trainer& t) :
    cg(),
    trainer(&t),


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

    expressions()

{}

void nn_scorer::register_last_expression(const Edge * e)
{
  expressions[e] = last_expression;
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



double nn_scorer::compute_lexical_score(int position, const MetaProduction* mp)
{
    auto r = static_cast<const Production*>(mp);

    //std::cerr << *r << std::endl;

    cnn::expr::Expression i = cnn::expr::concatenate({cnn::expr::lookup(*cg,_p_word,r->get_rhs0()),
                                                      cnn::expr::lookup(*cg,_p_nts,r->get_lhs())});

    cnn::expr::Expression W = cnn::expr::parameter(*cg, _p_W_lex);
    cnn::expr::Expression b = cnn::expr::parameter(*cg, _p_b_lex);
    cnn::expr::Expression o = cnn::expr::parameter(*cg, _p_o_lex);

    last_expression = o * cnn::expr::tanh(W*i + b);

    //std::cerr << "before as_scalar" << std::endl;
    cg->incremental_forward();
    double v = as_scalar(cg->get_value(last_expression.i));
    //std::cerr << v << std::endl;
    //std::cerr << "after as_scalar" << std::endl;

    if (gold and not anchored_lexicals.count(std::make_tuple(position,*r)))
      v += 1.0;

    return v;
  }

double nn_scorer::compute_unary_score(int begin, int end, const MetaProduction* mp)
  {

    auto r = static_cast<const Production*>(mp);

    cnn::expr::Expression i = cnn::expr::concatenate({cnn::expr::lookup(*cg,_p_nts,r->get_rhs0()),
                                                       cnn::expr::lookup(*cg,_p_nts,SymbolTable::instance_nt().get_symbol_count()),
                                                       cnn::expr::lookup(*cg,_p_nts,r->get_lhs())});


    cnn::expr::Expression W = cnn::expr::parameter(*cg, _p_W_int);
    cnn::expr::Expression b = cnn::expr::parameter(*cg, _p_b_int);
    cnn::expr::Expression o = cnn::expr::parameter(*cg, _p_o_int);

    last_expression = o * cnn::expr::tanh(W*i + b);

    cg->incremental_forward();
    // return 0.0;
    double v = as_scalar(cg->get_value(last_expression.i));

    if (gold and not anchored_unaries.count(std::make_tuple(begin,end,*r)))
      v += 1.0;


   return v;
  }

double nn_scorer::compute_binary_score(int s, int e, int m, const MetaProduction* mp)
  {
    auto r = static_cast<const Production*>(mp);

    cnn::expr::Expression i = cnn::expr::concatenate({cnn::expr::lookup(*cg,_p_nts,r->get_rhs0()),
                                                      cnn::expr::lookup(*cg,_p_nts,r->get_rhs1()),
                                                      cnn::expr::lookup(*cg,_p_nts,r->get_lhs())});

    cnn::expr::Expression W = cnn::expr::parameter(*cg, _p_W_int);
    cnn::expr::Expression b = cnn::expr::parameter(*cg, _p_b_int);
    cnn::expr::Expression o = cnn::expr::parameter(*cg, _p_o_int);

    last_expression = o * cnn::expr::tanh(W*i + b);

    cg->incremental_forward();
    double v = as_scalar(cg->get_value(last_expression.i));

    if (gold and not anchored_binaries.count(std::make_tuple(s,e,m,*r)))
      v += 1.0;

    return v;
    //return 0.0;
  }
