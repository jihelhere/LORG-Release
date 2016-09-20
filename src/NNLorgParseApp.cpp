
#include "NNLorgParseApp.h"

// #include "utils/Tagger.h"
#include "ParseSolution.h"
#include "utils/LorgConstants.h"

#include "utils/data_parsers/RuleInputParser.h"

#include "SimpleChartCKY.hpp"

#include "training/TreebankFactory.h"


#include "lexicon/WordSignatureFactory.h"


#include "parsers/ParserCKYNN.h"

#include "lexicon/BasicLexicon.h"


#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-local-typedef"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif
#include <boost/program_options.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#if defined(__clang__)
#pragma clang diagnostic pop
#pragma clang diagnostic pop
#endif


// template<>
// Grammar<Rule, Rule, Rule>::Grammar(const std::string& filename)
// {
//   RuleInputParser::read_rulefile(filename,
// 				 lexical_rules,
// 				 unary_rules,
// 				 binary_rules);
// }


NNLorgParseApp::NNLorgParseApp()
    : LorgParseApp(),
      train(true)
{}


NNLorgParseApp::~NNLorgParseApp()
{
}


void write_symboltable(const SymbolTable& st, const std::string& filename)
{
  std::ofstream stream(filename);
  for(unsigned i = 0; i < st.get_symbol_count(); ++i)
  {
    stream << st.get_label_string(i) << std::endl;
  }
}

void collect_rules(const std::vector<PtbPsTree>& trees,
                   std::unordered_set<Rule>& bin,
                   std::unordered_set<Rule>& un,
                   std::unordered_set<Rule>& lex)
{

  // collect counts
  std::unordered_map<Rule, int> rule_counts;
  std::unordered_map<int, int> lhs_counts;

  for (const auto& t : trees)
  {
    std::vector<Production> internals, lexicals;

    t.productions(internals,lexicals);

    for (const auto& i : internals)
    {
      Rule r(i,0,0.0, true);

      if (rule_counts.count(r))
        rule_counts[r]++;
      else
        rule_counts[r] = 1;

      if (lhs_counts.count(r.get_lhs()))
        lhs_counts[r.get_lhs()]++;
      else
        lhs_counts[r.get_lhs()] = 1;

    }

    for(const auto& l : lexicals)
    {
      Rule r(l,0,0.0,true);

      if (rule_counts.count(r))
        rule_counts[r]++;
      else
        rule_counts[r] = 1;

      if (lhs_counts.count(r.get_lhs()))
        lhs_counts[r.get_lhs()]++;
      else
        lhs_counts[r.get_lhs()] = 1;

    }
  }

  /// update scores
  for (const auto& ri : rule_counts)
  {
    auto r = ri.first;
    auto i = ri.second;

    r.set_probability(log (double(i)) - log(double(lhs_counts[r.get_lhs()])));
    if (r.get_rhs().size() == 1)
    {
      if (r.is_lexical())
        lex.insert(r);
      else
        un.insert(r);
    }
    else
      bin.insert(r);
  }
}



int NNLorgParseApp::run_train()
{
  if(verbose) std::clog << "Start learning process.\n";


  // initialise treebank from configuration
  // initialise training grammar from treebank

  Treebank<PtbPsTree> tb(tb_options, verbose);

  BasicLexicon bl(std::shared_ptr<WordSignature>(ws),5);


  auto& trees = tb.get_trees();
  bl.read_lexicon_from_Treebank(trees);
  if(verbose) std::clog << "Read treebank (" << trees.size() << " trees)"  << "\n";

  std::unordered_set<Rule> bin,un;
  std::unordered_set<Rule> lex;


  write_symboltable(SymbolTable::instance_nt(), "NTsymbols.txt");
  write_symboltable(SymbolTable::instance_word(), "Tsymbols.txt");


  collect_rules(trees, bin, un, lex);


  Grammar<Rule,Rule,Rule> grammar;
  grammar.set_rules(std::vector<Rule>(bin.begin(),bin.end()),
                    std::vector<Rule>(un.begin(),un.end()),
                    std::vector<Rule>(lex.begin(),lex.end())
                    );

  // std::cout << grammar << std::endl;


  ParserCKYNN parser(grammar);
  parser.set_word_signature(ws);
  Tagger tagger(&(parser.get_words_to_rules()));
  std::vector<bracketing> brackets;


  // for (const auto& r : grammar.unary_rules) std::cerr << r << std::endl;
  // for (const auto& r : grammar.lexical_rules) std::cerr << r << std::endl;

  cnn::Model m;

  //cnn::SimpleSGDTrainer trainer(&m);
  // trainer.eta_decay = 0.08;
  //cnn::MomentumSGDTrainer trainer(&m);
  cnn::AdagradTrainer trainer(&m);


  int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  // clock_t parse_start = (verbose) ? clock() : 0;


  ParserCKYNN::scorer network(m, trainer);
#define ITERATIONS 100

  for (unsigned iteration = 0; iteration < ITERATIONS; ++iteration)
  {
    if (verbose) std::cerr << "Iteration: " << iteration << std::endl;

    double loss = 0.0;

    std::stringstream oshyp;
    oshyp << "train-hyp-" << iteration << ".mrg";
    std::ofstream outhyp(oshyp.str());

    std::random_shuffle(std::begin(trees), std::end(trees));


    std::stringstream osref;
    osref << "train-ref-" << iteration << ".mrg";
    std::ofstream outref(osref.str());

    for (const auto& tree : trees)
    {
      cnn::ComputationGraph cg;
      network.set_cg(cg);
      network.expressions.clear();


      const auto s = tree.yield();

      PtbPsTree*  best_tree = nullptr;
      if (s.size() <= max_length)
      {
        if (verbose) std::cerr << tree << std::endl;
        //std::cerr << "new tree" << std::endl;
        outref << tree << std::endl;

        std::vector<anchored_binrule_type> anchored_binaries;
        std::vector<anchored_unirule_type> anchored_unaries;
        std::vector<anchored_lexrule_type> anchored_lexicals;

        tree.anchored_productions(anchored_binaries, anchored_unaries, anchored_lexicals);

        network.set_gold(anchored_binaries, anchored_unaries, anchored_lexicals);

        std::vector<Word> words;
        int i = -1;
        std::transform(s.begin(), s.end(), std::back_inserter(words),
                       [&](const std::string& w) {
                         ++i;
                         //if (verbose) std::cerr << w << " ";
                         return Word(w,i,i+1);});
        //if (verbose) std::cerr << '\n';

        // if (verbose)
        //   for (const auto t : anchored_lexicals)
        //   {
        //     std::cerr << std::get<0>(t) << ' ' << std::get<1>(t) << std::endl;
        //   }
        // if (verbose) std::cerr << '\n';

        // tag
        tagger.tag(words, *(parser.get_word_signature()));




        // create and initialise chart
        ParserCKYNN::Chart chart(words,parser.get_nonterm_count(), brackets, network);

        parser.parse(chart, network);

        // get results
        std::cerr << "getting the best tree..." << std::endl;
        best_tree = chart.get_best_tree(start_symbol, 0);

        if (not best_tree)
        {
          std::cerr << "NO SOLUTION" << std::endl;
          outhyp << "(())" << std::endl;
                              continue;
        }

        std::cerr << *best_tree << '\n';
        outhyp << *best_tree << std::endl;

        auto edges = chart.get_best_edges(start_symbol);


        std::cerr << "db1" << std::endl;

        std::vector<cnn::expr::Expression> errs;
        // std::transform(edges.begin(), edges.end(), std::back_inserter(errs),
        //                [&](const decltype(edges)::value_type& ep){return network.expressions[ep];});


        // unsigned c = 0;
        std::for_each(edges.begin(), edges.end(),
                      [&](const decltype(edges)::value_type& ep)
                      {
                        if (network.expressions.count(ep))
                        {
                          errs.push_back(network.expressions[ep]);
                          //                          ++c;
                        }
                      }
                      );

        // std::cerr << "c " << c << std::endl;


        std::cerr << "db2" << std::endl;


        // c = 0;
        std::transform(anchored_lexicals.begin(), anchored_lexicals.end(),
                       std::back_inserter(errs),
                       [&](const anchored_lexrule_type& al)
                       {
                         (void) network.compute_lexical_score(std::get<0>(al),
                                                              &std::get<1>(al));
                         // ++c;
                         return - network.last_expression;
                       }
                       );

        // std::cerr << "c " << c << std::endl;
        std::cerr << "db3" << std::endl;


        std::transform(anchored_unaries.begin(), anchored_unaries.end(),
                       std::back_inserter(errs),
                       [&](const anchored_unirule_type& au)
                       {
                         (void) network.compute_unary_score(std::get<0>(au),
                                                            std::get<1>(au),
                                                            &std::get<2>(au));
                         return - network.last_expression;
                       }
                       );


        std::cerr << "db4" << std::endl;


        std::transform(anchored_binaries.begin(), anchored_binaries.end(),
                       std::back_inserter(errs),
                       [&](const anchored_binrule_type& ab)
                       {
                         (void) network.compute_binary_score(std::get<0>(ab),
                                                             std::get<1>(ab),
                                                             std::get<2>(ab),
                                                             &std::get<3>(ab));
                         return - network.last_expression;
                       }
                       );


        std::cerr << "db5" << std::endl;
        //
        cnn::expr::Expression s = cnn::expr::sum(errs);
        std::cerr << "db5a" << std::endl;
        network.cg->incremental_forward();
        loss += as_scalar(network.cg->get_value(s.i));
        std::cerr << "db5b" << std::endl;
        network.cg->backward(s.i);
        std::cerr << "db5c" << std::endl;

        network.trainer->update(1.0);

        std::cerr << "db6" << std::endl;

        delete best_tree;
        std::cerr << "db7" << std::endl;
      }
      // else
      // {
      //   outhyp << "(())" << std::endl;
      // }
    }
    std::cerr << "ending iteration " << loss << std::endl;
    trainer.status();

    trainer.update_epoch();

    std::stringstream mout;
    mout << "model-" << iteration;
    std::ofstream ms(mout.str());
    boost::archive::text_oarchive oa(ms);

                        oa << m;

  }


  return 0;

}


int NNLorgParseApp::run()
{
  if (train) return run_train();

  if(verbose) std::clog << "Start parsing process.\n";

  // std::vector<Word> s;
  // std::vector< bracketing > brackets;
  // std::string test_sentence;
  // int count = 0;

  // int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  // clock_t parse_start = (verbose) ? clock() : 0;

  // std::vector<std::string> comments;
  //     while(tokeniser->tokenise(*in,test_sentence,s,brackets, comments)) {
  //   clock_t sent_start = (verbose) ? clock() : 0;

  //   // should be "extra-verbose"
  //   // if(verbose) {
  //   //   std::clog << "Tokens: ";
  //   //   for(std::vector<Word>::const_iterator i = s.begin(); i != s.end(); i++)
  //   // 	std::clog << "<" << i->form << ">";
  //   //   std::clog << "\n";
  //   // }


  //   // the pointer to the solution
  //   PtbPsTree*  best_tree = nullptr;

  //   // check length of input sentence
  //   if(s.size() <= max_length) {

  //     // tag
  //     tagger->tag(s, *(parser->get_word_signature()));

  //     // create and initialise chart
  //     ParserCKYBest::Chart chart(s,parser->get_nonterm_count(),brackets);

  //     // parse
  //     parser->parse(chart);

  //     // get results
  //     best_tree = chart.get_best_tree(start_symbol, 0);
  //   }

  //   //FIXME get real prob
  //   std::vector<std::pair<PtbPsTree*, double> > solutions(1, std::make_pair<>(best_tree, 1));

  //   parse_solution p = parse_solution(test_sentence,
  //                                     ++count,
  //                                     s.size(),
  //                                     solutions,
  //                                     (verbose) ? (clock() - sent_start) / double(CLOCKS_PER_SEC) : 0,
  //                                     verbose, comments, false);

  //   *out << unix_parse_solution(p)	 << '\n';

  //   delete best_tree;
  //   s.clear();
  //   brackets.clear();
  // }

  // if(verbose){
  //   std::clog << "overall time: " << (clock() - parse_start) / double(CLOCKS_PER_SEC) << "s\n";
  // }

  return 0; //everything's fine
}



LorgOptions NNLorgParseApp::get_options() const
{

  LorgOptions options(LorgParseApp::get_options());
  options.add_simple_parser_options();
  options.add_treebank_processing_options();
  options.add_lexicon_options();
  options.add_grammar_positionals();




  return options;

}


bool NNLorgParseApp::read_config(ConfigTable& configuration)
{

  if(not LorgParseApp::read_config(configuration))
    return false;

  // if (configuration.exists("train")) {
  //   train = true;

    if (configuration.exists("treebank")) {
      tb_options = TreebankFactory::read_config(configuration);
    }
    else {
        std::cerr << "treebank was not set. Exiting...\n";
        return false;
    }

  // }
  // else {
  //   train = false;
  // }


  // get training grammar
  if(configuration.exists("grammar")) {
    const std::string& training_filename = configuration.get_value<std::string>("grammar");
    if(verbose)
      std::cerr << "Setting grammar to " << training_filename << ".\n";

    //create grammar of rules and associated probs from training file
    //grammar = new Grammar<Rule,Rule,Rule>(training_filename);
    // should check if the grammar is correctly built

  }
  else {
    //std::cerr << "grammar wasn't set." << std::endl;
    //return false;
  }


  // parser = new ParserCKYBest(grammar);

  //tagger = std::auto_ptr<Tagger>(new Tagger(nullptr));
  //tagger->set_word_rules(&(parser->get_words_to_rules()));

  //parser->set_word_signature(WordSignatureFactory::create_wordsignature(configuration));

  ws = WordSignatureFactory::create_wordsignature(configuration);
  return true;
}
