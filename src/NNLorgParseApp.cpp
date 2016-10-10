
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
  unix_parse_solution::init();
  json_parse_solution::init();

  if(verbose) std::clog << "Start learning process.\n";

  // initialise treebank from configuration
  // initialise training grammar from treebank

  Treebank<PtbPsTree> tb(tb_options, verbose);

  BasicLexicon bl(std::shared_ptr<WordSignature>(ws),cutoff);

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

  for (const auto& b : grammar.binary_rules)
  {
    std::cout << b << std::endl;
  }
  for (const auto& u : grammar.unary_rules)
  {
    std::cout << u << std::endl;
  }
  for (const auto& l : grammar.lexical_rules)
  {
    std::cout << l << std::endl;
  }


  ParserCKYNN parser(grammar);

#ifdef USE_THREADS
  parser.set_nbthreads(nbthreads);
#endif

  parser.set_word_signature(ws);
  Tagger tagger(&(parser.get_words_to_rules()));
  std::vector<bracketing> brackets;


  // for (const auto& r : grammar.unary_rules) std::cerr << r << std::endl;
  // for (const auto& r : grammar.lexical_rules) std::cerr << r << std::endl;

  cnn::Model m;

  //cnn::SimpleSGDTrainer trainer(&m);
  // trainer.eta_decay = 0.08;
  //cnn::MomentumSGDTrainer trainer(&m);
  cnn::AdamTrainer trainer(&m);


  int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  // clock_t parse_start = (verbose) ? clock() : 0;


  ParserCKYNN::scorer network(m);

#define ITERATIONS 30
#define MINI_BATCH 10

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




    auto tree_idx = 0U;
    while (tree_idx < trees.size())
    {

      cnn::ComputationGraph cg;
      network.set_cg(cg);
      network.clear();

      //std::cerr << "set rule scores" << std::endl;
      network.precompute_rule_expressions(grammar.binary_rules, grammar.unary_rules);

      std::vector<cnn::expr::Expression> errs;

      for (auto mini_batch_idx = 0; mini_batch_idx < MINI_BATCH and tree_idx < trees.size(); ++ mini_batch_idx)
      {
        const auto& tree = trees[tree_idx];
        ++ tree_idx;

        const auto s = tree.yield();

        PtbPsTree*  best_tree = nullptr;
        if (s.size() <= max_length)
        {
          if (verbose) std::cerr << tree << std::endl;
          //std::cerr << "new tree" << std::endl;
          outref << tree << std::endl;

          std::vector<anchored_binrule_type> anchored_binaries_vec;
          std::vector<anchored_unirule_type> anchored_unaries_vec;
          std::vector<anchored_lexrule_type> anchored_lexicals_vec;

          tree.anchored_productions(anchored_binaries_vec, anchored_unaries_vec, anchored_lexicals_vec);
          auto anchored_binaries = std::unordered_set<decltype(anchored_binaries_vec)::value_type>(anchored_binaries_vec.begin(),
                                                                                                   anchored_binaries_vec.end());
          auto anchored_unaries = std::unordered_set<decltype(anchored_unaries_vec)::value_type>(anchored_unaries_vec.begin(),
                                                                                                 anchored_unaries_vec.end());
          auto anchored_lexicals = std::unordered_set<decltype(anchored_lexicals_vec)::value_type>(anchored_lexicals_vec.begin(),
                                                                                                   anchored_lexicals_vec.end());


          network.set_gold(anchored_binaries_vec, anchored_unaries_vec, anchored_lexicals_vec);

          std::vector<Word> words;
          int i = -1;
          std::transform(s.begin(), s.end(), std::back_inserter(words),
                         [&](const std::string& w) {
                           ++i;
                           return Word(w,i,i+1);});


          //std::cerr << "tag" << std::endl;
          tagger.tag(words, *(parser.get_word_signature()));


          //std::cerr << "set words" << std::endl;
          network.set_words(words);

          //std::cerr << "set embeddings" << std::endl;
          network.precompute_embeddings();
          //network.precompute_span_expressions(grammar.lhs_int_set);

          // create and initialise chart
          //std::cerr << "chart" << std::endl;
          ParserCKYNN::Chart chart(words,parser.get_nonterm_count(), brackets, network);

          //std::cerr << "parse" << std::endl;
          parser.parse(chart, network);

          // get results
          std::cerr << "getting the best tree..." << std::endl;
          best_tree = chart.get_best_tree(start_symbol, 0);


          std::vector<anchored_binrule_type> best_anchored_binaries_vec;
          std::vector<anchored_unirule_type> best_anchored_unaries_vec;
          std::vector<anchored_lexrule_type> best_anchored_lexicals_vec;

          if (not best_tree)
          {
            fprintf(stderr, "NO SOLUTION\n");
            outhyp << "(())" << std::endl;
          }
          else
          {
            std::cerr << *best_tree << '\n';
            outhyp << *best_tree << std::endl;

            best_tree->anchored_productions(best_anchored_binaries_vec,
                                            best_anchored_unaries_vec,
                                            best_anchored_lexicals_vec);
          }


          auto best_anchored_binaries = std::unordered_set<decltype(best_anchored_binaries_vec)::value_type>(best_anchored_binaries_vec.begin(),
                                                                                                             best_anchored_binaries_vec.end());
          auto best_anchored_unaries = std::unordered_set<decltype(best_anchored_unaries_vec)::value_type>(best_anchored_unaries_vec.begin(),
                                                                                                           best_anchored_unaries_vec.end());
          auto best_anchored_lexicals = std::unordered_set<decltype(best_anchored_lexicals_vec)::value_type>(best_anchored_lexicals_vec.begin(),
                                                                                                             best_anchored_lexicals_vec.end());


          //binary rules
          for (const auto& ref_anc_bin : anchored_binaries)
          {
            if (not best_anchored_binaries.count(ref_anc_bin))
            {
              errs.push_back( - network.rule_expression(std::get<3>(ref_anc_bin).get_lhs(),
                                                        std::get<3>(ref_anc_bin).get_rhs0(),
                                                        std::get<3>(ref_anc_bin).get_rhs1()
                                                        ));
            }
            // else
            // {
            //   std::cerr << std::get<3>(ref_anc_bin)
            //             << " (" << std::get<0>(ref_anc_bin) << "," << std::get<1>(ref_anc_bin) << "," << std::get<2>(ref_anc_bin) << ")"
            //             << "was correctly retrieved";
            // }
          }

          for (const auto& best_anc_bin : best_anchored_binaries)
          {
            if (not anchored_binaries.count(best_anc_bin))
            {
              errs.push_back(network.rule_expression(std::get<3>(best_anc_bin).get_lhs(),
                                                     std::get<3>(best_anc_bin).get_rhs0(),
                                                     std::get<3>(best_anc_bin).get_rhs1()
                                                     ));
            }
            // else
            // {
            //   std::cerr << std::get<3>(best_anc_bin)
            //             << " (" << std::get<0>(best_anc_bin) <<"," << std::get<1>(best_anc_bin) << "," << std::get<2>(best_anc_bin) << ")"
            //             << "was correctly retrieved";
            // }
          }


          //unary rules
          for (const auto& ref_anc_un : anchored_unaries)
          {
            if (not best_anchored_unaries.count(ref_anc_un))
            {
              errs.push_back( - network.rule_expression(std::get<2>(ref_anc_un).get_lhs(),
                                                        std::get<2>(ref_anc_un).get_rhs0(),
                                                        SymbolTable::instance_nt().get_symbol_count()
                                                        ));
            }
            // else
            // {
            //   std::cerr << std::get<2>(ref_anc_un)
            //             << " (" << std::get<0>(ref_anc_un) <<"," <<  std::get<1>(ref_anc_un) << ")"
            //             << "was correctly retrieved";
            // }
          }

          for (const auto& best_anc_un : best_anchored_unaries)
          {
            if (not anchored_unaries.count(best_anc_un))
            {
              errs.push_back(network.rule_expression(std::get<2>(best_anc_un).get_lhs(),
                                                     std::get<2>(best_anc_un).get_rhs0(),
                                                     SymbolTable::instance_nt().get_symbol_count()
                                                     ));
            }
            // else
            // {
            //   std::cerr << std::get<2>(best_anc_un)
            //             << " (" << std::get<0>(best_anc_un) <<"," <<  std::get<1>(best_anc_un) << ")"
            //             << "was correctly retrieved";
            // }
          }


          //lexical rules
          for (const auto& ref_anc_lex : anchored_lexicals)
          {
            if (not best_anchored_lexicals.count(ref_anc_lex))
            {
              errs.push_back( - network.lexical_rule_expression(std::get<1>(ref_anc_lex).get_lhs(),
                                                                std::get<0>(ref_anc_lex)
                                                                ));
            }
            // else
            // {
            //   std::cerr << std::get<1>(ref_anc_lex)
            //             << " (" << std::get<0>(ref_anc_lex) << ")"
            //             << "was correctly retrieved";
            // }
          }

          for (const auto& best_anc_lex : best_anchored_lexicals)
          {
            if (not anchored_lexicals.count(best_anc_lex))
            {
              errs.push_back(network.lexical_rule_expression(std::get<1>(best_anc_lex).get_lhs(),
                                                             std::get<0>(best_anc_lex)
                                                             ));
            }
            // else
            // {
            //   std::cerr << std::get<1>(best_anc_lex)
            //             << " (" << std::get<0>(best_anc_lex) << ")"
            //             << "was correctly retrieved";
            // }
          }


          // expp.push_back(network.span_expression(std::get<2>(au).get_lhs(),
          //                                        std::get<0>(au)));

          // expp.push_back(network.span_expression(std::get<3>(ab).get_lhs(),
          //                                        std::get<0>(ab)));

        }
        if (best_tree)
          delete best_tree;
      }

      if (not errs.empty())
      {
        cnn::expr::Expression s = cnn::expr::sum(errs);
        // std::cerr << "db5a" << std::endl;
        network.cg->incremental_forward();
        loss += as_scalar(network.cg->get_value(s.i));
        // std::cerr << "db5b" << std::endl;
        network.cg->backward(s.i);
        // std::cerr << "db5c" << std::endl;

        trainer.update(1.0);
      }
      // std::cerr << "db6" << std::endl;

      // std::cerr << "db7" << std::endl;
    }

    std::cerr << "ending iteration " << loss << std::endl;
    trainer.status();

    trainer.update_epoch();

    std::stringstream mout;
    mout << "model-" << iteration;
                        std::ofstream ms(mout.str());
                        boost::archive::text_oarchive oa(ms);

                        oa << m;

                              ///////////// process dev file
                              {

                                std::vector<Word> s;
                                std::vector< bracketing > brackets;
                                std::string test_sentence;
                                int count = 0;

                                int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

                                //clock_t parse_start = (verbose) ? clock() : 0;

                                std::stringstream devss;
                                devss << "dev-hyp-" << iteration;
                                std::ofstream devstream(devss.str());

                                //rewind test-input
                                delete in;
                                in = new std::ifstream(in_filename.c_str());

                                std::vector<std::string> comments;
                                while(tokeniser->tokenise(*in,test_sentence,s,brackets, comments)) {
                                  clock_t sent_start = (verbose) ? clock() : 0;

                                  // should be "extra-verbose"
                                  if(verbose) {
                                    std::clog << "Tokens: ";
                                    for(std::vector<Word>::const_iterator i = s.begin(); i != s.end(); i++)
                                      std::clog << "<" << i->get_form() << ">";
                                    std::clog << "\n";
                                  }

                                  cnn::ComputationGraph cg;
                                  network.set_cg(cg);
                                  network.unset_gold();
                                  network.clear();

                                  // the pointer to the solution
                                  PtbPsTree*  best_tree = nullptr;

                                  // check length of input sentence
                                  if(s.size() <= max_length) {

                                    // tag
                                    std::cerr << "tag" << std::endl;
                                    tagger.tag(s, *(parser.get_word_signature()));
                                    network.set_words(s);
                                    network.precompute_embeddings();
                                    // should save values once and for all
                                    network.precompute_rule_expressions(grammar.binary_rules, grammar.unary_rules);
                                    //network.precompute_span_expressions(grammar.lhs_int_set);

                                    // create and initialise chart
                                    std::cerr << "chart" << std::endl;
                                    ParserCKYNN::Chart chart(s,parser.get_nonterm_count(),brackets, network);

                                    // parse
                                    std::cerr << "parse" << std::endl;
                                    parser.parse(chart, network);

                                    // get results
                                    std::cerr << "results" << std::endl;
                                    best_tree = chart.get_best_tree(start_symbol, 0);
                                  }

                                  //FIXME get real prob
                                  std::vector<std::pair<PtbPsTree*, double> > solutions = {{best_tree, 1.0}};


                                  parse_solution * p_typed =
                                      parse_solution::factory.create_object(parse_solution::UNIX,
                                                                            parse_solution(test_sentence, ++count,
                                                                                           s.size(), solutions,
                                                                                           (verbose) ? (clock() - sent_start) / double(CLOCKS_PER_SEC) : 0,
                                                                                           verbose, comments,false)
                                                                            );
                                  p_typed->print(devstream);
                                  delete p_typed;


                                  delete best_tree;
                                  s.clear();
                                  brackets.clear();


                                }
                              }
  }

        return 0;

}


  int NNLorgParseApp::run()
  {
    if (train) return run_train();
    // else
    //   return run_parse();



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


  cutoff = configuration.get_value<unsigned>("unknown-word-cutoff");

  ws = WordSignatureFactory::create_wordsignature(configuration);
  return true;
}
