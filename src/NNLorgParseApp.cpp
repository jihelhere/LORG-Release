#include <algorithm>
#include <thread>



#include "NNLorgParseApp.h"


#include "ParseSolution.h"
#include "utils/LorgConstants.h"
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

NNLorgParseApp::NNLorgParseApp()
    : LorgParseApp(), train(true)
{}

NNLorgParseApp::~NNLorgParseApp()
{}


void write_symboltable(const SymbolTable& st, const std::string& filename)
{
  std::ofstream stream(filename);
  for(unsigned i = 0; i < st.get_symbol_count(); ++i)
  {
    stream << st.get_label_string(i) << std::endl;
  }
}


// extract a MLE pcfg from (binarized) trees
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



// TODO: separate parsing and training

std::pair<
  std::pair<std::vector<dynet::expr::Expression>,std::vector<dynet::expr::Expression>>,
  std::pair<std::string,std::string>>
NNLorgParseApp::train_instance(const PtbPsTree& tree,
                               const ParserCKYNN& parser,
                               const Tagger& tagger,
                               ParserCKYNN::scorer& network,
                               int start_symbol)
{
  network.span_scores.clear();

  std::vector<dynet::expr::Expression> local_corrects, local_errs;
  std::stringstream ssref,sshyp;

  const auto s = tree.yield();

  PtbPsTree*  best_tree = nullptr;
  if (s.size() <= max_length)
  {
    //std::cerr << tree << std::endl;
    ssref << tree << std::endl;

    std::unordered_set<anchored_binrule_type> anchored_binaries;
    std::unordered_set<anchored_unirule_type> anchored_unaries;
    std::unordered_set<anchored_lexrule_type> anchored_lexicals;

    tree.anchored_productions(anchored_binaries, anchored_unaries, anchored_lexicals);

    network.set_gold(anchored_binaries, anchored_unaries, anchored_lexicals);

    std::vector<Word> words;
    for (auto i = 0U; i < s.size(); ++i)
    {
      words.push_back(Word(s[i],i,i+1));
    }

    //std::cerr << "tag" << std::endl;
    tagger.tag(words, *(parser.get_word_signature()));


    //std::cerr << "set words" << std::endl;
    network.set_words(words);

    //std::cerr << "set embeddings" << std::endl;
    network.precompute_embeddings();

    if (span_level > 0) network.precompute_span_expressions();

    // create and initialise chart
    //std::cerr << "chart" << std::endl;
    std::vector<bracketing> brackets;

    ParserCKYNN::Chart chart(words,parser.get_nonterm_count(), brackets, network);


    //std::cerr << "parse" << std::endl;
    parser.parse(chart, network);

    best_tree = chart.get_best_tree(start_symbol, 0);

    std::unordered_set<anchored_binrule_type> best_anchored_binaries;
    std::unordered_set<anchored_unirule_type> best_anchored_unaries;
    std::unordered_set<anchored_lexrule_type> best_anchored_lexicals;


    if (not best_tree)
    {
      //std::cerr << "(())" << std::endl;
      sshyp << "(())" << std::endl;
    }
    else
    {
      //std::cerr << *best_tree << std::endl;
      sshyp << *best_tree << std::endl;

      best_tree->anchored_productions(best_anchored_binaries,
                                      best_anchored_unaries,
                                      best_anchored_lexicals);
    }


    //binary rules
    for (const auto& ref_anc_bin : anchored_binaries)
    {
      if (not best_anchored_binaries.count(ref_anc_bin))
      {
        local_corrects.push_back( network.rule_expression(std::get<3>(ref_anc_bin).get_lhs(),
                                                          std::get<3>(ref_anc_bin).get_rhs0(),
                                                          std::get<3>(ref_anc_bin).get_rhs1()
                                                          ));

        if((span_level > 0) and (std::get<1>(ref_anc_bin) - std::get<0>(ref_anc_bin) > 2))
          local_corrects.push_back( network.span_expression(std::get<3>(ref_anc_bin).get_lhs(),
                                                            std::get<0>(ref_anc_bin),
                                                            std::get<1>(ref_anc_bin) -1
                                                            ));
      }
    }

    for (const auto& best_anc_bin : best_anchored_binaries)
    {
      if (not anchored_binaries.count(best_anc_bin))
      {
        local_errs.push_back(network.rule_expression(std::get<3>(best_anc_bin).get_lhs(),
                                                     std::get<3>(best_anc_bin).get_rhs0(),
                                                     std::get<3>(best_anc_bin).get_rhs1()
                                                     ));

        if((span_level > 0) and (std::get<1>(best_anc_bin) - std::get<0>(best_anc_bin) > 2))
          local_errs.push_back( network.span_expression(std::get<3>(best_anc_bin).get_lhs(),
                                                        std::get<0>(best_anc_bin),
                                                        std::get<1>(best_anc_bin) - 1
                                                        ));
      }
    }


    //unary rules
    for (const auto& ref_anc_un : anchored_unaries)
    {
      if (not best_anchored_unaries.count(ref_anc_un))
      {
        local_corrects.push_back( network.rule_expression(std::get<2>(ref_anc_un).get_lhs(),
                                                          std::get<2>(ref_anc_un).get_rhs0(),
                                                          SymbolTable::instance_nt().get_symbol_count()
                                                          ));
      }
    }

    for (const auto& best_anc_un : best_anchored_unaries)
    {
      if (not anchored_unaries.count(best_anc_un))
      {
        local_errs.push_back(network.rule_expression(std::get<2>(best_anc_un).get_lhs(),
                                                     std::get<2>(best_anc_un).get_rhs0(),
                                                     SymbolTable::instance_nt().get_symbol_count()
                                                     ));
      }
    }


    //lexical rules
    for (const auto& ref_anc_lex : anchored_lexicals)
    {
      if (not best_anchored_lexicals.count(ref_anc_lex))
      {
        local_corrects.push_back( network.lexical_rule_expression(std::get<1>(ref_anc_lex).get_lhs(),
                                                                  std::get<0>(ref_anc_lex)
                                                                  ));
      }
    }

    for (const auto& best_anc_lex : best_anchored_lexicals)
    {
      if (not anchored_lexicals.count(best_anc_lex))
      {
        local_errs.push_back(network.lexical_rule_expression(std::get<1>(best_anc_lex).get_lhs(),
                                                             std::get<0>(best_anc_lex)
                                                             ));
      }
    }
  }
  if (best_tree) delete best_tree;

  return std::make_pair(
      std::make_pair(local_corrects,local_errs),
      std::make_pair(ssref.str(), sshyp.str()));
}



int NNLorgParseApp::run_train()
{
  unix_parse_solution::init();
  json_parse_solution::init();

  if(verbose) std::clog << "Start learning process.\n";

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

  ParserCKYNN parser(grammar);

#ifdef USE_THREADS
  parser.set_nbthreads(nbthreads);
#endif

  parser.set_word_signature(ws);
  Tagger tagger(&(parser.get_words_to_rules()));

  dynet::Model m;

  //dynet::SimpleSGDTrainer trainer(&m);
  // trainer.eta_decay = 0.08;
  //dynet::MomentumSGDTrainer trainer(&m);
  dynet::AdamTrainer trainer(&m);


  int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  // clock_t parse_start = (verbose) ? clock() : 0;


    std::vector<ParserCKYNN::scorer> networks;
    for (unsigned thidx = 0; thidx < nbthreads; ++ thidx)
      networks.push_back(ParserCKYNN::scorer(m, lstm_level, span_level));

  //  auto lhs_int_vec =std::vector<int>(grammar.lhs_int_set.begin(), grammar.lhs_int_set.end());

  for (unsigned iteration = 0; iteration < iterations; ++iteration)
  {
    if (verbose) std::cerr << "Iteration: " << iteration << std::endl;

    double loss = 0.0;

    std::stringstream oshyp;
    oshyp << train_output_name << "-train-hyp-" << iteration << ".mrg";
    std::ofstream outhyp(oshyp.str());

    std::random_shuffle(std::begin(trees), std::end(trees));

    std::stringstream osref;
    osref << train_output_name << "-train-ref-" << iteration << ".mrg";
    std::ofstream outref(osref.str());

    unsigned nb_chunks = trees.size() / batch_size;
    if (nb_chunks * batch_size < trees.size()) ++ nb_chunks;

    for (unsigned chunk = 0; chunk < nb_chunks; ++ chunk)
    {

      if (verbose) std::cerr << "Mini-batch: " << chunk << "/" << nb_chunks << std::endl;

      // computation graph for the mini batch
      dynet::ComputationGraph cg;
      // collect errors for the mini batch
      std::vector<dynet::expr::Expression> errs;

      for (unsigned thidx = 0; thidx < nbthreads; ++thidx)
      {
        networks[thidx].set_cg(cg);
        networks[thidx].clear();

        networks[thidx].precompute_rule_expressions(grammar.binary_rules, grammar.unary_rules);
      }
      auto lower_bound = chunk * batch_size;
      auto upper_bound = std::min<unsigned long>(trees.size(), (chunk+1)*batch_size);

      auto segment = (upper_bound - lower_bound) / nbthreads;


      std::vector<std::vector<dynet::expr::Expression>> thread_corrects(nbthreads);
      std::vector<std::vector<dynet::expr::Expression>> thread_errs(nbthreads);
      std::vector<std::vector<std::pair<std::string,std::string>>> thread_treestrings(nbthreads);
      std::vector<std::thread> threads;



        auto f =
            [&](unsigned i, unsigned mi, unsigned ma)
            {
              for(auto idx = mi; idx < ma; ++idx)
              {
                auto res = train_instance(trees[idx], parser,
                                          tagger, networks[i],
                                          start_symbol
                                          );
                if (verbose) std::cerr << '|' ;
                auto& local_corrects = res.first.first;
                auto& local_errs = res.first.second;
                auto& strpair = res.second;

                if (not local_corrects.empty())
                  thread_corrects[i].insert(thread_corrects[i].end(), local_corrects.begin(), local_corrects.end());
                if (not local_errs.empty())
                  thread_errs[i].insert(thread_errs[i].end(), local_errs.begin(), local_errs.end());


                thread_treestrings[i].push_back(strpair);
              }
            };



      for (unsigned thidx = 0; thidx < nbthreads; ++thidx)
      {
        auto mi = lower_bound + thidx * segment;
        auto ma = lower_bound + (thidx + 1) * segment;


        //f(thidx,mi,ma);
        threads.push_back(std::thread(f,thidx,mi,ma));
      }
      for (auto& thread : threads)
      {
        thread.join();
      }
      std::cerr << std::endl;


      for (unsigned thidx = 0; thidx < nbthreads; ++thidx)
      {
        if (not thread_corrects[thidx].empty())
          errs.push_back(- dynet::expr::sum(thread_corrects[thidx]));

        if (not thread_errs[thidx].empty())
          errs.push_back(dynet::expr::sum(thread_errs[thidx]));


        for (const auto& p : thread_treestrings[thidx])
        {
          outref << p.first;
          outhyp << p.second;
        }
      }

      if (not errs.empty())
      {
        //std::cerr << "here 1" << std::endl;
        dynet::expr::Expression s = dynet::expr::sum(errs);
        //std::cerr << "here 2" << std::endl;
        cg.incremental_forward(s);
        //std::cerr << "here 3" << std::endl;
        cg.backward(s.i);
        //std::cerr << "here 4" << std::endl;
        loss += as_scalar(cg.get_value(s.i));
        //std::cerr << "here 5" << std::endl;
        //std::cerr << loss << std::endl;
        //std::cerr << "here 6" << std::endl;
        trainer.update(1.0);
        //std::cerr << "here 7" << std::endl;
      }
    }

    std::cerr << "ending iteration " << loss << std::endl;
    trainer.status();

    trainer.update_epoch();

    std::stringstream mout;
    mout << train_output_name << "-model-" << iteration;
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
      devss << train_output_name << "-dev-hyp-" << iteration;
      std::ofstream devstream(devss.str());

      //rewind test-input
      delete in;
      in = new std::ifstream(in_filename.c_str());

      std::vector<std::string> comments;
      while(tokeniser->tokenise(*in,test_sentence,s,brackets, comments)) {
        clock_t sent_start = (verbose) ? clock() : 0;


        dynet::ComputationGraph cgdev;
        networks[0].set_cg(cgdev);
        networks[0].unset_gold();
        networks[0].clear();

        // the pointer to the solution
        PtbPsTree*  best_tree = nullptr;

        // check length of input sentence
        if(s.size() <= max_length) {

          // tag
          //std::cerr << "tag" << std::endl;
          tagger.tag(s, *(parser.get_word_signature()));
          networks[0].set_words(s);
          networks[0].precompute_embeddings();
          // should save values once and for all
          networks[0].precompute_rule_expressions(grammar.binary_rules, grammar.unary_rules);
          if (span_level > 0) networks[0].precompute_span_expressions();

          // create and initialise chart
          //std::cerr << "chart" << std::endl;
          ParserCKYNN::Chart chart(s,parser.get_nonterm_count(),brackets, networks[0]);

          // parse
          //std::cerr << "parse" << std::endl;
          parser.parse(chart, networks[0]);

          // get results
          //std::cerr << "results" << std::endl;
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
  options.add_nn_parser_options();

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

  // get training grammar
  if(configuration.exists("grammar")) {
    train_output_name = configuration.get_value<std::string>("grammar");
    if(verbose)
      std::cerr << "Setting grammar to " << train_output_name << ".\n";
  }
  else {
    //std::cerr << "grammar wasn't set." << std::endl;
    //return false;
  }
  cutoff = configuration.get_value<unsigned>("unknown-word-cutoff");

  ws = WordSignatureFactory::create_wordsignature(configuration);

  batch_size =  configuration.get_value<unsigned>("batch-size");
  iterations = configuration.get_value<unsigned>("iterations");
  lstm_level = configuration.get_value<unsigned>("lstm-level");
  span_level = configuration.get_value<unsigned>("span-level");

  return true;
}
