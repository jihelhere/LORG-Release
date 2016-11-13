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
    : LorgParseApp()
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

template<typename T>
void serialize_object(const T& obj, const std::string& filename)
{
  std::stringstream gout;
  std::ofstream gs(filename);
  boost::archive::text_oarchive goa(gs);
  goa << obj;
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





// Assume that rule scores are already precomputed
PtbPsTree *
NNLorgParseApp::parse_instance(const std::vector<Word>& words,
                               const ParserCKYNN& parser,
                               ParserCKYNN::scorer& network,
                               int start_symbol,
                               const std::vector<int>& lhs_int_vec)
{

  PtbPsTree*  best_tree = nullptr;

  //std::cerr << "set words" << std::endl;
  network.set_words(words);

  //std::cerr << "set embeddings" << std::endl;
  network.precompute_embeddings();

  if (span_level > 0)
  {
    network.span_expressions_bin.clear();
    network.span_expressions_un.clear();

    network.span_scores_bin.clear();
    network.span_scores_un.clear();

    network.precompute_span_expressions(lhs_int_vec);
  }

  network.lexical_expressions.clear();

  // create and initialise chart
  //std::cerr << "chart" << std::endl;
  std::vector<bracketing> brackets;

  ParserCKYNN::Chart chart(words,parser.get_nonterm_count(), brackets, network);

  //std::cerr << "parse" << std::endl;
  parser.parse(chart, network);

  best_tree = chart.get_best_tree(start_symbol, 0);

  return best_tree;
}

//////
// std::pair<
//   std::pair<std::vector<dynet::expr::Expression>,std::vector<dynet::expr::Expression>>,
//   std::pair<std::string,std::string>>
NNLorgParseApp::train_item
NNLorgParseApp::train_instance(const PtbPsTree& tree,
                               const ParserCKYNN& parser,
                               const Tagger& tagger,
                               ParserCKYNN::scorer& network,
                               int start_symbol,
                               const std::vector<int>& lhs_int_vec)
{
  std::vector<dynet::expr::Expression> local_corrects, local_errs;
  std::stringstream ssref,sshyp;

  std::unordered_set<anchored_binrule_type> anchored_binaries;
  std::unordered_set<anchored_unirule_type> anchored_unaries;
  std::unordered_set<anchored_lexrule_type> anchored_lexicals;

  tree.anchored_productions(anchored_binaries, anchored_unaries, anchored_lexicals);
  network.set_gold(anchored_binaries, anchored_unaries, anchored_lexicals);

  const auto s = tree.yield();
  PtbPsTree*  best_tree = nullptr;

  if (s.size() <= max_length)
  {
    //std::cerr << tree << std::endl;
    ssref << tree << std::endl;

    std::vector<Word> words;
    for (auto i = 0U; i < s.size(); ++i)
    {
      words.push_back(Word(s[i],i,i+1));
    }

    //std::cerr << "tag" << std::endl;
    tagger.tag(words, *(parser.get_word_signature()));

    best_tree = parse_instance(words,
                               parser,
                               network,
                               start_symbol,
                               lhs_int_vec);



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

    std::unordered_map<const dynet::expr::Expression*, int> exp_diff_count;

    //binary rules
    for (const auto& ref_anc_bin : anchored_binaries)
    {
      auto&& r = std::get<3>(ref_anc_bin);
      auto k = &nn_scorer::rule_expressions[nn_scorer::nt_triple_to_index(r.get_lhs(), r.get_rhs0(), r.get_rhs1())];
      exp_diff_count[k]--; // should be zero for non-existent k

      if (span_level > 0)
      {
        auto begin = std::get<0>(ref_anc_bin);
        auto end =   std::get<1>(ref_anc_bin) -1;
        auto split = std::get<2>(ref_anc_bin);

        k = &network.span_expression(r.get_lhs(), begin, end, split);
        exp_diff_count[k]--; // should be zero for non-existent k

        k = &network.span_init(r.get_lhs(), begin);
        exp_diff_count[k]--;

        k = &network.span_end(r.get_lhs(), end);
        exp_diff_count[k]--;

        // k = &network.span_split(r.get_lhs(), split);
        // exp_diff_count[k]--;
      }
    }

    for (const auto& best_anc_bin : best_anchored_binaries)
    {
      auto&& r = std::get<3>(best_anc_bin);
      auto k = &nn_scorer::rule_expressions[nn_scorer::nt_triple_to_index(r.get_lhs(), r.get_rhs0(), r.get_rhs1())];
      exp_diff_count[k]++; // should be zero for non-existent k

      if (span_level > 0)
      {
        auto begin = std::get<0>(best_anc_bin);
        auto end =   std::get<1>(best_anc_bin) -1;
        auto split = std::get<2>(best_anc_bin);

        k = &network.span_expression(r.get_lhs(), begin, end, split);
        exp_diff_count[k]++; // should be zero for non-existent k

        k = &network.span_init(r.get_lhs(), begin);
        exp_diff_count[k]++;

        k = &network.span_end(r.get_lhs(), end);
        exp_diff_count[k]++;

        // k = &network.span_split(r.get_lhs(), split);
        // exp_diff_count[k]++;
      }
    }


    //unary rules
    for (const auto& ref_anc_un : anchored_unaries)
    {
      auto&& r = std::get<2>(ref_anc_un);
      auto k = &nn_scorer::rule_expressions[nn_scorer::nt_triple_to_index(r.get_lhs(),r.get_rhs0(),-1)];
      exp_diff_count[k]--;

      if (span_level > 0)
      {
        auto begin = std::get<0>(ref_anc_un);
        auto end =   std::get<1>(ref_anc_un) -1;
        auto split = -1;

        k = &network.span_expression(r.get_lhs(), begin, end, split);
        exp_diff_count[k]--; // should be zero for non-existent k

        k = &network.span_init(r.get_lhs(), begin);
        exp_diff_count[k]--;

        k = &network.span_end(r.get_lhs(), end);
        exp_diff_count[k]--;
      }
    }

    for (const auto& best_anc_un : best_anchored_unaries)
    {
      auto&& r = std::get<2>(best_anc_un);
      auto k = &nn_scorer::rule_expressions[nn_scorer::nt_triple_to_index(r.get_lhs(),r.get_rhs0(), -1)];
      exp_diff_count[k]++;

      if (span_level > 0)
      {
        auto begin = std::get<0>(best_anc_un);
        auto end =   std::get<1>(best_anc_un) -1;
        auto split = -1;

        k = &network.span_expression(r.get_lhs(), begin, end, split);
        exp_diff_count[k]++; // should be zero for non-existent k

        k = &network.span_init(r.get_lhs(), begin);
        exp_diff_count[k]++;

        k = &network.span_end(r.get_lhs(), end);
        exp_diff_count[k]++;
      }
    }


    //lexical rules
    for (const auto& ref_anc_lex : anchored_lexicals)
    {
      auto k = &network.lexical_expressions[std::make_tuple(std::get<1>(ref_anc_lex).get_lhs(),
                                                            std::get<0>(ref_anc_lex))];
      exp_diff_count[k]--;
    }

    for (const auto& best_anc_lex : best_anchored_lexicals)
    {
      auto k = &network.lexical_expressions[std::make_tuple(std::get<1>(best_anc_lex).get_lhs(),
                                                                std::get<0>(best_anc_lex))];
      exp_diff_count[k]++;
    }

    // todo : why not return the table directly ??
    for(auto& kv : exp_diff_count)
    {
      if (kv.second > 0)
      {
        for (auto i = 0; i < kv.second; ++i)

          local_errs.emplace_back(*(kv.first));
      }
      else
      {
        for (auto i = 0; i < -kv.second; ++i)
          local_corrects.emplace_back(*(kv.first));
      }
    }
  }
  if (best_tree) delete best_tree;

  return {local_corrects,local_errs, ssref.str(), sshyp.str()};
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


  collect_rules(trees, bin, un, lex);


  Grammar<Rule,Rule,Rule> grammar;
  grammar.set_rules(std::vector<Rule>(bin.begin(),bin.end()),
                    std::vector<Rule>(un.begin(),un.end()),
                    std::vector<Rule>(lex.begin(),lex.end())
                    );

  write_symboltable(SymbolTable::instance_nt(), "NTsymbols.txt");
  //write_symboltable(SymbolTable::instance_word(), "Tsymbols.txt");
  std::stringstream wout;
  wout << train_output_name << ".ts";
  serialize_object(SymbolTable::instance_word(), wout.str());

  std::stringstream nout;
  nout << train_output_name << ".nts";
  serialize_object(SymbolTable::instance_nt(), nout.str());

  std::stringstream gout;
  gout << train_output_name << ".grammar";
  serialize_object(grammar, gout.str());

  ParserCKYNN parser(grammar);

#ifdef USE_THREADS
  parser.set_nbthreads(nbthreads);
#endif

  parser.set_word_signature(ws);
  Tagger tagger(&(parser.get_words_to_rules()));

  dynet::Model m;

  //dynet::SimpleSGDTrainer trainer(&m);
  //dynet::MomentumSGDTrainer trainer(&m);
  dynet::AdamTrainer trainer(&m);
  trainer.eta_decay = 0.08;


  int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  // clock_t parse_start = (verbose) ? clock() : 0;


  std::vector<ParserCKYNN::scorer> networks;
  for (unsigned thidx = 0; thidx < nbthreads; ++ thidx)
    networks.emplace_back(m, lstm_level, span_level,
                          word_embedding_size,
                          nt_embedding_size,
                          hidden_size,
                          lstm_hidden_size,
                          use_char_embeddings,
                          use_span_midpoints
                          );

  auto&& lhs_int_vec =std::vector<int>(grammar.lhs_int_set.begin(), grammar.lhs_int_set.end());

  for (unsigned iteration = 0; iteration < iterations; ++iteration)
  {
    if (verbose) std::cerr << "Iteration: " << iteration << std::endl;

    double loss = 0.0;

    std::stringstream oshyp;
    oshyp << train_output_name
          << "-b" << batch_size
          << "-s" << span_level
          << "-l" << lstm_level
          << "-w" << word_embedding_size
          << "-n" << nt_embedding_size
          << "-h" << hidden_size
          << "-lh" << lstm_hidden_size
          << "-d" << dropout
          << "-train-hyp-" << iteration << ".mrg";
    std::ofstream outhyp(oshyp.str());

    std::random_shuffle(std::begin(trees), std::end(trees));

    std::stringstream osref;
    osref << train_output_name
          << "-b" << batch_size
          << "-s" << span_level
          << "-l" << lstm_level
          << "-w" << word_embedding_size
          << "-n" << nt_embedding_size
          << "-h" << hidden_size
          << "-lh" << lstm_hidden_size
          << "-d" << dropout
          << "-train-ref-" << iteration << ".mrg";
    std::ofstream outref(osref.str());

    unsigned nb_chunks = trees.size() / batch_size;
    if (nb_chunks * batch_size < trees.size()) ++ nb_chunks;

    for (unsigned chunk = 0; chunk < nb_chunks; ++ chunk)
    {
      if (verbose) std::cerr << "Mini-batch: " << chunk << "/" << nb_chunks << std::endl;

      dynet::ComputationGraph cg;

      for (auto& network: networks)
      {
        network.clear();
        network.set_dropout(dropout);
      }

      nn_scorer::train_mode = true;
      nn_scorer::set_cg(cg);
      nn_scorer::precompute_rule_expressions(grammar.binary_rules, grammar.unary_rules);

      // collect errors for the mini batch
      std::vector<dynet::expr::Expression> errs;

      auto lower_bound = chunk * batch_size;
      auto upper_bound = std::min<unsigned long>(trees.size(), (chunk+1)*batch_size);

      auto segment = (upper_bound - lower_bound) / nbthreads;
      if ((lower_bound + nbthreads * segment) < upper_bound)
        ++ segment;


      std::vector<std::vector<dynet::expr::Expression>> thread_corrects(nbthreads);
      std::vector<std::vector<dynet::expr::Expression>> thread_errs(nbthreads);
      std::vector<std::vector<std::pair<std::string,std::string>>> thread_treestrings(nbthreads);
      std::vector<std::thread> threads;



        auto f =
            [&](unsigned i, unsigned mi, unsigned ma)
            {
              if (ma > upper_bound)
                ma = upper_bound;
              for(auto idx = mi; idx < ma; ++idx)
              {
                auto res = train_instance(trees[idx], parser,
                                          tagger, networks[i],
                                          start_symbol,
                                          lhs_int_vec
                                          );
                if (verbose) std::cerr << '|' ;
                auto&& local_corrects = res.correct_exprs;
                auto&& local_errs = res.error_exprs;

                if (not local_corrects.empty())
                  thread_corrects[i].insert(thread_corrects[i].end(), local_corrects.begin(), local_corrects.end());
                if (not local_errs.empty())
                  thread_errs[i].insert(thread_errs[i].end(), local_errs.begin(), local_errs.end());


                thread_treestrings[i].emplace_back(res.ref_tree_string, res.hyp_tree_string);
              }
            };



      for (unsigned thidx = 0; thidx < nbthreads; ++thidx)
      {
        auto mi = lower_bound + thidx * segment;
        auto ma = lower_bound + (thidx + 1) * segment;


        //f(thidx,mi,ma);
        threads.emplace_back(std::thread(f,thidx,mi,ma));
      }

      for (auto& thread : threads)
      {
        thread.join();
      }
      std::cerr << std::endl;


      for (unsigned thidx = 0; thidx < nbthreads; ++thidx)
      {
        if (not thread_corrects[thidx].empty())
          errs.emplace_back(- dynet::expr::sum(thread_corrects[thidx]));

        if (not thread_errs[thidx].empty())
          errs.emplace_back(dynet::expr::sum(thread_errs[thidx]));
        loss += thread_errs[thidx].size();

        for (const auto& p : thread_treestrings[thidx])
        {
          outref << p.first;
          outhyp << p.second;
        }
      }

      if (not errs.empty())
      {
        dynet::expr::Expression s = dynet::expr::sum(errs);
        cg.incremental_forward(s);
        cg.backward(s.i);
        loss += as_scalar(cg.get_value(s.i));
        trainer.update(1.0);
      }
    }

    std::cerr << "ending iteration " << loss << std::endl;
    trainer.status();

    trainer.update_epoch();

    std::stringstream mout;
    mout << train_output_name
         << "-b" << batch_size
         << "-s" << span_level
         << "-l" << lstm_level
         << "-w" << word_embedding_size
         << "-n" << nt_embedding_size
         << "-h" << hidden_size
         << "-lh" << lstm_hidden_size
         << "-d" << dropout
         << "-model-" << iteration;
    serialize_object(m, mout.str());

    ///////////// process dev file
    {
      int count = 0;

      int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

      //clock_t parse_start = (verbose) ? clock() : 0;

      std::stringstream devss;
      devss << train_output_name
            << "-b" << batch_size
            << "-s" << span_level
            << "-l" << lstm_level
            << "-w" << word_embedding_size
            << "-n" << nt_embedding_size
            << "-h" << hidden_size
            << "-lh" << lstm_hidden_size
            << "-d" << dropout
            << "-dev-hyp-" << iteration;
      std::ofstream devstream(devss.str());

      //rewind test-input
      delete in;
      in = new std::ifstream(in_filename.c_str());

      auto&& lhs_int_vec =std::vector<int>(grammar.lhs_int_set.begin(), grammar.lhs_int_set.end());


      for (auto& network : networks)
      {
        network.unset_dropout();
      }


      nn_scorer::train_mode = false;

      std::vector<std::vector<Word>> s(nbthreads);
      std::vector<std::vector< bracketing>> brackets(nbthreads);
      std::vector<std::vector<std::pair<PtbPsTree*,double>>> solutions(nbthreads, {{nullptr,1.0}});

      bool allvalid = true;
      while(allvalid)
      {
        dynet::ComputationGraph cgdev;

        std::vector<std::string> test_sentence(nbthreads);
        std::vector<std::vector<std::string>> comments(nbthreads);

        std::vector<bool> valid_sentence(nbthreads, true);

        std::vector<clock_t> clock_start(nbthreads,0);
        std::vector<clock_t> clock_end(nbthreads,0);

        std::vector<std::thread> threads;

        nn_scorer::train_mode = false;
        nn_scorer::set_cg(cgdev);
        nn_scorer::precompute_rule_expressions(grammar.binary_rules, grammar.unary_rules);

        for (unsigned i = 0; i < nbthreads; ++i)
        {
          networks[i].unset_gold();

          valid_sentence[i] = tokeniser->tokenise(*in,test_sentence[i],s[i],brackets[i], comments[i]);
          allvalid = allvalid and valid_sentence[i];

          auto f = [&](unsigned idx)
                   {
                     if (not valid_sentence[idx]) return;
                     clock_start[idx] = (verbose) ? clock() : 0;

                     // check length of input sentence
                     if(s[idx].size() <= max_length)
                     {
                       tagger.tag(s[idx], *(parser.get_word_signature()));
                       solutions[idx][0].first = parse_instance(s[idx], parser, networks[idx], start_symbol, lhs_int_vec);
                     }
                     clock_end[idx] = (verbose) ? clock() : 0;
                   };
          //f(i);
          threads.emplace_back(std::thread(f,i));
        }

        for (auto& thread : threads)
          thread.join();

        for (unsigned i = 0; i < nbthreads; ++i)
        {
          if (valid_sentence[i])
          {
            auto* p_typed =
                parse_solution::factory.create_object(parse_solution::UNIX,
                                                      parse_solution(test_sentence[i], ++count,
                                                                     s[i].size(), solutions[i],
                                                                     (verbose) ? clock_end[i] -clock_start[i] / double(CLOCKS_PER_SEC) : 0,
                                                                     verbose, comments[i],false)
                                                      );
            p_typed->print(devstream);
            delete p_typed;
          }


          delete solutions[i][0].first;
          solutions[i][0].first = nullptr;
          s[i].clear();
          brackets[i].clear();

        }
      }
    }
  }

  return 0;

}


  int NNLorgParseApp::run()
  {
    if (train_mode) return run_train();

    if(verbose) std::clog << "Start parsing process.\n";


    SymbolTable::instance_nt().load(train_output_name + ".nts");
    SymbolTable::instance_word().load(train_output_name + ".ts");

    Grammar<Rule,Rule,Rule> grammar;
    std::stringstream gin;
    std::ifstream gs(train_output_name + ".grammar");
    boost::archive::text_iarchive gia(gs);
    gia >> grammar;


    ParserCKYNN parser(grammar);

#ifdef USE_THREADS
    parser.set_nbthreads(nbthreads);
#endif

    parser.set_word_signature(ws);
    Tagger tagger(&(parser.get_words_to_rules()));

    dynet::Model m;
    std::ifstream ms(test_model_name);
    boost::archive::text_iarchive mia(ms);
    mia >> m;

    int count = 0;
    int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

    // clock_t parse_start = (verbose) ? clock() : 0;


    std::vector<ParserCKYNN::scorer> networks;
    for (unsigned thidx = 0; thidx < nbthreads; ++ thidx)
      networks.emplace_back(m, lstm_level, span_level,
                            word_embedding_size,
                            nt_embedding_size,
                            hidden_size,
                            lstm_hidden_size,
                            use_char_embeddings,
                            use_span_midpoints
                            );

    auto&& lhs_int_vec =std::vector<int>(grammar.lhs_int_set.begin(), grammar.lhs_int_set.end());

    // this is moved here and si sequential
    // bc a mutex is needed
    {
      dynet::ComputationGraph cg;
      networks[0].set_cg(cg);
      networks[0].precompute_rule_expressions(grammar.binary_rules, grammar.unary_rules);
      networks[0].unset_dropout();
    }
    for (unsigned i = 1; i < nbthreads; ++i)
    {
      networks[i].rule_scores = networks[0].rule_scores;
      networks[i].unset_dropout();
    }

    std::vector<std::vector<Word>> s(nbthreads);
    std::vector<std::vector< bracketing>> brackets(nbthreads);
    std::vector<std::vector<std::pair<PtbPsTree*,double>>> solutions(nbthreads, {{nullptr,1.0}});

    bool allvalid = true;
    while(allvalid)
    {
      dynet::ComputationGraph cgdev;

      std::vector<std::string> test_sentence(nbthreads);
      std::vector<std::vector<std::string>> comments(nbthreads);

      std::vector<bool> valid_sentence(nbthreads, true);

      std::vector<clock_t> clock_start(nbthreads,0);
      std::vector<clock_t> clock_end(nbthreads,0);

      std::vector<std::thread> threads;

      for (unsigned i = 0; i < nbthreads; ++i)
      {
        networks[i].set_cg(cgdev);
        networks[i].unset_gold();

        valid_sentence[i] = tokeniser->tokenise(*in,test_sentence[i],s[i],brackets[i], comments[i]);
        allvalid = allvalid and valid_sentence[i];

        auto f = [&](unsigned idx)
                 {
                   if (not valid_sentence[idx]) return;
                   clock_start[idx] = (verbose) ? clock() : 0;

                   // check length of input sentence
                   if(s[idx].size() <= max_length)
                   {
                     tagger.tag(s[idx], *(parser.get_word_signature()));
                     solutions[idx][0].first = parse_instance(s[idx], parser, networks[idx], start_symbol, lhs_int_vec);
                   }
                   clock_end[idx] = (verbose) ? clock() : 0;
                 };
        //f(i);
        threads.emplace_back(std::thread(f,i));
      }

      for (auto& thread : threads)
        thread.join();

      for (unsigned i = 0; i < nbthreads; ++i)
      {
        if (valid_sentence[i])
        {
          auto* p_typed =
              parse_solution::factory.create_object(parse_solution::UNIX,
                                                    parse_solution(test_sentence[i], ++count,
                                                                   s[i].size(), solutions[i],
                                                                   (verbose) ? clock_end[i] -clock_start[i] / double(CLOCKS_PER_SEC) : 0,
                                                                   verbose, comments[i],false)
                                                    );
          p_typed->print(*out);
          delete p_typed;
        }


        delete solutions[i][0].first;
        solutions[i][0].first = nullptr;
        s[i].clear();
        brackets[i].clear();

      }
    }

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

  if (configuration.exists("train")) {
    train_mode = configuration.get_value<bool>("train");
  }

  if (train_mode)
  {

    if (configuration.exists("treebank")) {
      tb_options = TreebankFactory::read_config(configuration);
    }
    else
    {
      std::cerr << "treebank was not set. Exiting...\n";
      return false;
    }
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

  word_embedding_size = configuration.get_value<unsigned>("word-embedding-size");
  nt_embedding_size = configuration.get_value<unsigned>("nt-embedding-size");
  hidden_size = configuration.get_value<unsigned>("hidden-size");
  lstm_hidden_size = configuration.get_value<unsigned>("lstm-hidden-size");
  dropout = configuration.get_value<float>("dropout");
  use_char_embeddings = configuration.get_value<bool>("use-char-embeddings");
  use_span_midpoints = configuration.get_value<bool>("use-span-midpoints");

  test_model_name = configuration.get_value<std::string>("test-model-name");

  return true;
}
