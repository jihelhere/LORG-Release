// -*- mode: c++ -*-
#include "TwoStageLorgParseApp.h"

#include "utils/Tagger.h"
#include "feature_extract/Extract.h"
#include "utils/PtbPsTree.h"
#include "utils/LorgConstants.h"
#include "parsers/ParserCKYAllFactory.h"

#include "parsers/ParserCKYAll.h"

#include "utils/tick_count.h"

#include "lexicon/WordSignatureFactory.h"

#include <thread>



TwoStageLorgParseApp::TwoStageLorgParseApp() : LorgParseApp(), parsers(1)
{
  unix_parse_solution::init();
  json_parse_solution::init();
}

TwoStageLorgParseApp::~TwoStageLorgParseApp()
{
  for(auto& p : this->parsers)
    if(p) delete p;
}

int TwoStageLorgParseApp::run()
{
  parse_solution::init_feature_extractor();

  if(verbose)  std::clog << "Start parsing process.\n";

  std::string raw_sentence; // the sentence string (only for pretty print solutions)
  std::vector<Word> sentence; // the vector of words used by the parser to initialise the chart
  std::vector< bracketing > brackets; // pre-bracketting for guided parsing
  std::vector<std::string> comments; // strings to store comments associated with a sentence
  int count = 0; // sentence count (debug & pretty-print)

  int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name); // axiom of the grammar
  tick_count parse_start = tick_count::now();

  //read input and fill raw_sentence, sentence and brackets
  while(tokeniser->tokenise(*in, raw_sentence, sentence, brackets, comments)) {

    // //Extra verbose
    // if(verbose) {
    //   std::clog << "Tokens: ";
    //   for(std::vector<Word>::const_iterator i(sentence.begin()); i != sentence.end(); ++i)
    //   {
    //     std::clog << "<" << i->get_form() << ">";
    //     std::clog << "(" << i->get_start() << ", " << i->get_end() << ")";
    //   }
    //   std::clog << "\n";

    //   std::clog << "brackets " << brackets.size() << std::endl;
    // }

    std::vector<std::vector<Word>> sentences(parsers.size(), sentence);


    std::function<void(int)> process_sentence =
        [&](int i){

      //      std::cout << i << std::endl;

      //std::cerr << "tag" << std::endl;
      taggers[i].tag(sentences[i], *(parsers[i]->get_word_signature()));
      //std::cerr << "init chart" << std::endl;
      parsers[i]->initialise_chart(sentences[i], brackets);
      //std::cerr << "parse" << std::endl;
      parsers[i]->parse(start_symbol);
      //std::cerr << "beam" << std::endl;
      parsers[i]->beam_c2f(start_symbol);

      if(parsers[i]->is_chart_valid(start_symbol))
      {
        //std::cerr << "extract" << std::endl;
        parsers[i]->extract_solution();
      }

    };

    std::vector<std::thread> threads;

    tick_count sent_start = tick_count::now();

    if(sentence.size() <=  max_length && sentence.size() > 0) {

      for(size_t i = 0; i < parsers.size(); ++i)
      {
        //        threads.push_back(std::thread(process_sentence,i));
        // std::cerr << i << std::endl;
        process_sentence(i); // for replicability
        // TODO : find why solutions are different in multithreaded code
        // std::cerr << i << std::endl;
      }

      // for(auto& thread : threads)
      // {
      //   thread.join();
      // }


      int k = 0;

      std::vector<std::pair<PtbPsTree *,double> > best_trees; // vector of (tree,score)

      if (parsers.size() > 1)
        k = find_consensus(best_trees);

      if(verbose)
        std::cerr << "k: " << k << std::endl;

      for (size_t i = 0; i < 1; ++i)
        //for (size_t i = 0; i < parsers.size(); ++i)
      {
        if(parsers[i]->is_chart_valid(start_symbol))
        {
          //                             BLOCKTIMING("get_parses");
          parsers[i]->get_parses(start_symbol, kbest, always_output_forms, output_annotations, best_trees);
          //std::cout << "getting " << i << std::endl;
        }
        parse_solution * p_typed =
            parse_solution::factory.create_object(output_format,
                                                  parse_solution(raw_sentence, ++count,
                                                                 sentence.size(), best_trees,
                                                                 (verbose) ? (tick_count::now() - sent_start).seconds() : 0,
                                                                 verbose, comments, extract_features)
                                                  );
        p_typed->print(*out);
        delete p_typed;

      }
      //sanity
      for(unsigned l = 0; l < best_trees.size(); ++l) { // delete solutions
        delete best_trees[l].first;
      }
      best_trees.clear();

      for (size_t i = 0; i < parsers.size(); ++i)
      {
        parsers[i]->clean();
      }


      ///*
      if(verbose && count % 50 == 0)
        std::clog << count << " parsed sentences in " << (tick_count::now() - parse_start).seconds() << " sec" << std::endl;
      //*/


    }
    sentence.clear();
    brackets.clear();
    comments.clear();
  }

  *out << std::flush;

  if(verbose) std::clog << "overall time: " << (tick_count::now() - parse_start).seconds() << "s" << std::endl;
  return 0; //everything's fine
}

LorgOptions TwoStageLorgParseApp::get_options() const
{
  LorgOptions options(LorgParseApp::get_options());
  options.add_2sparser_options();
  return options;
}



bool TwoStageLorgParseApp::read_config(ConfigTable& configuration)
{
  if(LorgParseApp::read_config(configuration) == false) return false;

  output_annotations = configuration.get_value<bool>("output-annotations");

  if(verbose) { std::clog << "creating the parser... ";}

  this->parsers = ParserCKYAllFactory::create_parser(configuration);
  if(this->parsers.empty()) return false;

#ifdef USE_THREADS
  for (const auto& p: this->parsers)
  {
    p->set_nbthreads(this->nbthreads);
  }
#endif

  if(verbose) {std::clog << "ok" << std::endl;}

  kbest = configuration.get_value<unsigned>("kbest");

  //creating taggers
  taggers.resize(parsers.size());
  for (size_t i = 0; i < taggers.size(); ++i)
  {
    taggers[i].set_word_rules(&(parsers[i]->get_words_to_rules()));
  }


  for(const auto& p : this->parsers)
  {
    if (p->get_word_signature() == nullptr)
    {
      // read parts of the configuration dedicated to word signature
      // and initialize word class
      std::cerr << "here" << std::endl;
      p->set_word_signature(WordSignatureFactory::create_wordsignature(configuration));
    }
  }

  extract_features = configuration.get_value<bool>("extract-features");

  output_format = parse_solution::format_from_string(configuration.get_value<std::string>("output-format"));

  return true;
}

//////////////////////////// DD /////////////////////////////

unsigned simplified_nt( unsigned id )
{
  std::string name = SymbolTable::instance_nt().get_label_string(id);

  static const boost::regex exp_artificial ("^\\[\\((.*)\\)>\\]$");
  static const boost::regex exp_funct ("^([A-Za-z]+\\$?)[=-].*");
  boost::cmatch matched;
  bool is_artificial = false;

  if(boost::regex_match(name.c_str(), matched, exp_artificial))
  {
    name = std::string(matched[1].first, matched[1].second);
    is_artificial = true;
  }

  if(boost::regex_match(name.c_str(),matched,exp_funct))
  {
    name = std::string(matched[1].first, matched[1].second);
  }

  if(is_artificial)
  {
    name = "[(" + name + ")>]";
  }

  return SymbolTable::instance_nt().get_label_id(name);
}



int TwoStageLorgParseApp::find_consensus(std::vector<std::pair<PtbPsTree *,double> >& best_trees)
{

  struct anchored_symbol
  {
    int symbol;
    int begin;
    int end;

    anchored_symbol(int n, int b, int e) : symbol(n), begin(b), end(e) {};

    bool operator<(const anchored_symbol& other) const
    {
      return
          (begin < other.begin)
          or
          (begin == other.begin and end < other.end)
          or
          (begin == other.begin and end == other.end  and symbol < other.symbol);
    }
  };


  int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name); // axiom of the grammar


  double lu = 0;
  bool valid = true;

  for (size_t i = 0; i < parsers.size(); ++i)
  {
    if(not parsers[i]->is_chart_valid(start_symbol))
      valid = false;
    else
      lu += parsers[i]->get_best_score(start_symbol);
  }

  if(!valid) return -1;

  unsigned k = 0;
  double c = 1;
  double t = 0;

  size_t fun_parser_size = 0;
  for(auto& p : parsers)
  {
    if(p->get_is_funct()) ++ fun_parser_size;
  }

  for (k = 0; k < 1000; ++k)
  {
    //std::cout << "k = " << k << std::endl;
    // std::cout << "lu = " << lu << std::endl;
    // std::cout << "t = " << t << std::endl;
    double delta_k = c / (t + 1);
    double nlu = 0;

    std::vector<std::set< anchored_symbol >> z_parser(this->parsers.size());
    std::map<anchored_symbol, int> z_average;
    std::vector<MAP<int,MAP<int,MAP<int,double>>>> z_lambdas(parsers.size());


    // fun
    std::vector<std::set< anchored_symbol >> f_parser(this->parsers.size());
    std::map<anchored_symbol, int> f_average;
    std::vector<MAP<int,MAP<int,MAP<int,double>>>> f_lambdas(parsers.size());


    for (size_t i = 0; i < parsers.size(); ++i)
    {
      for (auto& e: parsers[i]->get_vectorized_representation(start_symbol))
      {
        //std::cout << std::get<1>(e) << " " << std::get<2>(e) << ": ";

        const AnnotatedRule* r = std::get<0>(e);

        int l;
        l  = simplified_nt(r->get_lhs());

        if(SymbolTable::instance_nt().get_label_string(l)[0] != '[')
          // not an 'artificial node'
        {
          anchored_symbol vv(l, std::get<1>(e), std::get<2>(e));
          if(z_parser[i].insert(vv).second)
          {
            z_average[vv] += 1;
          }
        }

        if (parsers[i]->get_is_funct() && (SymbolTable::instance_nt().get_label_string(l)[0] != '['))
        {
          int f = r->get_lhs();

          anchored_symbol ff(f, std::get<1>(e), std::get<2>(e));
          if(f_parser[i].insert(ff).second)
          {
            f_average[ff] += 1;
          }
        }
      }
    }


    // check same solution
    bool same = true;
    for (auto& p : z_average)
    {
      if (p.second != int(parsers.size()))
      {
        //std::cout << p.second << std::endl;
        same = false;
        break;
      }
    }
    for (auto& p : f_average)
    {
      if (p.second != int(fun_parser_size))
      {
        //std::cout << p.second << std::endl;
        same = false;
        break;
      }
    }

    if(same)
    {
      break;
    }

    for(size_t i = 0; i < z_parser.size(); ++i)
    {
      z_lambdas[i].clear();

      // add update for things in ith solution
      for(const auto& v: z_parser[i])
      {
        z_lambdas[i][v.begin][v.end][v.symbol]
            = delta_k * (1.0 - double(z_average[v]) / parsers.size());
      }

      // add update for things missing in ith solution
      for(const auto& v: z_average)
      {
        if(not z_lambdas[i].count(v.first.begin)
           or not z_lambdas[i][v.first.begin].count(v.first.end)
           or not z_lambdas[i][v.first.begin][v.first.end].count(v.first.symbol)
           )
        {
          z_lambdas[i][v.first.begin][v.first.end][v.first.symbol]
              = delta_k * (- double(z_average[v.first]) / parsers.size());
        }
      }
    }


    for(size_t i = 0; i < f_parser.size(); ++i)
    {
      if(not parsers[i]->get_is_funct()) continue;

      f_lambdas[i].clear();

      // add update for things in ith solution
      for(const auto& v: f_parser[i])
      {
        f_lambdas[i][v.begin][v.end][v.symbol]
            = delta_k * (1.0 - double(f_average[v]) / fun_parser_size);
      }

      // add update for things missing in ith solution
      for(const auto& v: f_average)
      {
        if(not f_lambdas[i].count(v.first.begin)
           or not f_lambdas[i][v.first.begin].count(v.first.end)
           or not f_lambdas[i][v.first.begin][v.first.end].count(v.first.symbol)
           )
        {
          f_lambdas[i][v.first.begin][v.first.end][v.first.symbol]
              = delta_k * (- double(f_average[v.first]) / fun_parser_size);
        }
      }
    }



    // update relaxations
    std::vector<std::thread> threads;
    for (size_t i = 0; i < this->parsers.size(); ++i)
    {
      threads.push_back(
          std::thread([&](int j)
                      {
                        if(this->parsers[j]->is_chart_valid(start_symbol))
                        {
                          this->parsers[j]->update_relaxations(true,  z_lambdas[j]);
                          if(parsers[j]->get_is_funct())
                            this->parsers[j]->update_relaxations(false, f_lambdas[j]);
                        }

                        },i));
    }
    for(auto& thread : threads)
    {
      thread.join();
    }


    //    std::cout << "HERE" << std::endl;

    // update solutions
    threads.clear();
    for (size_t i = 0; i < parsers.size(); ++i)
    {
      threads.push_back(
          std::thread([&](int j)
                      {if(parsers[j]->is_chart_valid(start_symbol))
                        {
                          parsers[j]->simple_extract_solution();
                        }
                      },i)
                        );
    }

    for(auto& thread : threads)
    {
      thread.join();
    }


    for (size_t i = 0; i < parsers.size(); ++i)
    {
      if(parsers[i]->is_chart_valid(start_symbol))
      {
        nlu += parsers[i]->get_best_score(start_symbol);
        // if (i == 0)
        // {
        //   parsers[i]->get_parses(start_symbol, 1, always_output_forms, output_annotations, best_trees);
        // }
      }
    }
    //        std::cout << std::endl;


    //update stepsize
    // std::cout << "nlu = " << nlu << std::endl;
    // std::cout << "lu = " << lu << std::endl;
    if (nlu > lu )//|| (k % 20) == 0)
    {
      ++t;
    }
    lu = nlu;
    //          std::cerr << nlu << std::endl;
  }

  // if (k == 1000)
  // {
  //   std::vector<std::pair<PtbPsTree *,double> > trees;

  //   for (const auto& p : parsers)
  //   {
  //     p->get_parses(start_symbol, 1, always_output_forms, output_annotations, trees);
  //   }
  //   for(const auto& t : trees)
  //   {
  //     t.first->unbinarise();
  //     std::cout << *(t.first) << std::endl;
  //     delete t.first;
  //   }

  // }

  return k;
}
