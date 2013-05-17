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
    //     std::clog << "<" << i->get_form() << ">";
    //   std::clog << "\n";
    // }

    std::vector<std::vector<Word>> sentences(parsers.size(), sentence);


    std::function<void(int)> process_sentence =
        [&](int i){

      //      std::cout << i << std::endl;

      taggers[i].tag(sentences[i], *(parsers[i]->get_word_signature()));
      parsers[i]->initialise_chart(sentences[i], brackets);
      parsers[i]->parse(start_symbol);
      parsers[i]->beam_c2f(start_symbol);
      if(parsers[i]->is_chart_valid(start_symbol))
      {
        parsers[i]->extract_solution();
      }

    };

    std::vector<std::thread> threads;

    tick_count sent_start = tick_count::now();

    if(sentence.size() <=  max_length && sentence.size() > 0) {

      for(size_t i = 0; i < parsers.size(); ++i)
      {
        threads.push_back(std::thread(process_sentence,i));
      }

      for(auto& thread : threads)
      {
        thread.join();
      }


      int k = 0;
      if (parsers.size() > 1)
        k = find_consensus();

      std::cerr << "k: " << k << std::endl;

      for (size_t i = 0; i < 1; ++i)
        //for (size_t i = 0; i < parsers.size(); ++i)
      {
        std::vector<std::pair<PtbPsTree *,double> > best_trees; // vector of (tree,score)
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


        //sanity
        for(unsigned l = 0; l < best_trees.size(); ++l) { // delete solutions
          delete best_trees[l].first;
        }
        best_trees.clear();
      }
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





int TwoStageLorgParseApp::find_consensus()
{

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
  MAP<int,MAP<int,MAP<int, MAP<int, MAP<int,double>>>>> u;
  double c = 1;
  double t = 0;

  for (k = 0; k < 1000; ++k)
  {
    // std::cout << "k = " << k << std::endl;
    // std::cout << "lu = " << lu << std::endl;
    // std::cout << "t = " << t << std::endl;
    double delta_k = c / (t + 1);
    double nlu = 0;

    std::vector<SET<std::tuple<const AnnotatedRule*,int,int> > > sets(this->parsers.size());
    std::vector<std::vector< std::vector<int> >> sets_v(this->parsers.size());

    for (size_t i = 0; i < parsers.size(); ++i)
    {
      sets[i] = parsers[i]->get_vectorized_representation(start_symbol);
      // if (sets[i].empty())
      //   break;

      for (auto& e: sets[i])
      {
        //std::cout << std::get<1>(e) << " " << std::get<2>(e) << ": ";

        auto& case_u = u[std::get<1>(e)][std::get<2>(e)];

        const AnnotatedRule* r = std::get<0>(e);

        int l, r0, r1;

        if(r->is_binary())
        {
          const BRule * br = static_cast<const BRule*>(r);
          l  = simplified_nt(br->get_lhs());
          r0 = simplified_nt(br->get_rhs0());
          r1 = simplified_nt(br->get_rhs1());

          sets_v[i].push_back({std::get<1>(e), std::get<2>(e), l,r0,r1});
        }
        else
        {
          if(r->is_lexical())
          {
            const LexicalRule * lr = static_cast<const LexicalRule*>(r);
            l  = simplified_nt( lr->get_lhs());
            //r0 = lr->get_rhs0();
            r0 = -1; // for lexical rules
            r1 = -2; // for lexical rules
            sets_v[i].push_back({std::get<1>(e), std::get<2>(e), l,r0,r1});
          }
          else //unary
          {
            const URule * ur = static_cast<const URule*>(r);
            l  = simplified_nt( ur->get_lhs());
            r0 = simplified_nt( ur->get_rhs0());
            r1 = -1; // for unary rules
            sets_v[i].push_back({std::get<1>(e), std::get<2>(e), l,r0,r1});
          }
        }

        // update u
        if ((i % 2) == 0)
        {
          case_u[l][r0][r1] -= delta_k;
        }
        else
        {
          case_u[l][r0][r1] += delta_k;
        }

      }
    }


    // check same solution
    bool same = false;
    if (sets_v[0].size() == sets_v[1].size())
    {
      same = true;
      for(const auto&s : sets_v[0])
      {
        if (std::find(sets_v[1].begin(), sets_v[1].end(),s) == sets_v[1].end())
        {
          same = false;
          break;
        }
      }
    }

    if(same)
    {
      break;
    }


    // update relaxations
    std::vector<std::thread> threads;
    for (size_t i = 0; i < this->parsers.size(); ++i)
    {
      threads.push_back(
          std::thread([&](int i)
                      {
                        if(this->parsers[i]->is_chart_valid(start_symbol))
                        {
                          this->parsers[i]->update_relaxations(u, (i % 2) == 0);
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
          std::thread([&](int i)
                      {if(parsers[i]->is_chart_valid(start_symbol))
                        {
                          parsers[i]->simple_extract_solution();
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
      }
    }
    //        std::cout << std::endl;


    //update stepsize
    // std::cout << "nlu = " << nlu << std::endl;
    // std::cout << "lu = " << lu << std::endl;
    if (nlu > lu )//|| (k % 20) == 0)
    {
      t++;
    }
    lu = nlu;
    //          std::cerr << nlu << std::endl;

  }

  return k;
}
