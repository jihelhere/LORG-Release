// -*- mode: c++ -*-
#include "TwoStageLorgParseApp.h"

#include "utils/Tagger.h"
#include "feature_extract/Extract.h"
#include "utils/PtbPsTree.h"
#include "utils/LorgConstants.h"
#include "parsers/ParserCKYAllFactory.h"

#include "parsers/ParserCKYAll.h"

#ifdef USE_THREADS
#include "utils/tick_count.h"
#endif

#include "lexicon/WordSignatureFactory.h"

#include <thread>

#include "utils/WapitiWrapper.h"


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

#ifdef USE_THREADS
  tick_count parse_start = tick_count::now();
#endif

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

    //   std::clog << "length " << sentence.size() << std::endl;
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
      // //std::cerr << "parse" << std::endl;
      parsers[i]->parse(start_symbol);
      // //std::cerr << "beam" << std::endl;
      parsers[i]->beam_c2f(start_symbol);
      //std::cerr << "after beam" << std::endl;
      if(parsers[i]->is_chart_valid(start_symbol))
      {
        //std::cerr << "extract" << std::endl;
        parsers[i]->extract_solution();
        //std::cerr << "after extract" << std::endl;
      }
    };


#ifdef USE_THREADS
    tick_count sent_start = tick_count::now();
#endif

    std::vector<std::pair<PtbPsTree *,double> > best_trees; // vector of (tree,score)

    if(sentence.size() <=  max_length && sentence.size() > 0) {

      std::vector<std::thread> threads;

      for(size_t i = 0; i < parsers.size(); ++i)
      {
        threads.push_back(std::thread(process_sentence,i));
      }

      for (size_t i = 0; i < bi_crfs.size(); ++i)
      {
        threads.push_back(std::thread([&](int j) {bi_crfs[j].crf_tag();}, i));

        // for(const auto& s : bi_crfs[i].best_string_sequence)
        // {
        //   std::cout << s << " " ;
        // }
        // std::cout << std::endl;
      }

      for (size_t i = 0; i < pos_crfs.size(); ++i)
      {
        threads.push_back(std::thread([&](int j) {pos_crfs[j].crf_tag();}, i));

        // for(const auto& s : bi_crfs[i].best_string_sequence)
        // {
        //   std::cout << s << " " ;
        // }
        // std::cout << std::endl;
      }

      for(auto& thread : threads)
      {
        thread.join();
      }


      int k = 0;

      if (parsers.size() + bi_crfs.size() + pos_crfs.size() > 1)
        k = find_consensus(best_trees, sentence.size()); // doesn't
                                                         // work with lattices

      // if(verbose)
      //   std::cerr << "k: " << k << std::endl;

      ////

      for (size_t i = 0; i < ( 0 == bi_crfs.size() ? 0 : 1); ++i)
      {
        bi_crfs[i].crf_retag();

        for (size_t j = 0; j < bi_crfs[i].best_string_sequence.size(); ++j)
        {
          std::cout << sentence[j].get_form() << "\t" << bi_crfs[i].best_string_sequence[j] << std::endl;
        }
        std::cout << std::endl;
      }


      for (size_t i = 0; i < ( 0 == pos_crfs.size() ? 0 : 1); ++i)
      {
        pos_crfs[i].crf_retag();

        for (size_t j = 0; j < pos_crfs[i].best_string_sequence.size(); ++j)
        {
          std::cout << sentence[j].get_form() << "\t" << pos_crfs[i].best_string_sequence[j] << std::endl;
        }
        std::cout << std::endl;
      }


      ///


      for (size_t i = 0; i < 1; ++i)
        //for (size_t i = 0; i < parsers.size(); ++i)
      {
        if(parsers[i]->is_chart_valid(start_symbol))
        {
          //BLOCKTIMING("get_parses");
          parsers[i]->get_parses(start_symbol, kbest, best_trees);
          //std::cout << "getting " << i << std::endl;
        }
        parse_solution * p_typed =
            parse_solution::factory.create_object(output_format,
                                                  parse_solution(raw_sentence, ++count,
                                                                 sentence.size(), best_trees,
#ifdef USE_THREADS
                                                                 (verbose) ? (tick_count::now() - sent_start).seconds() : 0,
#else
                                                                 0,
#endif
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


#ifdef USE_THREADS
      ///*
      if(verbose && count % 50 == 0)
        std::clog << count << " parsed sentences in " << (tick_count::now() - parse_start).seconds() << " sec" << std::endl;
      //*/
#endif


    }
    else
    {
        parse_solution * p_typed =
            parse_solution::factory.create_object(output_format,
                                                  parse_solution(raw_sentence, ++count,
                                                                 sentence.size(), best_trees,
#ifdef USE_THREADS
                                                                 (verbose) ? (tick_count::now() - sent_start).seconds() : 0,
#else
                                                                 0,
#endif
                                                                 verbose, comments, extract_features)
                                                  );
        p_typed->print(*out);
        delete p_typed;
    }

    sentence.clear();
    brackets.clear();
    comments.clear();
    for (auto &crf : bi_crfs) {crf.clean_sentence();}
    for (auto &crf : pos_crfs) {crf.clean_sentence();}
  }

  *out << std::flush;

#ifdef USE_THREADS
  if(verbose) std::clog << "overall time: " << (tick_count::now() - parse_start).seconds() << "s" << std::endl;
#endif
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
      p->set_word_signature(WordSignatureFactory::create_wordsignature(configuration));
    }
  }

  extract_features = configuration.get_value<bool>("extract-features");

  output_format = parse_solution::format_from_string(configuration.get_value<std::string>("output-format"));

  if (configuration.exists("crf-model"))
  {
    auto crf_model_names = configuration.get_value<std::vector<std::string>>("crf-model");
    auto crf_types = configuration.get_value<std::vector<std::string>>("crf-type");
    auto crf_input_names = configuration.get_value<std::vector<std::string>>("crf-input");
    auto crf_coefs = configuration.get_value<std::vector<int>>("crf-coefficient");

    // TODO: use boost::zip operator ?
    for (auto i = 0U; i < crf_model_names.size(); ++i)
    {
      wapiti_wrapper * wap_pt;

      if (crf_types[i] == "ftb_bi")
      {
        bi_crfs.push_back(wapiti_wrapper_bi(crf_model_names[i]));
        wap_pt = & bi_crfs.back();
      }
      else if (crf_types[i] == "pos")
      {
        pos_crfs.push_back(wapiti_wrapper_pos(crf_model_names[i]));
        wap_pt = & pos_crfs.back();
      }
      else
      {
        // meaningful error message here
        return false;
      }

      wap_pt->set_file(crf_input_names[i]);
      wap_pt->set_coefficient(crf_coefs[i]);

    }
  }

  return true;
}

//////////////////////////// DD /////////////////////////////

// assume the pos is @[BIO]@Category
// returns positions of the BI transitions
std::vector<unsigned> crf_get_mwe_beginings(const std::vector<std::string>& tag_strings)
{
  std::vector<unsigned> results;


  //std::cerr << "crf b: ";

  for (size_t i = 0; i < tag_strings.size() - 1; ++i)
  {
    if ((tag_strings[i][1] == 'B') and (tag_strings[i+1][1] == 'I'))
    {
      //        std::cerr << i << " ";
      results.push_back(i);
    }
  }
  //    std::cerr << std::endl;
  return results;
}

// returns positions of the IB transitions
std::vector<unsigned> crf_get_mwe_ends(const std::vector<std::string>& tag_strings)
{
  std::vector<unsigned> results;
  bool i_before = false; // never an I at first

  //  std::cerr << "crf e: ";


  for (size_t i = 1; i < tag_strings.size(); ++i)
  {
    if (tag_strings[i][1] == 'B' and i_before)
    {
      results.push_back(i-1);
      //       std::cerr << (i-1) << " ";
    }

    if (tag_strings[i][1] == 'I')
    {
      i_before = true;
    }
    else
    {
      i_before = false;
    }
  }

  //  std::cerr << std::endl;
  return results;
}


std::pair<std::vector<unsigned>,std::vector<unsigned>>
    pcfg_get_mwe_positions(ParserCKYAll* pparser)
{
  std::pair<std::vector<unsigned>,std::vector<unsigned>> results;

  // axiom of the grammar
  static int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  auto vect = pparser->get_vectorized_representation(start_symbol);

  for (auto& triple : vect)
  {
    const AnnotatedRule* r = std::get<0>(triple);
    std::string s = SymbolTable::instance_nt().get_label_string(r->get_lhs());
    //if(s[0] == 'M' and s[1] == 'W') // MWE !
    if(s[s.size()-1] == '+')
    {
      results.first.push_back(std::get<1>(triple));
      results.second.push_back(std::get<2>(triple));
    }
  }

  return results;
}



int TwoStageLorgParseApp::find_consensus(std::vector<std::pair<PtbPsTree *,double> >& /*best_trees*/, size_t sent_length)
{
  struct anchored_symbol
  {
    int symbol;
    int begin;
    int end;

    anchored_symbol(int n, int b, int e) : symbol(n), begin(b), end(e)
    {
      // std::cout << b << " " << e << " "
      //           << SymbolTable::instance_nt().get_label_string(n)
      //           << std::endl;
    };

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


  static int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name); // axiom of the grammar
  static std::unordered_map<int,int> simplification_map = SymbolTable::instance_nt().build_simplification_map();

  double lu = 0;
  bool valid = true;

  for (size_t i = 0; i < parsers.size(); ++i)
  {
    if(not parsers[i]->is_chart_valid(start_symbol))
      valid = false;
    else
      lu += parsers[i]->get_best_score(start_symbol);
  }

  for (size_t i = 0; i < bi_crfs.size(); ++i)
  {
    lu += bi_crfs[i].score * bi_crfs[i].coefficient;
  }

  for (size_t i = 0; i < pos_crfs.size(); ++i)
  {
    lu += pos_crfs[i].score * pos_crfs[i].coefficient;
  }

  if(!valid) return -1;

  unsigned k = 0;
  double c = 1;
  double t = 0;

  // size_t fun_parser_size = 0;
  // for(auto& p : parsers)
  // {
  //   if(p->get_is_funct()) ++fun_parser_size;
  // }

  bool use_ftb_bi = bi_crfs.size() > 0;
  bool use_pos = pos_crfs.size() > 0;


  for (k = 0; k < 1000; ++k)
  {
    //std::cout << "k = " << k << std::endl;
    //std::cerr << "lu = " << lu << std::endl;
    // std::cout << "t = " << t << std::endl;
    double delta_k = c / (t + 1);
    double nlu = 0;

    //std::cerr << "delta_k " << delta_k << std::endl;


    std::vector<std::set< anchored_symbol >> z_parser(this->parsers.size());
    std::map<anchored_symbol, int> z_average;
    std::vector<MAP<int,MAP<int,MAP<int,double>>>> z_lambdas(parsers.size());


    // // fun
    // std::vector<std::set< anchored_symbol >> f_parser(this->parsers.size());
    // std::map<anchored_symbol, int> f_average;
    // std::vector<MAP<int,MAP<int,MAP<int,double>>>> f_lambdas(parsers.size());


    // MW bi_crfs
    std::vector<std::vector<unsigned>> mwe_starts(this->parsers.size() + this->bi_crfs.size());
    std::vector<std::vector<unsigned>> mwe_ends(this->parsers.size() + this->bi_crfs.size());

    std::map<unsigned,double> mwe_starts_average , mwe_ends_average;
    std::vector<MAP<unsigned, double>> mwe_starts_lambdas(this->parsers.size() + this->bi_crfs.size());
    std::vector<MAP<unsigned, double>> mwe_ends_lambdas(this->parsers.size() + this->bi_crfs.size());


    // POS
    std::vector<std::set<anchored_symbol>>  sol_pos(this->parsers.size() + this->pos_crfs.size());
    std::map<anchored_symbol, int> sol_pos_average;
    std::vector<MAP<int,MAP<int,MAP<int, double>>>> pos_lambdas(this->parsers.size() + this->pos_crfs.size());


    for (size_t i = 0; i < parsers.size(); ++i)
    {
      for (auto& e: parsers[i]->get_vectorized_representation(start_symbol))
      {
        //std::cout << std::get<1>(e) << " " << std::get<2>(e) << ": ";

        const AnnotatedRule* r = std::get<0>(e);

        int l;
        l  = simplification_map.at(r->get_lhs());

        if(SymbolTable::instance_nt().get_label_string(l)[0] != '['
           &&
           (!r->is_lexical() || !use_pos) // dont process lexicals twice
           )
          // not an 'artificial node'
        {
          anchored_symbol vv(l, std::get<1>(e), std::get<2>(e));
          if(z_parser[i].insert(vv).second)
          {
            z_average[vv] += 1;
          }
        }

        // if (parsers[i]->get_is_funct() && (SymbolTable::instance_nt().get_label_string(l)[0] != '['))
        // {
        //   int f = r->get_lhs();

        //   anchored_symbol ff(f, std::get<1>(e), std::get<2>(e));
        //   if(f_parser[i].insert(ff).second)
        //   {
        //     f_average[ff] += 1;
        //   }
        // }

        //BI_CRFS
        if (use_ftb_bi)
        {
          std::string s = SymbolTable::instance_nt().get_label_string(r->get_lhs());
          //if(s[0] == 'M' and s[1] == 'W') // MWE !
          if(s[s.size()-1] == '+')
          {
            // std::cerr
            //     << "parser: "
            //     << "(" << std::get<1>(e)
            //     << "," << std::get<2>(e)
            //     << ") " << sent_length << std::endl;

            mwe_starts[i].push_back(std::get<1>(e));
            mwe_starts_average[std::get<1>(e)] +=1;


            if (std::get<2>(e) >= int(sent_length) -1)
            {
              //std::cerr << "SENT_LENGTH" << std::endl;
            }
            else
            {
              mwe_ends[i].push_back(std::get<2>(e));
              mwe_ends_average[std::get<2>(e)] +=1;
            }

            //std::cout << i <<  " : " << s << " " << "(" << std::get<1>(e) << "," << std::get<2>(e) << ")" << std::endl;
          }
        }


        // POS
        if (use_pos)
        {
          if (r->is_lexical())
          {
            anchored_symbol sym(l/*r->get_lhs()*/, std::get<1>(e), std::get<2>(e));
            sol_pos[i].insert(sym);
            sol_pos_average[sym] += 1;

            // std::string st = SymbolTable::instance_nt().get_label_string(l);
            // std::cerr << i <<  " : " << st << " " << "(" << std::get<1>(e) << "," << std::get<2>(e) << ")" << std::endl;
          }
        }

      }
    }


    if (use_ftb_bi)
    {
      for (size_t i = 0; i < bi_crfs.size(); ++i)
      {
        //      double score = k == 0 ? bi_crfs[i].crf_tag() : bi_crfs[i].crf_retag();

        const std::vector<std::string>& res = bi_crfs[i].best_string_sequence;
        mwe_starts[this->parsers.size() + i] = crf_get_mwe_beginings(res);
        mwe_ends[this->parsers.size() + i] = crf_get_mwe_ends(res);
        for (size_t idx = 0; idx < mwe_starts[this->parsers.size() + i].size(); ++idx)
          mwe_starts_average[mwe_starts[this->parsers.size() + i][idx]] +=1*bi_crfs[i].coefficient;
        for (size_t idx = 0; idx < mwe_ends[this->parsers.size() + i].size(); ++idx)
          mwe_ends_average[mwe_ends[this->parsers.size() + i][idx]] +=1*bi_crfs[i].coefficient;
      }
    }

    if (use_pos)
    {
      for (auto i = 0U; i < pos_crfs.size(); ++i)
      {
        const std::vector<std::string>& res = pos_crfs[i].best_string_sequence;

        for (auto j = 0U; j < res.size(); ++j)
        {
          anchored_symbol s(SymbolTable::instance_nt().get_label_id(res[j]), j, j);
          sol_pos[parsers.size()+i].insert(s);
          sol_pos_average[s] += 1*pos_crfs[i].coefficient;
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


    // for (auto& p : f_average)
    // {
    //   if (p.second != int(fun_parser_size))
    //   {
    //     //std::cout << p.second << std::endl;
    //     same = false;
    //     break;
    //   }
    // }


    int sum_coefficient_bi_crf = std::accumulate(bi_crfs.begin(), bi_crfs.end(), 0, [] (int acc, const wapiti_wrapper& c) {return acc + c.coefficient;});
    int sum_coefficient_pos_crf = std::accumulate(pos_crfs.begin(), pos_crfs.end(), 0, [] (int acc, const wapiti_wrapper& c) {return acc + c.coefficient;});
    //    std::cerr << "sum_coefficient_crf " << sum_coefficient_crf << std::endl;


    if (use_ftb_bi)
    {
      for (auto& p : mwe_starts_average)
      {
        if (p.second != this->parsers.size() + sum_coefficient_bi_crf)
        {
          same = false;
          break;
        }
      }
      for (auto& p : mwe_ends_average)
      {
        if (p.second != this->parsers.size() + sum_coefficient_bi_crf)
        {
          same = false;
          break;
        }
      }
    }

    if (use_pos)
    {
      for (auto& p : sol_pos_average)
      {
        if (p.second != int(this->parsers.size()) + sum_coefficient_pos_crf)
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

    for(size_t i = 0; i < z_parser.size(); ++i)
    {
      z_lambdas[i].clear();

      // add update for spans in ith solution
      for(const auto& v: z_parser[i])
      {
        z_lambdas[i][v.begin][v.end][v.symbol]
            = delta_k * (1.0 - double(z_average[v]) / parsers.size());
      }

      // add update for spans missing in ith solution
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


    // for(size_t i = 0; i < f_parser.size(); ++i)
    // {
    //   if(not parsers[i]->get_is_funct()) continue;

    //   f_lambdas[i].clear();

    //   // add update for things in ith solution
    //   for(const auto& v: f_parser[i])
    //   {
    //     f_lambdas[i][v.begin][v.end][v.symbol]
    //         = delta_k * (1.0 - double(f_average[v]) / fun_parser_size);
    //   }

    //   // add update for things missing in ith solution
    //   for(const auto& v: f_average)
    //   {
    //     if(not f_lambdas[i].count(v.first.begin)
    //        or not f_lambdas[i][v.first.begin].count(v.first.end)
    //        or not f_lambdas[i][v.first.begin][v.first.end].count(v.first.symbol)
    //        )
    //     {
    //       f_lambdas[i][v.first.begin][v.first.end][v.first.symbol]x
    //           = delta_k * (- double(f_average[v.first]) / fun_parser_size);
    //     }
    //   }
    // }



    // MWE \delta 's
    if (use_ftb_bi)
    {
      for (size_t i = 0; i < this->parsers.size() + this->bi_crfs.size(); ++i)
      {
        mwe_starts_lambdas[i].clear();
        mwe_ends_lambdas[i].clear();

        for(auto& e : mwe_starts[i])
        {
          mwe_starts_lambdas[i][e] = delta_k * (1.0 - double(mwe_starts_average[e]) / (parsers.size() + sum_coefficient_bi_crf));
        }
        for(auto& e : mwe_ends[i])
        {
          mwe_ends_lambdas[i][e] = delta_k * (1.0  - double(mwe_ends_average[e]) / (parsers.size() + sum_coefficient_bi_crf));
        }

        // add update for spans missing in ith solution
        for(const auto& v: mwe_starts_average)
        {
          if(not mwe_starts_lambdas[i].count(v.first))
          {
            mwe_starts_lambdas[i][v.first] = delta_k * (- double(mwe_starts_average[v.first]) / (parsers.size() + sum_coefficient_bi_crf));
          }
        }
        for(const auto& v: mwe_ends_average)
        {
          if(not mwe_ends_lambdas[i].count(v.first))
          {
            mwe_ends_lambdas[i][v.first] = delta_k * (- double(mwe_ends_average[v.first]) / (parsers.size() + sum_coefficient_bi_crf));
          }
        }
      }
    }

    // POS \delta 's
    if (use_pos)
    {
      for (size_t i = 0; i < this->parsers.size() + this->pos_crfs.size(); ++i)
      {
        pos_lambdas[i].clear();

        for (auto& e : sol_pos[i])
        {
          pos_lambdas[i][e.begin][e.end][e.symbol] =  delta_k * (1.0 - double(sol_pos_average[e]) / (parsers.size() + sum_coefficient_pos_crf));
        }


        // add update for lexical spans missing in ith solution
        for(const auto& v: sol_pos_average)
        {
          if(not pos_lambdas[i].count(v.first.begin)
             or not pos_lambdas[i][v.first.begin].count(v.first.end)
             or not pos_lambdas[i][v.first.begin][v.first.end].count(v.first.symbol)
             )
          {
            pos_lambdas[i][v.first.begin][v.first.end][v.first.symbol] = delta_k * (- double(sol_pos_average[v.first]) / (parsers.size() + sum_coefficient_pos_crf));
          }
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
                          this->parsers[j]->update_relaxations(true,  z_lambdas[j], simplification_map);
                          // if(parsers[j]->get_is_funct())
                          //   this->parsers[j]->update_relaxations(false, f_lambdas[j], simplification_map);
                        }
                        },i));
    }
    for(auto& thread : threads)
    {
      thread.join();
    }

    // MWE update
    if (use_ftb_bi)
    {
      threads.clear();
      for (size_t i = 0; i < this->parsers.size(); ++i)
      {
        threads.push_back(
            std::thread([&](int j)
                        {
                          if(this->parsers[j]->is_chart_valid(start_symbol))
                          {
                            this->parsers[j]->update_relaxations_starts(mwe_starts_lambdas[j]);
                            this->parsers[j]->update_relaxations_ends(mwe_ends_lambdas[j]);
                          }

                        },i));
      }
      // for(auto& thread : threads)
      // {
      //   thread.join();
      // }


      // MWEs
      //threads.clear();
      for (size_t i = 0; i < this->bi_crfs.size(); ++i)
      {
        threads.push_back(
            std::thread([&](int j)
                        {
                          bi_crfs[j].update_relaxations(mwe_starts_lambdas[this->parsers.size() + j], 'B', 'I', 1);
                          bi_crfs[j].update_relaxations(mwe_ends_lambdas[this->parsers.size() + j], 'I', 'B', 1);
                        },i));

        // bi_crfs[i].update_relaxations(mwe_starts_lambdas[this->parsers.size() + i], 'B', 'I', 1);
        // bi_crfs[i].update_relaxations(mwe_ends_lambdas[this->parsers.size() + i], 'I', 'B', 1);

      }
      for(auto& thread : threads)
      {
        thread.join();
      }
    }

    // POS
    if (use_pos)
    {
      threads.clear();
      for (size_t i = 0; i < this->parsers.size(); ++i)
      {
        threads.push_back(
            std::thread([&](int j)
                        {
                          if(this->parsers[j]->is_chart_valid(start_symbol))
                          {
                            this->parsers[j]->update_relaxations(true,  pos_lambdas[j], simplification_map);
                          }

                        },i));
      }
      // for(auto& thread : threads)
      // {
      //   thread.join();
      // }


      //threads.clear();
      for (size_t i = 0; i < this->pos_crfs.size(); ++i)
      {
        threads.push_back(
            std::thread([&](int j)
                        {
                          pos_crfs[j].update_relaxations(pos_lambdas[this->parsers.size() + j]);
                        },i));

        // bi_crfs[i].update_relaxations(mwe_starts_lambdas[this->parsers.size() + i], 'B', 'I', 1);
        // bi_crfs[i].update_relaxations(mwe_ends_lambdas[this->parsers.size() + i], 'I', 'B', 1);

      }

      for(auto& thread : threads)
      {
        thread.join();
      }
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




    for (size_t i = 0; i < bi_crfs.size(); ++i)
    {
      threads.push_back(
          std::thread([&](int j)
                      {
                        bi_crfs[j].crf_retag();
                      }, i)
                        );
    }


    for (size_t i = 0; i < pos_crfs.size(); ++i)
    {
      threads.push_back(
          std::thread([&](int j)
                      {
                        pos_crfs[j].crf_retag();
                      }, i)
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
        //std::cerr << "p: " <<  parsers[i]->get_best_score(start_symbol) << std::endl;
        nlu += parsers[i]->get_best_score(start_symbol);
        // if (i == 0)
        // {
        //   parsers[i]->get_parses(start_symbol, 1, always_output_forms, output_annotations, best_trees);
        // }
      }
    }

    for (size_t i = 0; i < bi_crfs.size(); ++i)
    {
      //std::cerr << "c_bi: " <<  bi_crfs[i].score << std::endl;
      nlu += bi_crfs[i].score * bi_crfs[i].coefficient;
    }

    for (size_t i = 0; i < pos_crfs.size(); ++i)
    {
      //std::cerr << "c_pos: " <<  pos_crfs[i].score << std::endl;
      nlu += pos_crfs[i].score * pos_crfs[i].coefficient;
    }

    //        std::cout << std::endl;


    //update stepsize
    // std::cerr << "nlu = " << nlu << std::endl;
    // std::cerr << "lu = " << lu << std::endl;
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
