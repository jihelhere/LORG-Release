// -*- mode: c++ -*-
#include "WapitiDDApp.h"

#include "utils/LorgConstants.h"

#ifdef USE_THREADS
#include "utils/tick_count.h"
#endif

#include <thread>
#include "utils/WapitiWrapper.h"
#include "Word.h"
#include "Bracketing.h"
#include "utils/TokeniserFactory.h"


WapitiDDApp::WapitiDDApp() : LorgParseApp()
{}

WapitiDDApp::~WapitiDDApp()
{}


int WapitiDDApp::run()
{
  if(verbose)  std::clog << "Start tagging process.\n";

  std::string raw_sentence; // the sentence string (only for pretty print solutions)
  std::vector<Word> sentence; // the vector of words used by the parser to initialise the chart
  std::vector< bracketing > brackets; // pre-bracketting for guided parsing
  std::vector<std::string> comments; // strings to store comments associated with a sentence
  int count = 0; // sentence count (debug & pretty-print)

#ifdef USE_THREADS
  tick_count parse_start = tick_count::now();
#endif

  //read input and fill raw_sentence, sentence and brackets
  while(tokeniser->tokenise(*in, raw_sentence, sentence, brackets, comments)) {


    // // //Extra verbose
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


#ifdef USE_THREADS
    tick_count sent_start = tick_count::now();
#endif

    if(sentence.size() <=  max_length && sentence.size() > 0) {

      std::vector<std::thread> threads;


      // std::vector<std::vector<std::string>> crf_results(crfs.size());
      for (size_t i = 0; i < crfs.size(); ++i)
      {
        crfs[i].crf_tag();

        //std::cerr << crfs[i] << std::endl;
      }
      //      std::cerr << std::endl;



      int k = 0;

      if (crfs.size() > 1)
        k = find_consensus();

      if(verbose)
        std::cerr << "k: " << k << std::endl;

      ////

      for (size_t i = 0; i < 1 /*crfs.size()*/; ++i)
      {
        crfs[i].crf_retag();

        for (size_t j = 0; j < crfs[i].best_string_sequence.size(); ++j)
        {
          *out << sentence[j].get_form() << "\t" << crfs[i].best_string_sequence[j] << std::endl;
        }

        // for(const auto& s : crfs[i].best_string_sequence)
        // {
        //   *out << s << " " ;
        // }
        *out << std::endl;
      }

      ///


// #ifdef USE_THREADS
//       ///*
//       if(verbose && count % 50 == 0)
//         std::clog << count << " parsed sentences in " << (tick_count::now() - parse_start).seconds() << " sec" << std::endl;
//       //*/
// #endif


    }
    sentence.clear();
    brackets.clear();
    comments.clear();
    for (auto &crf : crfs) {crf.clean_sentence();}
  }

  *out << std::flush;

#ifdef USE_THREADS
  if(verbose) std::clog << "overall time: " << (tick_count::now() - parse_start).seconds() << "s" << std::endl;
#endif
  return 0; //everything's fine
}

LorgOptions WapitiDDApp::get_options() const
{
  LorgOptions options(LorgParseApp::get_options());
  options.add_2sparser_options();
  return options;
}


bool WapitiDDApp::read_config(ConfigTable& configuration)
{
  if(LorgParseApp::read_config(configuration) == false) return false;

  if(verbose) { std::clog << "creating the crf... ";}

  if (configuration.exists("crf-model"))
  {
    std::vector<std::string> crf_model_names = configuration.get_value<std::vector<std::string>>("crf-model");

    for (const auto& name : crf_model_names)
    {
      crfs.push_back(wapiti_wrapper(name));
    }


    std::vector<std::string> crf_input_names = configuration.get_value<std::vector<std::string>>("crf-input");
    for (auto i = 0U; i < crfs.size(); ++i)
    {
      crfs[i].set_file(crf_input_names[i]);
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

  for (size_t i = 0; i < tag_strings.size() - 1; ++i)
  {
    if ((tag_strings[i][1] == 'B') and (tag_strings[i+1][1] == 'I'))
    {
      results.push_back(i);
    }
  }
  return results;
}

// returns positions of the IB transitions
std::vector<unsigned> crf_get_mwe_ends(const std::vector<std::string>& tag_strings)
{
  std::vector<unsigned> results;
  bool i_before = false; // never an I at first

  for (size_t i = 1; i < tag_strings.size(); ++i)
  {
    if (tag_strings[i][1] == 'B' and i_before)
    {
      results.push_back(i-1);
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
  return results;
}


int WapitiDDApp::find_consensus()
{
  double lu = 0;
  bool valid = true;

  for (size_t i = 0; i < crfs.size(); ++i)
  {
    lu += crfs[i].score;
  }

  unsigned k = 0;
  double c = 1;
  double t = 0;

  for (k = 0; k < 1000; ++k)
  {
    //std::cout << "k = " << k << std::endl;
    // std::cout << "lu = " << lu << std::endl;
    // std::cout << "t = " << t << std::endl;
    double delta_k = c / (t + 1);
    double nlu = 0;

    // MW crfs
    std::vector<std::vector<unsigned>> mwe_starts(this->crfs.size());
    std::vector<std::vector<unsigned>> mwe_ends(this->crfs.size());

    std::map<unsigned,double> mwe_starts_average, mwe_ends_average;
    std::vector<MAP<unsigned, double>> mwe_starts_lambdas( this->crfs.size());
    std::vector<MAP<unsigned, double>> mwe_ends_lambdas(this->crfs.size());

    // for (size_t i = 0; i < parsers.size(); ++i)
    // {
    //   for (auto& e: parsers[i]->get_vectorized_representation(start_symbol))
    //   {
    //     //std::cout << std::get<1>(e) << " " << std::get<2>(e) << ": ";

    //     const AnnotatedRule* r = std::get<0>(e);

    //     int l;
    //     l  = simplification_map.at(r->get_lhs());

    //     if(SymbolTable::instance_nt().get_label_string(l)[0] != '[')
    //       // not an 'artificial node'
    //     {
    //       anchored_symbol vv(l, std::get<1>(e), std::get<2>(e));
    //       if(z_parser[i].insert(vv).second)
    //       {
    //         z_average[vv] += 1;
    //       }
    //     }

    //     // if (parsers[i]->get_is_funct() && (SymbolTable::instance_nt().get_label_string(l)[0] != '['))
    //     // {
    //     //   int f = r->get_lhs();

    //     //   anchored_symbol ff(f, std::get<1>(e), std::get<2>(e));
    //     //   if(f_parser[i].insert(ff).second)
    //     //   {
    //     //     f_average[ff] += 1;
    //     //   }
    //     // }

    //     //CRFS
    //     std::string s = SymbolTable::instance_nt().get_label_string(r->get_lhs());
    //     //if(s[0] == 'M' and s[1] == 'W') // MWE !
    //     if(s[s.size()-1] == '+')
    //     {
    //       mwe_starts[i].push_back(std::get<1>(e));
    //       mwe_ends[i].push_back(std::get<2>(e));

    //       //std::cout << i <<  " : " << s << " " << "(" << std::get<1>(e) << "," << std::get<2>(e) << ")" << std::endl;

    //       mwe_starts_average[std::get<1>(e)] +=1;
    //       mwe_ends_average[std::get<2>(e)] +=1;
    //     }
    //   }
    // }

    for (size_t i = 0; i < crfs.size(); ++i)
    {
      //      double score = k == 0 ? crfs[i].crf_tag() : crfs[i].crf_retag();

      const std::vector<std::string>& res = crfs[i].best_string_sequence;
      mwe_starts[i] = crf_get_mwe_beginings(res);
      mwe_ends[i] = crf_get_mwe_ends(res);
      for (size_t idx = 0; idx < mwe_starts[i].size(); ++idx)
        mwe_starts_average[mwe_starts[i][idx]] +=1;
      for (size_t idx = 0; idx < mwe_ends[i].size(); ++idx)
        mwe_ends_average[mwe_ends[i][idx]] +=1;
    }

    // check same solution
    bool same = true;

    for (auto& p : mwe_starts_average)
    {
      if (p.second != int(this->crfs.size()))
      {
        same = false;
        break;
      }
    }

    for (auto& p : mwe_ends_average)
    {
      if (p.second != int(this->crfs.size()))
      {
        same = false;
        break;
      }
    }


    if(same)
    {
      break;
    }


    for (size_t i = 0; i < this->crfs.size(); ++i)
    {
       mwe_starts_lambdas[i].clear();
      mwe_ends_lambdas[i].clear();

      for(auto& e : mwe_starts[i])
      {
        mwe_starts_lambdas[i][e] = delta_k * (1.0 - double(mwe_starts_average[e]) / (crfs.size()));
      }
      for(auto& e : mwe_ends[i])
      {
        mwe_ends_lambdas[i][e] = delta_k * (1.0 - double(mwe_ends_average[e]) / (crfs.size()));
      }

      // add update for spans missing in ith solution
      for(const auto& v: mwe_starts_average)
      {
        if(not mwe_starts_lambdas[i].count(v.first))
        {
          mwe_starts_lambdas[i][v.first] = delta_k * (- double(mwe_starts_average[v.first]) / (crfs.size()));
        }
      }
      for(const auto& v: mwe_ends_average)
      {
        if(not mwe_ends_lambdas[i].count(v.first))
        {
          mwe_ends_lambdas[i][v.first] = delta_k * (- double(mwe_ends_average[v.first]) / (crfs.size()));
        }
      }
    }



    // update relaxations
    std::vector<std::thread> threads;
    // MWE update
    threads.clear();
    for (size_t i = 0; i < this->crfs.size(); ++i)
    {
      threads.push_back(
          std::thread([&](int j)
                      {
                        crfs[j].update_relaxations(mwe_starts_lambdas[j], 'B', 'I', 1);
                        crfs[j].update_relaxations(mwe_ends_lambdas[j], 'I', 'B', 1);
                        },i));
    }
    for(auto& thread : threads)
    {
      thread.join();
    }

    threads.clear();
    for (size_t i = 0; i < crfs.size(); ++i)
    {
      threads.push_back(
          std::thread([&](int j)
                      {
                        crfs[j].crf_retag();
                      }, i)
                        );
    }

    for(auto& thread : threads)
    {
      thread.join();
    }


    for (size_t i = 0; i < crfs.size(); ++i)
    {
      std::cerr << "c: " << crfs[i].score << std::endl;
      nlu += crfs[i].score;
    }


    //        std::cout << std::endl;


    //update stepsize
    std::cerr << "nlu = " << nlu << std::endl;
    std::cerr << "lu = " << lu << std::endl;
    if (nlu > lu )//|| (k % 20) == 0)
    {
      ++t;
    }
    lu = nlu;
    //          std::cerr << nlu << std::endl;


    // for (const auto& crf : crfs)
    // {
    //   std::cerr << crf << std::endl;
    // }
    // std::cerr << std::endl;


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
