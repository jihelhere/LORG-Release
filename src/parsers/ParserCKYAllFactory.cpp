// -*- mode: c++ -*-

#include "ParserCKYAllFactory.h"

#include "utils/ConfigTable.h"
#include "viterbi/ParserCKYAllViterbi.h"
#include "maxrule/ParserCKYAllMaxRule1B.h"
#include "maxrule/ParserCKYAllMaxRuleKB.h"
#include "maxrule/ParserCKYAllMaxRuleMultiple.h"
#include "mindivkb/ParserCKYAllMinDivKB.h"
#include "variational/ParserCKYAllVariational.h"

#include "grammars/GrammarAnnotated.hpp"
#include "utils/data_parsers/AnnotHistoriesParser.h"

#include "lexicon/WordSignatureFactory.h"

#include <fstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "utils/tuple_serialization.hpp"


typedef boost::archive::text_oarchive oarchive;
typedef boost::archive::text_iarchive iarchive;



ParserCKYAllFactory::Parsing_Algorithm
ParserCKYAllFactory::string_to_pa(const std::string& s)
{
    if (s == "max") return ParserCKYAllFactory::MaxRule;
    if (s == "vit") return ParserCKYAllFactory::Viterbi;

    if (s == "maxn") return ParserCKYAllFactory::MaxN;
    if (s == "kmax") return ParserCKYAllFactory::KMaxRule;

    if (s == "mind") return ParserCKYAllFactory::MinDiv;

    if (s == "var") return ParserCKYAllFactory::Variational;

    std::clog << "Unknown parser type: " << s << std::endl;

    throw std::out_of_range(s);
}


ParserCKYAllFactory::MaxParsing_Calculation
ParserCKYAllFactory::string_to_mpc(const std::string& s)
{
    if (s == "product") return ParserCKYAllFactory::Product;
    if (s == "sum") return ParserCKYAllFactory::Sum;
    if (s == "prodsum") return ParserCKYAllFactory::ProdSum;

    std::clog << "Unknown calculation type: " << s << std::endl;

    throw std::out_of_range(s);
}


ParserCKYAll * create_parser(std::vector<ParserCKYAll::AGrammar*> cgs,
                             ParserCKYAllFactory::Parsing_Algorithm pt,
                             ParserCKYAllFactory::MaxParsing_Calculation c,
                             const std::vector<double>& p, double b_t,
                             const std::vector< std::vector<ParserCKYAll::AGrammar*> >& fgs,
                             const std::vector< annot_descendants_type >& all_annot_descendants,
                             bool accurate, unsigned min_beam, int stubborn,
                             unsigned k)
{
    using namespace ParserCKYAllFactory;
    switch(pt) {
      case Viterbi :
        return new ParserCKYAllViterbi(cgs, p, b_t, all_annot_descendants[0], accurate, min_beam, stubborn);
      case MaxRule :
        return new ParserCKYAllMaxRule1B(c, cgs, p, b_t, all_annot_descendants[0], accurate, min_beam, stubborn);
      case MaxN :
        return new ParserCKYAllMaxRuleMultiple(c, cgs, p, b_t, fgs, all_annot_descendants, accurate, min_beam, stubborn, k);
      case KMaxRule :
        return new ParserCKYAllMaxRuleKB(c, cgs, p, b_t, all_annot_descendants[0], accurate, min_beam, stubborn, k);
      case MinDiv :
        return new ParserCKYAllMinDivKB(cgs, p, b_t, all_annot_descendants[0], accurate, min_beam, stubborn, k);
      case Variational :
        return new ParserCKYAllVariational(cgs, p, b_t, all_annot_descendants[0], accurate, min_beam, stubborn, k);
      default :
        return NULL;
    }
}



std::vector<ParserCKYAll::AGrammar*>
create_intermediates(ParserCKYAll::AGrammar& grammar, const annot_descendants_type& annot_descendants);

annot_descendants_type
create_annot_descendants(const std::map<short, Tree<unsigned> >& annot_histories);


std::vector<ParserCKYAll::AGrammar*> create_grammars(const std::string& filename, bool verbose)
{
  std::vector<ParserCKYAll::AGrammar*> grammars;

  std::vector<annot_descendants_type> all_annot_descendants;

  // compute expected intermediate grammars

  //get grammar
  if(verbose) std::cerr << "Setting grammar to " << filename << ".\n";
  ParserCKYAll::AGrammar * cg = new ParserCKYAll::AGrammar(filename);
  if(verbose) std::cerr << "Grammar set\n";


  if(verbose) std::clog << "create intermediate grammars" << std::endl;
  annot_descendants_type annot_descendants = create_annot_descendants(cg->get_history_trees());
  grammars = create_intermediates(*cg, annot_descendants);
  //  std::clog << "init" << std::endl;
  cg->init();
  //    std::clog << "smooth" << std::endl;
  // TODO: options to set this from command line
  cg->linear_smooth(0.01,0.1);
  //    std::clog << "clean" << std::endl;
  cg->remove_unlikely_annotations_all_rules(1e-10);
  grammars.push_back(cg);
  //    std::clog << "creation finished" << std::endl;

  ///////

  if(verbose)
  {
    std::tuple<decltype(grammars)&,
               decltype(SymbolTable::instance_nt())&,
               decltype(SymbolTable::instance_word())&> triple(grammars, SymbolTable::instance_nt(), SymbolTable::instance_word());
    std::ofstream ofs(filename + ".bin");
    oarchive oa(ofs);
    oa << triple;
  }
  //////////


  return grammars;
}

std::vector<ParserCKYAll *>
ParserCKYAllFactory::create_parser(ConfigTable& config)
{
    bool verbose = config.exists("verbose");

#ifdef USE_THREADS
    unsigned nbthreads = config.get_value<unsigned>("nbthreads");
#endif

    if(verbose)
    {
      std::clog << "using " <<
#ifndef USE_THREADS
          1
#else
          (nbthreads == 0 ? tbb::task_scheduler_init::default_num_threads() : nbthreads)
#endif
                << " thread(s) to parse" << std::endl;
    }

    // the vector of grammars to be returned
    std::vector<ParserCKYAll*> results;


    double beam_threshold = config.get_value<double>("beam-threshold");
    if(verbose) std::clog << "using threshold " << beam_threshold << " for chart construction" << std::endl;

    bool accurate = config.exists("accurate");

    if(verbose)
    {
        if(accurate)
            std::clog << "using accurate c2f thresholds" << std::endl;
        else
            std::clog << "using standard c2f thresholds" << std::endl;
    }

    unsigned min = config.get_value<unsigned>("min-length-beam");

    /// utils functions
    std::function<void(std::vector<ParserCKYAll::AGrammar*>&, const SymbolTable&, const SymbolTable&)> align_grammar
        =[](std::vector<ParserCKYAll::AGrammar*>& gs, const SymbolTable& words, const SymbolTable& nts)
        {
          for (auto& g : gs)
          {
            for (auto& br: g->binary_rules)
            {
              int l = br.get_lhs();
              std::string ls = nts.get_label_string(l);
              int l2 = SymbolTable::instance_nt().insert(ls);
              br.set_lhs(l2);

              int r0 = br.get_rhs0();
              std::string r0s = nts.get_label_string(r0);
              int r02 = SymbolTable::instance_nt().insert(r0s);
              br.set_rhs0(r02);

              int r1 = br.get_rhs1();
              std::string r1s = nts.get_label_string(r1);
              int r12 = SymbolTable::instance_nt().insert(r1s);
              br.set_rhs1(r12);
            }

            for (auto& ur: g->unary_rules)
            {
              int l = ur.get_lhs();
              std::string ls = nts.get_label_string(l);
              int l2 = SymbolTable::instance_nt().insert(ls);
              ur.set_lhs(l2);

              int r0 = ur.get_rhs0();
              std::string r0s = nts.get_label_string(r0);
              int r02 = SymbolTable::instance_nt().insert(r0s);
              ur.set_rhs0(r02);
            }

            for (auto& lr: g->lexical_rules)
            {
              int l = lr.get_lhs();
              std::string ls = nts.get_label_string(l);
              int l2 = SymbolTable::instance_nt().insert(ls);
              lr.set_lhs(l2);

              int r0 = lr.get_rhs0();
              std::string r0s = words.get_label_string(r0);
              int r02 = SymbolTable::instance_word().insert(r0s);
              lr.set_rhs0(r02);
            }


            // align history trees and num_annotations
            auto& h = g->get_history_trees();
            //use decltype(h) :  how to remove & ??
            std::map< short, Tree<unsigned> > h2;
            for(const auto& p: h)
            {
              auto pfs = nts.get_label_string(p.first);
              auto pf2 = SymbolTable::instance_nt().get_label_id(pfs);
              h2[pf2] = p.second;
            }
            h = h2;

            //
            auto& ai = g->get_annotations_info();
            std::vector<short unsigned> v(SymbolTable::instance_nt().get_size(),0);
            for (size_t i = 0; i < ai.get_number_of_unannotated_labels(); ++i)
            {
              auto is = nts.get_label_string(i);
              auto i2 = SymbolTable::instance_nt().get_label_id(is);
              v[i2] = ai[i];
            }

            g->get_annotations_info().set_num_annotations_map(v);


            //unary chains
            PathMatrix& old = g->get_unary_decoding_paths();
            PathMatrix new_path;
            for(auto i : old)
            {
              auto is = nts.get_label_string(i.first.first);
              auto i2 = SymbolTable::instance_nt().get_label_id(is);
              for(auto j : i.second)
              {

                auto js = nts.get_label_string(j.first.first);
                auto j2 = SymbolTable::instance_nt().get_label_id(js);

                auto ks = nts.get_label_string(j.second.first);
                auto k2 = SymbolTable::instance_nt().get_label_id(ks);

                new_path[{i2,i.first.second}][{j2,j.first.second}] = {k2,j.second.second};
              }
            }

            old = new_path;
          }
        };


    std::function<std::vector<ParserCKYAll::AGrammar*>(std::vector<annot_descendants_type>&, int)>
        build_grammar =
        [&config, &verbose,&align_grammar](std::vector<annot_descendants_type>& my_all_annot_descendants, int parser_idx)
        {
          std::vector<ParserCKYAll::AGrammar*> my_grammars;
          // get my_grammars

          // hackish
          std::string str_idx("1");
          str_idx[0] += parser_idx;

          if(config.exists("grammar" + str_idx))
          {
            const std::string& filename = config.get_value< std::string >("grammar" + str_idx);
            my_grammars = create_grammars(filename, verbose);
          }
          else
          {
            if(config.exists("archive-grammar"+str_idx))
            {
                if(verbose) std::cerr << "loading binary file" << std::endl;
                std::tuple<decltype(my_grammars), SymbolTable, SymbolTable> triple;
                std::ifstream ifs(config.get_value<std::string>("archive-grammar" + str_idx));
                iarchive(ifs) >> triple;
                my_grammars = std::get<0>(triple);

                if(parser_idx == 0)
                {// grammars for first parser
                  SymbolTable::instance_nt().load(std::get<1>(triple));
                  SymbolTable::instance_word().load(std::get<2>(triple));
                }
                else // parser_idx != 0
                {
                  align_grammar(my_grammars, std::get<2>(triple), std::get<1>(triple));
                }
            }
          }
          if(my_grammars.size() == 0)
          {
            std::cerr << "Grammar wasn't set. Exit program." << std::endl;
            throw std::runtime_error("Grammar not set");
          }

          annot_descendants_type annot_descendants = create_annot_descendants(my_grammars.back()->get_history_trees());
          my_all_annot_descendants.push_back(annot_descendants);

          return my_grammars;
        };

    std::function<std::vector<std::vector<ParserCKYAll::AGrammar*>>(std::vector<annot_descendants_type>&, int)>
        build_alternate =
        [&config,&verbose,&align_grammar](std::vector<annot_descendants_type>& my_all_annot_descendants, int parser_idx)
        {
          std::vector<std::vector<ParserCKYAll::AGrammar*>> my_alts;

          // hackish
          std::string str_idx("1");
          str_idx[0] += parser_idx;

          // if(config.exists("alternate-grammar" + str_idx) or config.exists("archive-alternategrammars" + str_idx))
          // {
          const std::vector<std::string>& filenames = config.exists("alternate-grammar" + str_idx) ?
          config.get_value<std::vector<std::string> >("alternate-grammar" + str_idx) :
          config.get_value<std::vector<std::string> >("archive-alternategrammars" + str_idx);
          if(string_to_pa(config.get_value<std::string>("parser-type")) != MaxN && filenames.size() > 0)
          {
            throw std::runtime_error("incorrect parsing algorithm");
          }

          if(config.exists("alternate-grammar" + str_idx))
          {
            for(size_t i = 0; i < filenames.size(); ++i)
            {
              if(verbose) std::cerr << "Setting alternate grammar to " << filenames[i] << ".\n";
              my_alts.push_back(create_grammars(filenames[i], verbose));
              annot_descendants_type annot_descendants = create_annot_descendants(my_alts.back().back()->get_history_trees());
              my_all_annot_descendants.push_back(annot_descendants);
            }
          }
          if(config.exists("archive-alternategrammars" + str_idx)) {
            for(size_t i = 0; i < filenames.size(); ++i)
            {
              if(verbose) std::cerr << "Setting alternate grammar to " << filenames[i] << ".\n";

              std::vector<ParserCKYAll::AGrammar*> grammars;

              std::tuple<decltype(grammars), SymbolTable, SymbolTable> triple;
              std::ifstream ifs(filenames[i]);
              iarchive(ifs) >> triple;
              grammars = std::get<0>(triple);

              if(parser_idx != 0)
              {
                align_grammar(grammars, std::get<2>(triple), std::get<1>(triple));
              }

              my_alts.push_back(grammars);
              annot_descendants_type annot_descendants = create_annot_descendants(grammars.back()->get_history_trees());
              my_all_annot_descendants.push_back(annot_descendants);
            }
          }

          return my_alts;
        };

    ////
    {
      std::vector<double> priors;
      std::vector<annot_descendants_type> all_annot_descendants;
      std::vector<ParserCKYAll::AGrammar*> grammars;

      grammars = build_grammar(all_annot_descendants, 0);
      // compute priors for base grammar
      if(verbose) std::cerr << "computing priors" << std::endl;
      priors = grammars[0]->compute_priors();

      std::vector< std::vector<ParserCKYAll::AGrammar*> > alt_gs;

      if(config.exists("alternate-grammar1") or config.exists("archive-alternategrammars1"))
      {
        alt_gs = build_alternate(all_annot_descendants, 0);
      }

      results.push_back( create_parser(grammars, string_to_pa(config.get_value<std::string>("parser-type")),
                                       string_to_mpc(config.get_value<std::string>("max-type")), priors, beam_threshold,
                                       alt_gs, all_annot_descendants, accurate, min, config.get_value<int>("stubbornness"),
                                       config.get_value<unsigned>("kbest")));

      if (grammars.back()->gram_conf.count("lex_unknown_map"))
      {
        std::cerr << "gram_conf contains lex_unknown_map info" << std::endl;
        if (results.back()->get_word_signature() == nullptr)
        {
          std::cerr << "overwriting unknown_map from command-line (if you don't want this, edit the grammar)" << std::endl;
          results.back()->set_word_signature(
              WordSignatureFactory::create_wordsignature(
                  WordSignature::string_2_lex_unknown_map(grammars.back()->gram_conf.at("lex_unknown_map")),
                  true));
        }
      }
      if (grammars.back()->gram_conf.count("remove_functions"))
      {
        if(grammars.back()->gram_conf.at("remove_functions") == "0")
        {
          results.back()->set_is_funct(true);
        }
      }
    }

    {
      std::vector<double> priors2;
      std::vector<annot_descendants_type> all_annot_descendants2;
      std::vector<ParserCKYAll::AGrammar*> grammars2;
      std::vector< std::vector<ParserCKYAll::AGrammar*> > alt_gs2;
      // get grammars
      if(config.exists("grammar2") or config.exists("archive-grammar2"))
      {
        grammars2 = build_grammar(all_annot_descendants2, 1);
        // compute priors for base grammar
        priors2 = grammars2[0]->compute_priors();

        ///////////////////
        if(config.exists("alternate-grammar2") or config.exists("archive-alternategrammars2"))
        {
          alt_gs2 = build_alternate(all_annot_descendants2,1);
        }

        ///////////////////


        results.push_back(create_parser(grammars2, string_to_pa(config.get_value<std::string>("parser-type")),
                                        string_to_mpc(config.get_value<std::string>("max-type")), priors2, beam_threshold,
                                        alt_gs2, all_annot_descendants2, accurate, min, config.get_value<int>("stubbornness"),
                                        config.get_value<unsigned>("kbest")));

        if (grammars2.back()->gram_conf.count("lex_unknown_map"))
        {
          if (results.back()->get_word_signature() == nullptr)
          {
            std::cerr << "overwriting unknown_map from command-line (if you don't want this, edit the grammar)" << std::endl;
            results.back()->set_word_signature(
                WordSignatureFactory::create_wordsignature(
                    WordSignature::string_2_lex_unknown_map(grammars2.back()->gram_conf.at("lex_unknown_map")),
                    true));
          }
        }
        if (grammars2.back()->gram_conf.count("remove_functions"))
        {
          if(grammars2.back()->gram_conf.at("remove_functions") == "0")
          {
            results.back()->set_is_funct(true);
          }
        }
      }
    }


    ////////////////
    {
      std::vector<double> priors3;
      std::vector<annot_descendants_type> all_annot_descendants3;
      std::vector<ParserCKYAll::AGrammar*> grammars3;
      std::vector< std::vector<ParserCKYAll::AGrammar*> > alt_gs3;
      // get grammars
      if(config.exists("grammar3") or config.exists("archive-grammar3")) {

        grammars3 = build_grammar(all_annot_descendants3, 2);
        // compute priors for base grammar
        priors3 = grammars3[0]->compute_priors();

        ///////////////////
        if(config.exists("alternate-grammar3") or config.exists("archive-alternategrammars3"))
        {
          alt_gs3 = build_alternate(all_annot_descendants3,2);
        }

        ///////////////////
        results.push_back(create_parser(grammars3, string_to_pa(config.get_value<std::string>("parser-type")),
                                        string_to_mpc(config.get_value<std::string>("max-type")), priors3, beam_threshold,
                                        alt_gs3, all_annot_descendants3, accurate, min, config.get_value<int>("stubbornness"),
                                        config.get_value<unsigned>("kbest")));

        if (grammars3.back()->gram_conf.count("lex_unknown_map"))
        {
          if (results.back()->get_word_signature() == nullptr)
          {
            std::cerr << "overwriting unknown_map from command-line (if you don't want this, edit the grammar)" << std::endl;
            results.back()->set_word_signature(
                WordSignatureFactory::create_wordsignature(
                    WordSignature::string_2_lex_unknown_map(grammars3.back()->gram_conf.at("lex_unknown_map")),
                    true));
          }
        }
        if (grammars3.back()->gram_conf.count("remove_functions"))
        {
          if(grammars3.back()->gram_conf.at("remove_functions") == "0")
          {
            results.back()->set_is_funct(true);
          }
        }
      }
    }


    ////////////////
    {
      std::vector<double> priors4;
      std::vector<annot_descendants_type> all_annot_descendants4;
      std::vector<ParserCKYAll::AGrammar*> grammars4;
      std::vector< std::vector<ParserCKYAll::AGrammar*> > alt_gs4;
      // get grammars
      if(config.exists("grammar4") or config.exists("archive-grammar4")) {

        grammars4 = build_grammar(all_annot_descendants4, 3);
        // compute priors for base grammar
        priors4 = grammars4[0]->compute_priors();

        ///////////////////
        if(config.exists("alternate-grammar4") or config.exists("archive-alternategrammars4"))
        {
          alt_gs4 = build_alternate(all_annot_descendants4,3);
        }

        ///////////////////
        results.push_back(create_parser(grammars4, string_to_pa(config.get_value<std::string>("parser-type")),
                                        string_to_mpc(config.get_value<std::string>("max-type")), priors4, beam_threshold,
                                        alt_gs4, all_annot_descendants4, accurate, min, config.get_value<int>("stubbornness"),
                                        config.get_value<unsigned>("kbest")));

        if (grammars4.back()->gram_conf.count("lex_unknown_map"))
        {
          if (results.back()->get_word_signature() == nullptr)
          {
          std::cerr << "overwriting unknown_map from command-line (if you don't want this, edit the grammar)" << std::endl;
          results.back()->set_word_signature(
              WordSignatureFactory::create_wordsignature(
                  WordSignature::string_2_lex_unknown_map(grammars4.back()->gram_conf.at("lex_unknown_map")),
                  true));
          }
        }
        if (grammars4.back()->gram_conf.count("remove_functions"))
        {
          if(grammars4.back()->gram_conf.at("remove_functions") == "0")
          {
            results.back()->set_is_funct(true);
          }
        }
      }
    }


    ////////////////
    {
      std::vector<double> priors5;
      std::vector<annot_descendants_type> all_annot_descendants5;
      std::vector<ParserCKYAll::AGrammar*> grammars5;
      std::vector< std::vector<ParserCKYAll::AGrammar*> > alt_gs5;
      // get grammars
      if(config.exists("grammar5") or config.exists("archive-grammar5")) {

        grammars5 = build_grammar(all_annot_descendants5, 3);
        // compute priors for base grammar
        priors5 = grammars5[0]->compute_priors();

        ///////////////////
        if(config.exists("alternate-grammar5") or config.exists("archive-alternategrammars5"))
        {
          alt_gs5 = build_alternate(all_annot_descendants5,4);
        }

        ///////////////////
        results.push_back(create_parser(grammars5, string_to_pa(config.get_value<std::string>("parser-type")),
                                        string_to_mpc(config.get_value<std::string>("max-type")), priors5, beam_threshold,
                                        alt_gs5, all_annot_descendants5, accurate, min, config.get_value<int>("stubbornness"),
                                        config.get_value<unsigned>("kbest")));

        if (grammars5.back()->gram_conf.count("lex_unknown_map"))
        {
          if (results.back()->get_word_signature() == nullptr)
          {
          std::cerr << "overwriting unknown_map from command-line (if you don't want this, edit the grammar)" << std::endl;
          results.back()->set_word_signature(
              WordSignatureFactory::create_wordsignature(
                  WordSignature::string_2_lex_unknown_map(grammars5.back()->gram_conf.at("lex_unknown_map")),
                  true));
          }
        }
        if (grammars5.back()->gram_conf.count("remove_functions"))
        {
          if(grammars5.back()->gram_conf.at("remove_functions") == "0")
          {
            results.back()->set_is_funct(true);
          }
        }
      }
    }


    ////////////////
    {
      std::vector<double> priors6;
      std::vector<annot_descendants_type> all_annot_descendants6;
      std::vector<ParserCKYAll::AGrammar*> grammars6;
      std::vector< std::vector<ParserCKYAll::AGrammar*> > alt_gs6;
      // get grammars
      if(config.exists("grammar6") or config.exists("archive-grammar6")) {

        grammars6 = build_grammar(all_annot_descendants6, 5);
        // compute priors for base grammar
        priors6 = grammars6[0]->compute_priors();

        ///////////////////
        if(config.exists("alternate-grammar6") or config.exists("archive-alternategrammars6"))
        {
          alt_gs6 = build_alternate(all_annot_descendants6,3);
        }

        ///////////////////
        results.push_back(create_parser(grammars6, string_to_pa(config.get_value<std::string>("parser-type")),
                                        string_to_mpc(config.get_value<std::string>("max-type")), priors6, beam_threshold,
                                        alt_gs6, all_annot_descendants6, accurate, min, config.get_value<int>("stubbornness"),
                                        config.get_value<unsigned>("kbest")));

        if (grammars6.back()->gram_conf.count("lex_unknown_map"))
        {
          if (results.back()->get_word_signature() == nullptr)
          {
          std::cerr << "overwriting unknown_map from command-line (if you don't want this, edit the grammar)" << std::endl;
          results.back()->set_word_signature(
              WordSignatureFactory::create_wordsignature(
                  WordSignature::string_2_lex_unknown_map(grammars6.back()->gram_conf.at("lex_unknown_map")),
                  true));
          }
        }
        if (grammars6.back()->gram_conf.count("remove_functions"))
        {
          if(grammars6.back()->gram_conf.at("remove_functions") == "0")
          {
            results.back()->set_is_funct(true);
          }
        }
      }
    }




    return results;
}



///////////
// build the mapping 'gram_idx->nt_symb->annot->annots in gram_idx + 1'
// needed for c2f pruning
///////////
/**
 *     \brief create an efficient data structure for c2f parsing
 *   from a vector of trees denoting annotation histories
 *   \param annot_histories the trees
 */
annot_descendants_type
    create_annot_descendants(const std::map<short, Tree<unsigned> >& annot_histories)
{
  annot_descendants_type result;

  // how many intermediate grammars should we create ?
  unsigned height = 0;
  for (auto& e : annot_histories)
  {
    if(e.second.height() > height)
      height = e.second.height();
  }
  //    std::clog << height << std::endl;
  result.resize(height);

  // how many symbols (some symbols may not be in the grammar)
  unsigned nb_symbols = SymbolTable::instance_nt().get_symbol_count();

  for (unsigned nt_idx = 0; nt_idx < nb_symbols; ++nt_idx)
  {
    //      std::clog << annot_histories[nt_idx] << std::endl;

    for(unsigned gram_idx = 0; gram_idx < height; ++gram_idx)
    {
      //    std::clog << gram_idx << " : " << std::endl;
      result[gram_idx].resize(nb_symbols);

      if (annot_histories.count(nt_idx))
      {
         std::vector<Tree<unsigned>::const_depth_first_iterator> subs = annot_histories.at(nt_idx).get_at_depth(gram_idx);
         result[gram_idx][nt_idx].resize(subs.size());

        for(std::vector<Tree<unsigned>::const_depth_first_iterator>::const_iterator i(subs.begin()); i != subs.end(); ++i)
        {

          std::vector<unsigned> descs;
          Tree<unsigned>::const_depth_first_iterator copy = *i;
          copy.down_first();
          while(copy != annot_histories.at(nt_idx).dfend())
          {
            //    std::clog << '\t' << gram_idx << "->" << nt_idx << "->" << **i << "->" << *copy << std::endl;
            descs.push_back(*copy);
            copy.right();
          }
          result[gram_idx][nt_idx][**i]=descs;
        }
      }
    }
  }

  return result;
}

std::vector<ParserCKYAll::AGrammar*>
create_intermediates(ParserCKYAll::AGrammar& grammar, const annot_descendants_type& annot_descendants)
{
  std::vector<ParserCKYAll::AGrammar*> result(annot_descendants.size() -1);

  //std::clog << "before transition" << std::endl;
  uomap<int, uomap<unsigned, uomap<int, uomap<unsigned, double > > > >transition_probabilities;
  grammar.compute_transition_probabilities(transition_probabilities);
  //  std::clog << "after transition" << std::endl;

  //std::clog << "before expected counts" << std::endl;
  std::vector<std::vector<double> > expected_counts;
  calculate_expected_counts(transition_probabilities, grammar.get_annotations_info(), expected_counts);
  //  std::clog << "after expected counts" << std::endl;
  // for(size_t i = 0; i < expected_counts.size(); ++i)
  //   for(size_t j = 0; j < expected_counts[i].size(); ++j)
  //     std::cout << expected_counts[i][j] << std::endl;

#ifndef USE_THREADS
  for (unsigned i = 0; i < annot_descendants.size() - 1; ++i) {

    //      std::clog << "before mapping " << i << std::endl;
    std::vector<std::vector<std::vector<unsigned> > > annotation_mapping = compute_mapping(i, annot_descendants.size() - 1, annot_descendants);
    //      std::clog << "after mapping " << i << std::endl;

    //      std::clog << "before create_projection" << std::endl;
    ParserCKYAll::AGrammar * cg = grammar.create_projection(expected_counts, annotation_mapping);
    //      std::clog << "after create_projection" << std::endl;


    //      std::clog << "finishing creation" << std::endl;
    //add unary chains
    //      std::clog << "init" << std::endl;
    cg->init();
    //smooth
    //      std::clog << "smoothing" << std::endl;
    cg->linear_smooth(0.01,0.1);
    // remove zeros
    //      std::clog << "removing zeros" << std::endl;
    cg->remove_unlikely_annotations_all_rules(1e-10);
    // done!
    result[i] = cg;
    //      std::clog << "creation finished" << std::endl;
  }
#else
  tbb::parallel_for(tbb::blocked_range<unsigned>(0, annot_descendants.size() -1),
                    [&]
                    (tbb::blocked_range<unsigned>& r)
                    {
                      for(unsigned i = r.begin(); i < r.end(); ++i) {
                        //      std::clog << "before mapping " << i << std::endl;
                        std::vector<std::vector<std::vector<unsigned> > > annotation_mapping = compute_mapping(i, annot_descendants.size() - 1, annot_descendants);
                        //      std::clog << "after mapping " << i << std::endl;

                        //      std::clog << "before create_projection" << std::endl;
                        // std::clog << "create_intermediates i=" << i << std::endl;
                        ParserCKYAll::AGrammar * cg = grammar.create_projection(expected_counts, annotation_mapping);
                        //      std::clog << "after create_projection" << std::endl;


                        //      std::clog << "finishing creation" << std::endl;
                        //add unary chains
                        //      std::clog << "init" << std::endl;
                        cg->init();
                        //smooth
                        //      std::clog << "smoothing" << std::endl;
                        cg->linear_smooth(0.01,0.1);
                        // remove zeros
                        //      std::clog << "removing zeros" << std::endl;
                        cg->remove_unlikely_annotations_all_rules(1e-10);
                        // done!
                        result[i] = cg;
                        //      std::clog << "creation finished" <<
                        //      std::endl;
                      }
                    }
                    );
#endif

  return result;
}
