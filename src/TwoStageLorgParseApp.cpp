// -*- mode: c++ -*-
#include "TwoStageLorgParseApp.h"

#include "utils/Tagger.h"
#include "feature_extract/Extract.h"
#include "utils/PtbPsTree.h"
#include "utils/LorgConstants.h"
#include "parsers/ParserCKYAllFactory.h"

#include "parsers/ParserCKYAll.h"

#include "utils/tick_count.h"


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

        tick_count sent_start = tick_count::now();

        if(sentence.size() <=  max_length && sentence.size() > 0) {

            /*Extra verbose
            if(verbose) {
                std::clog << "Tokens: ";
                for(std::vector<Word>::const_iterator i(sentence.begin()); i != sentence.end(); ++i)
                    std::clog << "<" << i->get_form() << ">";
                std::clog << "\n";
            }*/

          for (size_t i = 0; i < parsers.size(); ++i)
          {
            std::vector<std::pair<PtbPsTree *,double> > best_trees; // vector of (tree,score)

            //std::cout << "iteration " << i << std::endl;
            //tag sentence
            {
              //                             BLOCKTIMING("tagger");
              //std::cout << "tagging " << i << std::endl;
              taggers[i].tag(sentence);
            }
            // create and initialise chart
            {
              //                             BLOCKTIMING("initialise_chart");
              //std::cout << "intialisation " << i << std::endl;
              parsers[i]->initialise_chart(sentence, brackets);
            }


            // parse, aka create the coarse forest
            {
              //                             BLOCKTIMING("parse");
              //std::cout << "parsing " << i << std::endl;
              parsers[i]->parse(start_symbol);
            }

            //use intermediate grammars to prune the chart
            {
              //                             BLOCKTIMING("beam_c2f");
              //std::cout << "beaming " << i << std::endl;
              parsers[i]->beam_c2f(start_symbol);
            }

            // extract best solution with the finest grammar
            if(parsers[i]->is_chart_valid(start_symbol))
            {
              //                             BLOCKTIMING("extract_solution");

              //std::cout << "extracting " << i << std::endl;
              parsers[i]->extract_solution();
            }

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
            for(unsigned i = 0; i < best_trees.size(); ++i) { // delete solutions
              delete best_trees[i].first;
            }
            best_trees.clear();
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

    parsers[0]->set_nbthreads(this->nbthreads);

    if(verbose) {std::clog << "ok" << std::endl;}

    kbest = configuration.get_value<unsigned>("kbest");

    //creating taggers
    taggers.resize(parsers.size());
    for (size_t i = 0; i < taggers.size(); ++i)
    {
      taggers[i].set_word_rules(&(parsers[i]->get_words_to_rules()));
    }


    extract_features = configuration.get_value<bool>("extract-features");

    output_format = parse_solution::format_from_string(configuration.get_value<std::string>("output-format"));

    return true;
}
