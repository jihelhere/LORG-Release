#include "LorgTrainerApp.h"

#include <iostream>

//#define BOOST_FILESYSTEM_VERSION 2
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "TreebankFactory.h"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/phoenix/statement/sequence.hpp>
#include <boost/phoenix/object/construct.hpp>
#include <boost/phoenix/stl.hpp>

#include <boost/phoenix/bind/bind_member_function.hpp>
#include <boost/phoenix/bind/bind_member_variable.hpp>
#pragma clang diagnostic pop

#include "utils/RandomGenerator.h"

#include "utils/hash_impl.h"

#include "EMTrainer.h"

#include "edges/BinaryTrainingTree.h"
#include "utils/RuleVect2Map.h"

#include <tbb/tick_count.h>
using namespace tbb;


// HELPER FUNCTIONS IN ANONYMOUS NAMESPACE
namespace {
    void ptbpstrees_to_bttrees(const std::vector<PtbPsTree>& input,
                               TrainingGrammar& grammar, std::vector<BinaryTrainingTree>& result)
    {
        MAP<std::pair<int,std::pair<int,int> >,BRuleTraining*> brulemap;
        MAP<std::pair<int,int>,URuleTraining*> urulemap;
        MAP<std::pair<int,int>,LexicalRuleTraining*> lrulemap;

        //creating the maps from the vectors
        rulevect2map<std::pair<int,std::pair<int,int> >,BRuleTraining> vect2bmap(brulemap, grammar.get_binary_rules());

        rulevect2map<std::pair<int,int>,URuleTraining> vect2umap(urulemap, grammar.get_unary_rules());

        rulevect2map<std::pair<int,int>,LexicalRuleTraining> vect2lmap(lrulemap, grammar.get_lexical_rules());

        //building training trees from ptbpstrees and maps
        for(const auto& p : input)
        {
          result.emplace_back(p,brulemap,urulemap,lrulemap);
        }




    }
}

LorgTrainerApp::LorgTrainerApp()
    : lexicon(nullptr),
      filter_level(0), baseoutputname(), split_size(0), n_iterations(0),
      merge_em(0), split_em(0), turn_off_merge(false), base_grammar_only(false),
      smooth_grammar(0), smooth_lexicon(0),
      split_randomness(0), merge_percentage(0), prob_threshold(0.0),
      modified_treebank_name(), keep_fun(false)
{}

LorgTrainerApp::~LorgTrainerApp()
{}



// create a zero-length file to test output
void LorgTrainerApp::test_grammar(const std::string & suffix, int num)
{
    if(baseoutputname != "")
    {
        std::ostringstream op;
        op << baseoutputname << suffix ;
        if(num >=0 ) op << '_' << num;

        std::string opath_string = op.str();
        fs::path output_path = fs::system_complete(fs::path(opath_string.c_str()));
        fs::ofstream file_out(output_path,std::ios_base::out | std::ios_base::trunc);

        if(!file_out) {
            std::cerr << "Could not open file " << opath_string << std::endl;
            //is it clean ?
            exit(1);
        }
        file_out.close();
        fs::remove(output_path);
    }
}



//create a new file baseoutputname+suffix, and open it for writing
void LorgTrainerApp::write_grammar(const TrainingGrammar& grammar, const std::string & suffix, int num)
{
    if(this->baseoutputname != "")
    {
        std::ostringstream op;
        op << baseoutputname << suffix ;
        if(num >=0 ) op << '_' << num;

        std::string opath_string = op.str();
        fs::path output_path = fs::system_complete(fs::path(opath_string.c_str()));
        fs::ofstream file_out(output_path,std::ios_base::out | std::ios_base::trunc);

        if(!file_out) {
            std::cerr << "Could not open file " << opath_string << std::endl;
        }


        file_out << grammar.get_lexicon()->header_string() ;

        // TODO find a proper way to do that
        if (this->keep_fun)
        {
          file_out << std::string("// conf: remove_functions 0") << std::endl;
        }
        else
        {
          file_out << std::string("// conf: remove_functions 1") << std::endl;
        }

        file_out << std::endl;

        file_out << grammar << std::flush;
        file_out.close();
    }
    else
        std::cout << grammar << std::flush;
}



int LorgTrainerApp::run()
{
  // initialise treebank from configuration
  //initialise training grammar from treebank

  Treebank<PtbPsTree> tb(tb_options, verbose);
  if(verbose)
    std::clog << "Treebank is ready ("
              << tb.get_trees().size()
              << " trees)"  << "\n";

  TrainingGrammar em_grammar(tb, lexicon);

  //output treebank trees modified by the lexicon -- if requested
  if(modified_treebank_name != "")
  {
    fs::ofstream file_out(modified_treebank_name, std::ios_base::out | std::ios_base::trunc);
    if(!file_out)
    {
      std::cerr << "Could not open file " << modified_treebank_name << std::endl;
    }

    tb.output_unbinarised(file_out);
    file_out.close();
  }

  //output base grammar
  if (verbose)
  {
    test_grammar("_base");

    std::clog << "outputting baseline grammar " << std::endl;
    //    write_grammar(em_grammar,"_base");
    {
      TrainingGrammar gout;
      assert(em_grammar.get_lexicon() != NULL);
      gout = em_grammar;
      assert(em_grammar.get_lexicon() != NULL);

      //      std::cout << "creating additional rules" << std::endl;
      gout.create_additional_rules();
      //      std::cout << "created " << gout.additional_rules.size() << " additional rules" << std::endl;

      gout.lexical_smoothing();
      write_grammar(gout,"_base");
      assert(gout.get_lexicon() != em_grammar.get_lexicon());
    }
  }


  if (base_grammar_only){ return 0;}


  if (verbose)
  {
    for(unsigned iter = 1; iter <= n_iterations; ++iter)
      test_grammar("_smoothed",iter);
  }
  test_grammar("_final",-1);


  // SPLIT/EM starts here !

  //transform ptb trees into edges for faster access
  std::clog << "converting ptbpstrees to binary training trees" << std::endl;

  // the collection representing the trees for inside/outside computations
  std::vector<BinaryTrainingTree> training_trees;
  ptbpstrees_to_bttrees(tb.get_trees(), em_grammar, training_trees);

  //hack for now -- needs to be rewritten
  keep_fun = not(tb.get_options().func);

  // create the trainer object
  EMTrainer em_trainer(split_size,prob_threshold,verbose);

  if(verbose) {
    std::clog.precision(22);
    std::clog << "\n" <<"baseline training set log-likelihood is: "
              << em_trainer.calculate_likelihood(training_trees) << std::endl;
  }

  unsigned size = 1;
  tick_count t0 = tick_count::now();
  for(unsigned iter = 1; iter <= n_iterations; ++iter) {
    //tick_count itersm_start = tick_count::now();
    std::clog << "Split/Merge cycle number " << iter << std::endl;


    if(verbose){
      size *= split_size;
      std::clog << "maximum training tree node  annotation size: " << size << std::endl;
      std::clog << "splitting rules annotations in " << split_size << std::endl;
    }

    ///////////////////////////////////////
    //split edges
        ////////////////////////////////////////

    em_trainer.split(training_trees, em_grammar, split_randomness);

    // if(verbose) {
    //   std::clog << "outputting split grammar " << std::endl;
    //   write_grammar(em_grammar,"_split",iter);

    // }

    if(verbose)  std::clog << "Starting EM iterations for SPLIT GRAMMAR " << iter << std::endl;

    em_trainer.do_em(training_trees, em_grammar, split_em, 0.0, 0.0, true, smooth_method);

    assert(em_grammar.get_lexicon() != NULL);


    // if(verbose) {
    //   std::clog << "outputting grammar after split +EM" << std::endl;
    //   write_grammar(em_grammar,"_splitEM",iter);
    //   {
    // 	assert(em_grammar.get_lexicon() != NULL);
    // 	TrainingGrammar gout;
    // 	gout = em_grammar;
    // 	assert(em_grammar.get_lexicon() != NULL);
    // 	gout.lexical_smoothing();
    // 	em_trainer.smooth_grammar_rules(gout, 0, smooth_lexicon, smooth_method);
    // 	write_grammar(gout,"_splitEM_ls",iter);
    // 	assert(gout.get_lexicon() != em_grammar.get_lexicon());
    //   }
    // }

    ///////////////////////////////
    //Merging step
    ////////////////////////////////
    if (!turn_off_merge && merge_percentage > 0) {


      assert(em_grammar.get_lexicon() != NULL);


      std::clog << "Merging grammar..." << std::endl;
      em_trainer.merge(training_trees, em_grammar, merge_percentage, final_lex_em);

      std::clog << "Starting EM iterations for MERGED GRAMMAR " << iter <<  std::endl;
      em_trainer.do_em(training_trees, em_grammar,
                       merge_em, 0.0, 0.0, true, smooth_method);

      // if(verbose) {
      // 	std::clog << "outputting grammar after merge + EM" << std::endl;
      // 	write_grammar(em_grammar, "_mergedEM", iter);
      // 	{
      // 	  TrainingGrammar gout;
      // 	  gout = em_grammar;
      // 	  gout.lexical_smoothing();
      // 	  em_trainer.smooth_grammar_rules(gout, 0, smooth_lexicon, smooth_method);
      // 	  write_grammar(gout,"_mergedEM_ls",iter);
      // 	}
      // }
    }



    ///////////////////////////////
    //Smoothing step
    ////////////////////////////////
    if(smooth_grammar > 0.0 || smooth_lexicon > 0.0) {

      if(verbose) std::clog << "Smoothing grammar" << std::endl;

      if(verbose) std::clog << "Starting EM iterations for SMOOTHED GRAMMAR " << iter << std::endl;

      em_trainer.smooth_grammar_rules(em_grammar, 0, smooth_lexicon, smooth_method);
      //	std::cout << em_trainer.calculate_likelihood(training_trees) << std::endl;

      em_trainer.do_em(training_trees, em_grammar, 10,
                       smooth_grammar, smooth_lexicon, true, smooth_method);

      if(verbose){
        std::clog << "outputting grammar after smooth + EM" << std::endl;


        if(lexicon_type == LexiconFactory::Basic)
          write_grammar(em_grammar, "_smoothed", iter);
        else
        {
          TrainingGrammar gout;
          gout = em_grammar;

          //std::cout << "creating additional rules" << std::endl;
          gout.create_additional_rules();
          //std::cout << "created " << gout.additional_rules.size() << " additional rules" << std::endl;

          gout.lexical_smoothing();
          //	  em_trainer.smooth_grammar_rules(gout, smooth_grammar, smooth_lexicon, smooth_method);
          write_grammar(gout,"_smoothed",iter);
        }
      }
    }
    //std::cout << "### itersm:\t" << iter << "\t" << (tick_count::now() - itersm_start).seconds() << std::endl;
  } //End of for split-merge cycles


  if (lexicon_type != LexiconFactory::Basic)
  {
    //one iteration with lexical rule smoothing turned on, relevant for berkeley sophisticated lexicon
    em_trainer.do_em(training_trees, em_grammar, 2, 0.0, 0.0, true, smooth_method);

    //  std::cout << "creating additional rules" << std::endl;
    em_grammar.create_additional_rules();
    //  std::cout << "created " << em_grammar.additional_rules.size() << " additional rules" << std::endl;

    em_grammar.lexical_smoothing();
    em_trainer.smooth_grammar_rules(em_grammar, smooth_grammar, smooth_lexicon,smooth_method);
  }
  std::clog << "outputting final grammar" << std::endl;
  write_grammar(em_grammar, "_final");


  for (auto& t: training_trees)
    t.free_nodes();

  if(verbose)
    std::cerr << "overall time: " << (tick_count::now() - t0).seconds() << std::endl;

  return 0;
}


LorgOptions LorgTrainerApp::get_options() const
{
    LorgOptions options;
    options.add_grammar_options();
    options.add_treebank_processing_options();
    options.add_lexicon_options();
    options.add_unknown_word_mapping_options();
    options.add_grammar_positionals();
    return options;
}

bool LorgTrainerApp::read_config(ConfigTable& conf)
{
    if(LorgApp::read_config(conf) == false)
        return false;

    filter_level = conf.get_value<int>("filter-level");
    if (verbose) std::clog << "Filter level set to " << filter_level << "\n";


    // output grammar
    if(conf.exists("output")) { baseoutputname = conf.get_value<std::string>("output");}

    split_size = conf.get_value<unsigned>("split-size");
    n_iterations = conf.get_value<unsigned>("split-merge-cycles");

    merge_em = conf.get_value<unsigned>("merge-em");
    split_em = conf.get_value<unsigned>("split-em");

    turn_off_merge = conf.get_value<bool>("turn-off-merge");
    merge_percentage = conf.get_value<unsigned>("merge-percentage");
    if(merge_percentage > 100 || merge_percentage < 0) {
        std::clog << "A percentage of " << merge_percentage  << " ! Are you insane ?" << std::endl;
        return false;
    }
    smooth_grammar = conf.get_value<double>("smooth-grammar");
    smooth_lexicon = conf.get_value<double>("smooth-lexicon");

    split_randomness = conf.get_value<unsigned>("split-random");
    prob_threshold = conf.get_value<double>("prob-threshold");
    if(verbose) std::clog << "Annoted rules with probability below "
        << prob_threshold << " will be removed." << std::endl;
    base_grammar_only =  conf.get_value<bool>("base-grammar-only");

    lexical_smoothing = conf.exists("lexical_smoothing");

    unsigned seed = conf.get_value<unsigned>("random-seed");
    if(verbose) std::clog << "random seed is " << seed << std::endl;
    RandomGenerator::set_global_seed(seed);

    //sub modules can read configuration too

    //tb = TreebankFactory::BuildEmptyTreebank<PtbPsTree>(conf);
    if (conf.exists("treebank")) {
      tb_options = TreebankFactory::read_config(conf);
    }
    else {
        std::cerr << "treebank was not set. Exiting...\n";
        return false;
    }

    lexicon = std::shared_ptr<Lexicon>(LexiconFactory::create_lexicon(conf));

    lexicon_type = LexiconFactory::string_2_lex_type(conf.get_value<std::string>("lexicon-type"));


    modified_treebank_name = conf.get_value<std::string>("modified-treebank");


    //  type of smoothing

    std::string conf_smooth(conf.get_value<std::string>("smooth-method"));
    if(conf_smooth != "linear" && conf_smooth != "weighted") {
        std::clog << "smooth-method is not known: " << conf_smooth << std::endl;
        return false;
    }

    smooth_method = (conf_smooth == "linear") ? TrainingGrammar::LinearSmooth : TrainingGrammar::WeightedSmooth;

    final_lex_em = conf.get_value<bool>("final-lex-em");

    return true;
}
