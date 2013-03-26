// -*- mode: c++ -*-
#ifndef _PARSERCKYALLFACTORY_CPP_
#define _PARSERCKYALLFACTORY_CPP_

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
create_annot_descendants(const std::vector< Tree<unsigned> >& annot_histories);


std::vector<ParserCKYAll::AGrammar*> create_grammars(const std::string& filename, bool verbose)
{
  std::vector<ParserCKYAll::AGrammar*> grammars;

  std::vector<annot_descendants_type> all_annot_descendants;

  // compute expected intermediate grammars

  //get grammar
  if(verbose) std::cerr << "Setting grammar to " << filename << ".\n";
  ParserCKYAll::AGrammar * cg = new ParserCKYAll::AGrammar(filename);
  if(verbose) std::cerr << "Grammar set\n";

  //perform some sanity checks
  if(cg->get_history_trees().empty() ||
     cg->get_annotations_info().get_number_of_unannotated_labels() != cg->get_history_trees().size()) {

    std::cerr << "Problem in the grammar file." << std::endl
              << "Annotation history and number of annotations are inconsistent." << std::endl
              << "Aborting now !" << std::endl;
    delete cg;
    exit(1);
  }

  if(verbose)
    std::clog << "create intermediate grammars" << std::endl;

  annot_descendants_type annot_descendants = create_annot_descendants(cg->get_history_trees());
  grammars = create_intermediates(*cg, annot_descendants);

  // compute priors for base grammar
  //        std::clog << "before priors" << std::endl;
  //  priors = grammars[0]->compute_priors();
  cg->init();
  //    std::clog << "smooth" << std::endl;
  // TODO: options to set this from command line
  cg->linear_smooth(0.01,0.1);
  //    std::clog << "clean" << std::endl;
  cg->remove_unlikely_annotations_all_rules(1e-10);
  grammars.push_back(cg);
  //    std::clog << "creation finished" << std::endl;

  return grammars;
}

std::vector<ParserCKYAll *>
ParserCKYAllFactory::create_parser(ConfigTable& config)
{
    bool verbose = config.exists("verbose");

    #ifdef USE_THREADS
    unsigned nbthreads = config.get_value<unsigned>("nbthreads");

    #endif
    if(verbose){
      std::clog << "using " <<
      #ifndef USE_THREADS
      1
      #else
      (nbthreads == 0 ? tbb::task_scheduler_init::default_num_threads() : nbthreads)
      #endif
      << " thread(s) to parse" << std::endl;
    }

    std::vector<ParserCKYAll*> results;

    std::vector<double> priors;
    std::vector<annot_descendants_type> all_annot_descendants;

    std::vector<ParserCKYAll::AGrammar*> grammars;
    // get grammars
    if(config.exists("grammar")) {
        const std::string& filename = config.get_value< std::string >("grammar");

        grammars = create_grammars(filename, verbose);
        // compute priors for base grammar
        //        std::clog << "before priors" << std::endl;
        priors = grammars[0]->compute_priors();
    }
    else {
        std::cerr << "Grammar wasn't set. Exit program." << std::endl;
        return results;
    }


    annot_descendants_type annot_descendants = create_annot_descendants(grammars.back()->get_history_trees());
    all_annot_descendants.push_back(annot_descendants);


    std::cout << "here" << std::endl;

    std::vector< std::vector<ParserCKYAll::AGrammar*> > alt_gs;
    if(config.exists("alternate-grammar")) {

        const std::vector<std::string>& filenames = config.get_value<std::vector<std::string> >("alternate-grammar");
        if(string_to_pa(config.get_value<std::string>("parser-type")) != MaxN && filenames.size() > 0) {
            std::cerr << "Wrong parsing algorithm. Exit program." << std::endl;
            return results;
        }

        for(unsigned i = 0; i < filenames.size(); ++i)
        {
          if(verbose) std::cerr << "Setting alternate grammar to " << filenames[i] << ".\n";
          alt_gs.push_back(create_grammars(filenames[i], verbose));


          annot_descendants_type annot_descendants = create_annot_descendants(alt_gs.back().back()->get_history_trees());
          all_annot_descendants.push_back(annot_descendants);
        }
    }

    double beam_threshold = config.get_value<double>("beam-threshold");
    if(verbose)
        std::clog << "using threshold " << beam_threshold << " for chart construction" << std::endl;

    bool accurate = config.exists("accurate");

    if(verbose) {
        if(accurate)
            std::clog << "using accurate c2f thresholds" << std::endl;
        else
            std::clog << "using standard c2f thresholds" << std::endl;
    }


    unsigned min = config.get_value<unsigned>("min-length-beam");

    results.push_back(
      create_parser(grammars,
        string_to_pa(config.get_value<std::string>("parser-type")),
        string_to_mpc(config.get_value<std::string>("max-type")),
        priors, beam_threshold,
        alt_gs, all_annot_descendants, accurate, min,
        config.get_value<int>("stubbornness"),
        config.get_value<unsigned>("kbest"))
    );



    std::vector<double> priors2;
    std::vector<annot_descendants_type> all_annot_descendants2;
    std::vector<ParserCKYAll::AGrammar*> grammars2;
    std::vector< std::vector<ParserCKYAll::AGrammar*> > alt_gs2;
    // get grammars
    if(config.exists("grammar2")) {
        const std::string& filename = config.get_value< std::string >("grammar");

        grammars2 = create_grammars(filename, verbose);
        // compute priors for base grammar
        //        std::clog << "before priors" << std::endl;
        priors2 = grammars2[0]->compute_priors();

        annot_descendants_type annot_descendants = create_annot_descendants(grammars.back()->get_history_trees());
        all_annot_descendants2.push_back(annot_descendants);

        results.push_back(
          create_parser(grammars2,
            string_to_pa(config.get_value<std::string>("parser-type")),
            string_to_mpc(config.get_value<std::string>("max-type")),
            priors2, beam_threshold,
            alt_gs2, all_annot_descendants2, accurate, min,
            config.get_value<int>("stubbornness"),
            config.get_value<unsigned>("kbest"))
        );

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
create_annot_descendants(const std::vector< Tree<unsigned> >& annot_histories)
{
  annot_descendants_type result;

  unsigned height = 0;
  for (unsigned i = 0; i < annot_histories.size(); ++i)
  {
    if(annot_histories[i].height() > height)
      height = annot_histories[i].height();
  }
  //    std::clog << height << std::endl;

  result.resize(height);

  for (unsigned nt_idx = 0; nt_idx < annot_histories.size(); ++nt_idx)
  {
    //      std::clog << annot_histories[nt_idx] << std::endl;

    for(unsigned gram_idx = 0; gram_idx < height; ++gram_idx)
    {
      //    std::clog << gram_idx << " : " << std::endl;
      result[gram_idx].resize(annot_histories.size());

      std::vector<Tree<unsigned>::const_depth_first_iterator> subs = annot_histories[nt_idx].get_at_depth(gram_idx);
      result[gram_idx][nt_idx].resize(subs.size());

      for(std::vector<Tree<unsigned>::const_depth_first_iterator>::const_iterator i(subs.begin()); i != subs.end(); ++i) {

        std::vector<unsigned> descs;
        Tree<unsigned>::const_depth_first_iterator copy = *i;
        copy.down_first();
        while(copy != annot_histories[nt_idx].dfend())
        {
          //    std::clog << '\t' << gram_idx << "->" << nt_idx << "->" << **i << "->" << *copy << std::endl;
          descs.push_back(*copy);
          copy.right();
        }
        result[gram_idx][nt_idx][**i]=descs;
      }
    }
  }

  return result;
}

std::vector<ParserCKYAll::AGrammar*>
create_intermediates(ParserCKYAll::AGrammar& grammar, const annot_descendants_type& annot_descendants)
{
  std::vector<ParserCKYAll::AGrammar*> result(annot_descendants.size() -1);

  //  std::clog << "before transition" << std::endl;
  uomap<int, uomap<unsigned, uomap<int, uomap<unsigned, double > > > >transition_probabilities;
  grammar.compute_transition_probabilities(transition_probabilities);
  //  std::clog << "after transition" << std::endl;

  //  std::clog << "before expected counts" << std::endl;
  std::vector<std::vector<double> > expected_counts;
  calculate_expected_counts(transition_probabilities, grammar.get_annotations_info(), expected_counts);
  //  std::clog << "after expected counts" << std::endl;

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
    cg->remove_unlikely_annotations_all_rules(1e-6);
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
                        cg->remove_unlikely_annotations_all_rules(1e-6);
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

#endif /* _PARSERCKYALLFACTORY_CPP_ */
