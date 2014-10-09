// -*- mode: c++ -*-

#ifndef TREEBANKFACTORY_H
#define TREEBANKFACTORY_H


#include "Treebank.h"
#include "utils/ConfigTable.h"

#include "utils/PtbPsTree.h"

#include <boost/tokenizer.hpp>


namespace TreebankFactory {

  template<class Tree>
  Treebank<Tree> * BuildEmptyTreebank(ConfigTable& config);
}


template<class Tree>
Treebank<Tree> * TreebankFactory::BuildEmptyTreebank(ConfigTable& conf)
{
  bool verbose = conf.exists("verbose");

  Bin_Direction bindir = LEFT;

  std::string bin =  conf.get_value<std::string>("binarisation");
  if (verbose)
    std::clog << "binarisation mode set to " << bin << ".\n";

  if(bin == "none")
  {
    bindir = NONE;
    if (verbose)
      std::clog << "Grammar will not be binarised\n";
  }
  else if (bin == "right")
  {
    bindir = RIGHT;
    if (verbose)
      std::clog << "Grammar will be right-binarised\n";
  }
  else if (bin == "left")
  {
    bindir = LEFT;
    if (verbose)
      std::clog << "Grammar will be left-binarised\n";
  }
  else if (bin == "random")
  {
    bindir = RAND;
    if (verbose)
      std::clog << "Grammar will be randomly binarised\n";
  }
  else
  {
    throw std::out_of_range("binarisation direction");
  }

  HorizMarkov binmark =  -1 ; //infinite -> "exact" binarization
  binmark =  conf.get_value<HorizMarkov>("hm");


  if(verbose) {
    if(binmark < 0)
      std::clog << "Horizontal markovisation set to infinite\n";
    else
      std::clog << "Horizontal markovisation set to " <<  binmark << "\n";
  }

  //reading labels to remove from trees
  std::unordered_set<std::string> labels_to_remove;
  std::string labels = conf.get_value<std::string>("nodes-to-remove");
  boost::tokenizer<boost::char_separator<char> > tokeniser(labels,boost::char_separator<char>(" "));
  for(const auto& c : tokeniser)
  {
    if(verbose) std::clog << "Will remove label: " << c << " from trees\n";
    labels_to_remove.insert(c);
  }

  //reading regex to recognize numbers
  boost::regex num_regex(conf.get_value<std::string>("number-regex"));

  // read treebank processing options
  treebank_options tb_options(labels_to_remove,
                              conf.get_value<bool>("remove-functions"),
                              conf.get_value<bool>("replace-numbers"),
                              num_regex,
                              conf.get_value<unsigned>("vm")== 0 ? 0 :conf.get_value<unsigned>("vm") -1,
                              //difference vm/additional layer : vm = 2 correponds to adding 1 layer
                              conf.get_value<bool>("pos-vm"),
                              conf.get_value<bool>("remove-same-unary"),
                              bindir,
                              binmark,
                              conf.get_value<unsigned>("sentence-max-length"));

  return  new Treebank<Tree>(tb_options,verbose);


}



#endif // TREEBANKFACTORY_H
