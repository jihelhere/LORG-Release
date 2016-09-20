// -*- mode: c++ -*-

#pragma once
#include <boost/filesystem.hpp>


#include "Treebank.h"
#include "utils/ConfigTable.h"

#include "utils/PtbPsTree.h"

#include <boost/tokenizer.hpp>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#include <boost/phoenix.hpp>
#pragma clang diagnostic pop




namespace {
    namespace px = boost::phoenix;
    namespace px_args = px::arg_names;
    namespace fs = boost::filesystem;

    //helper function: gets all the files in path, recursively
    void collect_files(const fs::path& path, std::vector<std::string>& files)
    {
        if(fs::exists(path)) {
            if(fs::is_directory(path)) {
                fs::directory_iterator end_iter;
                for (fs::directory_iterator dir_itr(path); dir_itr != end_iter; ++dir_itr )
                    collect_files(*dir_itr,files);
            }
            else files.push_back(path.string());
        }
    }

    //helper function: gets all the files in a collection of paths
    void collect_files_from_vector(const std::vector<std::string>& vector_names, std::vector<std::string>& filenames,bool verbose)
    {
        for(std::vector<std::string>::const_iterator iter = vector_names.begin(); iter != vector_names.end(); ++iter) {
            if(verbose) {std::clog << "\"" << *iter << "\" ";}
            fs::path original_path = fs::system_complete(fs::path(iter->c_str()));
            collect_files(original_path,filenames);
        }
        if(verbose) {std::clog << "\n";}
    }
}


namespace TreebankFactory {

  template<class Tree>
  Treebank<Tree> BuildEmptyTreebank(const ConfigTable& config);

  treebank_options read_config(const ConfigTable& conf);
}


treebank_options TreebankFactory::read_config(const ConfigTable& conf)
{
  bool verbose = conf.exists("verbose");

  std::vector<std::string> treebank_args = conf.get_value< std::vector<std::string> >("treebank");
  std::vector<std::string> treebank_filenames;
  if(verbose)
    std::clog << "treebank set to: ";
  collect_files_from_vector(treebank_args,  treebank_filenames, verbose);


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
  return
  treebank_options(treebank_filenames,
                   labels_to_remove,
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
}


template<class Tree>
Treebank<Tree> TreebankFactory::BuildEmptyTreebank(const ConfigTable& conf)
{
  bool verbose = conf.exists("verbose");
  return Treebank<Tree>(read_config(conf), verbose);
}
