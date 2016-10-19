// -*- mode: c++ -*-
#pragma once
#include "LorgParseApp.h"

#include "parsers/ParserCKYNN.h"


#include "training/Treebank.h"

#include "lexicon/WordSignature.h"

class NNLorgParseApp : public LorgParseApp
{
public:
  NNLorgParseApp();
  ~NNLorgParseApp();
  int run();

private:
  bool read_config(ConfigTable& configuration);
  LorgOptions get_options() const;

  //Grammar<Rule,Rule,Rule>* grammar;
  // ParserCKYBest* parser;
  // std::auto_ptr<Tagger> tagger;

  WordSignature * ws;


  bool train;
  treebank_options tb_options;

  unsigned cutoff;
  unsigned batch_size;
  unsigned iterations;
  unsigned lstm_level;
  unsigned span_level;

  std::string train_output_name;


  int run_train();

  std::pair<
    std::pair<std::vector<dynet::expr::Expression>,std::vector<dynet::expr::Expression>>,
    std::pair<std::string,std::string>>
  train_instance(const PtbPsTree& tree,
                 const ParserCKYNN& parser,
                 const Tagger& tagger,
                 ParserCKYNN::scorer& network,
                 int start_symbol,
                 const std::vector<int>& lhs_int_vec
                 );

  PtbPsTree *
  parse_instance(const std::vector<Word>& words,
                 const ParserCKYNN& parser,
                 ParserCKYNN::scorer& network,
                 int start_symbol,
                 const std::vector<int>& lhs_int_vec);



  //not defined -> forbidden
  NNLorgParseApp(const NNLorgParseApp&);
  NNLorgParseApp& operator=(const NNLorgParseApp);
};
