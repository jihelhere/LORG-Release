// -*- mode: c++ -*-
#pragma once
#include "LorgParseApp.h"

//#include "grammars/Grammar.h"
//#include "rules/Rule.h"
//#include "parsers/ParserCKYBest.h"

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

  int run_train();

  //not defined -> forbidden
  NNLorgParseApp(const NNLorgParseApp&);
  NNLorgParseApp& operator=(const NNLorgParseApp);
};
