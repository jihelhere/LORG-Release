// -*- mode: c++ -*-
#pragma once

#include "../Word.h"
#include "SymbolTable.h"

#include <vector>
#include <boost/regex.hpp>

class Tagger
{
public:
  Tagger(const std::vector< std::vector<const MetaProduction*> >* word_2_rule = NULL);

  void tag( std::vector< Word >& sentence, const WordSignature& ws) const;
  void tag( Word& w, const WordSignature& ws ) const;


  void set_word_rules(const std::vector< std::vector<const MetaProduction*> >* word_2_rule);

private:
  const std::vector< std::vector<const MetaProduction*> >* word_2_rule_;
};
