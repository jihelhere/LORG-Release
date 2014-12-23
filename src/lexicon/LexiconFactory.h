// -*- mode: c++ -*-

#ifndef LEXICONFACTORY_H
#define LEXICONFACTORY_H

#include "utils/ConfigTable.h"

#include "Lexicon.h"
#include "BasicLexicon.h"
#include "BerkeleyLexicon.h"

#include "WordSignatureFactory.h"


namespace LexiconFactory {

  enum lex_type {Basic, Bsophisticated};


  lex_type string_2_lex_type(const std::string& s);

  std::string lex_type_2_string(const lex_type& l);

  Lexicon * create_lexicon(lex_type type, const std::shared_ptr<WordSignature>& ws,
                           unsigned cutoff, bool verbose);

  Lexicon * create_lexicon(ConfigTable& config);
}









#endif // LEXICONFACTORY_H
