// -*- mode: c++ -*-
#pragma once

#include <vector>
#include <iostream>
#include "utils/SymbolTable.h"
#include "lexicon/WordSignature.h"

class MetaProduction;

/**
  \class Word
  \brief represents a word in a sentence
  \author Wolfgang Seeker & others
*/
class Word
{
  friend class Tagger;
  friend class Tokeniser;

public:
  friend std::ostream& operator<<(std::ostream& out, const Word& word);

public:
// Word& operator=(const Word&);
// Word& operator=(Word&&);
// Word(Word&&);

  /**
     \brief Constructor
     \param str the given word form
     \param start_idx start posisition in the input sentence (may be a lattice)
     \param end_idx end position in the input lattice (-1 if sentence)
     \param pos the assigned pos tags
  */
  Word(const std::string& str, int start_idx, int end_idx = -1, const std::vector<std::string>& pos  = std::vector<std::string>());

  bool is_tagged() const;

  short get_given_tag(unsigned i) const;

  int get_id() const;

  const std::string& get_form() const {return form;};

  const std::vector<const MetaProduction*>& get_rules() const;

  int get_start() const;
  int get_end() const;

  void initialize_id(const WordSignature& wordsignature);

  void untag();


public:

protected:
  std::string form;                         ///< the actual word form

  int start;
  int end;

  int id;                                          ///< the assigned id
  int sigid;
  std::vector<short> tags;
  std::vector<const MetaProduction*> rules;   ///< the possible lexical rules with this word as rhs

private:
  Word();
};



inline
bool Word::is_tagged() const
{
  return !this->tags.empty();
}

inline
short Word::get_given_tag(unsigned i) const
{
  assert(i < this->tags.size());
  return this->tags[i];
}

inline
const std::vector<const MetaProduction*>& Word::get_rules() const
{
  return this->rules;
}

inline
int Word::get_id() const
{
  return this->id;
}

inline
int Word::get_start() const
{
  return this->start;
}

inline
int Word::get_end() const
{
  return this->end;
}


inline
void Word::untag()
{
  this->tags.clear();
}
