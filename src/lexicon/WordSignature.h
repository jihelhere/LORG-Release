// -*- mode: c++ -*-
#ifndef WORDSIGNATURE_H
#define WORDSIGNATURE_H

#include "utils/SymbolTable.h"
#include <string>
#include <sstream>


  enum unknown_word_map {Generic,
			 BerkeleyEnglish, BaselineFrench, Arabic,
			 EnglishIG, FrenchIG, ArabicIG, ItalianIG};

  unknown_word_map string_2_lex_unknown_map(const std::string& s);
  std::string lex_unknown_map_2_string(const unknown_word_map& m);


class WordSignature{
 private:
  unknown_word_map type;

protected:
  static bool is_upper_case_letter(char c){ return ((c >= 'A') && (c <= 'Z'));}

  static bool is_lower_case_letter(char c){ return ((c >= 'a') && (c <= 'z'));}

  static std::string to_lower_case(const std::string& word)
  {
    std::string new_word(word);
    for(std::string::iterator c = new_word.begin(); c != new_word.end(); ++c){
      *c = tolower(*c);
    }
    //std::cout << new_word << std::endl;
    return new_word;
  }

public:
  WordSignature(unknown_word_map t) : type(t) {};

  // static unknown_word_map string_2_lex_unknown_map(const std::string& s);
  // static std::string lex_unknown_map_2_string(const unknown_word_map& m);

static unknown_word_map
string_2_lex_unknown_map(const std::string& s)
{
  if(s == "BerkeleyEnglish")
    return BerkeleyEnglish;
  if(s == "BaselineFrench")
    return BaselineFrench;
  if(s == "EnglishIG")
    return EnglishIG;
  if(s == "FrenchIG")
    return FrenchIG;
  if(s == "Arabic")
    return Arabic;
  if(s == "ArabicIG")
    return ArabicIG;

  if(s == "ItalianIG")
    return ItalianIG;

  //if(s == "generic")
  return Generic;
}

static std::string
lex_unknown_map_2_string(const unknown_word_map& m)
{
  if(m == BerkeleyEnglish)
    return "BerkeleyEnglish";
  if(m == BaselineFrench)
    return "BaselineFrench";
  if(m == EnglishIG)
    return "EnglishIG";
  if(m == FrenchIG)
    return "FrenchIG";
  if(m == Arabic)
    return "Arabic";
  if(m == ArabicIG)
    return "ArabicIG";

  if(m == ItalianIG)
    return "ItalianIG";

  return "generic";
}



  unknown_word_map get_type() const {return type;}

  virtual std::string get_unknown_mapping(const std::string& word, unsigned position) const =0;
  virtual ~WordSignature() {};
};


class TrivialWordSignature : public WordSignature
{
public:
  TrivialWordSignature() :
      WordSignature(Generic)
  {};


  std::string get_unknown_mapping(const std::string& /*word*/, unsigned /*position*/) const {
	//  std::cout << "where is this being called? " << std::endl;
	 // std::exit(1);
	  return SymbolTable::unknown_string;
  };
};






#endif // WORDSIGNATURE_H
