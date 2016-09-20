// -*- mode: c++ -*-
#pragma once

#include <vector>
#include <string>
#include <algorithm>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#pragma clang diagnostic pop

template <typename Bin, typename Un, typename Lex>
class Grammar
{
public:
  typedef Bin MyBinaryRule;
  typedef Un  MyUnaryRule;
  typedef Lex MyLexicalRule;

  std::vector<Bin> binary_rules;
  std::vector<Un>  unary_rules;
  std::vector<Lex> lexical_rules;

public:
  Grammar();
  virtual ~Grammar();

  /**
     \brief set rules
     \param binary_rules binary rules
     \param unary_rules unary rules
     \param lexical_rules lexical rules
  */
  void set_rules(const std::vector<Bin>&  binary_rules,
		 const std::vector<Un>&   unary_rules,
		 const std::vector<Lex>& lexical_rules);

  void init();

  Grammar(const std::string& filename);

  template <typename B, typename U, typename L>
  friend
  std::ostream& operator<<(std::ostream& out, const Grammar<B,U,L> & gram);

protected:

    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/)
    {
      ar & binary_rules;
      ar & unary_rules;
      ar & lexical_rules;
    }



};

#include <algorithm>


template<class Bin, class Un, class Lex>
Grammar<Bin, Un, Lex>::Grammar() :
  binary_rules(), unary_rules(), lexical_rules()
{}

template<class Bin, class Un, class Lex>
Grammar<Bin, Un, Lex>::~Grammar() {}


template<class Bin, class Un, class Lex>
void Grammar<Bin, Un, Lex>::init()
{
  // remove rules "overly low" probabilities
  binary_rules.erase(std::remove_if(binary_rules.begin(), binary_rules.end(),
                                    std::mem_fun_ref(&Bin::is_invalid)),
                     binary_rules.end());

  unary_rules.erase(std::remove_if(unary_rules.begin(), unary_rules.end(),
                                   std::mem_fun_ref(&Un::is_invalid)),
                    unary_rules.end());

  std::sort(binary_rules.begin(),binary_rules.end());
  std::sort(unary_rules.begin(),unary_rules.end());
  std::sort(lexical_rules.begin(),lexical_rules.end());
}


template<class Bin, class Un, class Lex>
inline
void Grammar<Bin, Un, Lex>::set_rules(const std::vector<Bin>& binary_rules_,
				      const std::vector<Un>&  unary_rules_,
				      const std::vector<Lex>& lexical_rules_)
{
  binary_rules = binary_rules_;
  unary_rules = unary_rules_;
  lexical_rules = lexical_rules_;;
}
