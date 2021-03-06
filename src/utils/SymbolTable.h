// -*- mode: c++ -*-
#ifndef SYMBOLTABLE_H
#define SYMBOLTABLE_H

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wredeclared-class-member"
#include <boost/bimap/bimap.hpp>
#pragma clang diagnostic pop

#include <string>
#include <exception>
#include <map>
#include "utils/LorgConstants.h"

#include <sstream>

#include <unordered_map>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-register"
#include <boost/regex.hpp>
#pragma clang diagnostic pop

#include <boost/serialization/map.hpp>

/**
   \class Miss
   \brief used as exception in class \ref SymbolTable
 */
class Miss : public std::exception
{
  std::string msg;

public:
  virtual const char* what() const throw()
  {
    return msg.c_str();
  };

public:

  Miss(const std::string& s)
    : msg()
  {
    msg = "can't find string \"" + s + "\" in symbol table";
  };

  Miss(const int i)
  : msg()
  {
    std::ostringstream op;
    op << i;
    msg = "can't find int \"" + op.str() + "\" in symbol table";
  };

  ~Miss()throw() {};

private:
  Miss();
  Miss& operator=(const Miss& other);
};


/**
    \class SymbolTable
    \brief makes the  string <--> int translation easier

    This class provides 3 methods
    to ease translations between strings and their encoding as integers.
*/
class SymbolTable
{
private:
  typedef boost::bimaps::bimap<std::string, unsigned int> symtab;

  symtab table;  ///< string to integer bijective mapping
  unsigned int cpt;       ///< numbers of strings inserted so far


  static std::shared_ptr<SymbolTable> NT_instancePtr;   ///< pointer to the one instance of SymbolTable for non terminals
  static std::shared_ptr<SymbolTable> word_instancePtr; ///< pointer to the one instance of SymbolTable for terminals


private:
  /**
    \brief private constructor
  */

  SymbolTable(SymbolTable const& ); //not defined, not copyable
  SymbolTable& operator=(SymbolTable& st) {if(this != &st) {table =st.table; cpt =st.cpt;} return *this;};



  friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/)
    {
      ar & table;
      ar & cpt;
    }
public:
  //should not be public
  SymbolTable();

  /**
       \brief Destructor
    */
  ~SymbolTable() {};

  int insert(const std::string& str) throw();
  /**
     \brief Depreciated! See get_string_label() translates the given integer into the associated symbol or returns
     the empty string
    \param i integer to translate
    \return the symbol associated with \param i or the empty string
    \todo raise error instead of returning empty strings
  */

  std::string translate(unsigned int i) const throw(Miss);

  /**
     \brief translates the given integer into the associated symbol or returns
     the empty string
    \param i integer to translate
    \return the symbol associated with \param i or the empty string
    \todo raise error instead of returning empty strings
  */

  std::string get_label_string(unsigned int i) const throw(Miss);


  /**
     \brief returns the integer for a given symbol
     \param str symbol
     \return the assigned integer to the given symbol
   */

  unsigned int get_label_id(const std::string& str) const throw(Miss);
  /**
     \brief Depreciated - see get_label_id returns the integer for a given symbol
     \param str symbol
     \return the assigned integer to the given symbol
   */

  unsigned int get(const std::string& str) const throw(Miss);

  bool token_exists(const std::string&);

  /**
    \brief Global point of access to SymbolTable for terminals.

    If this is the first call to instance() it will create a SymbolTable object and return a pointer
    to it.  Otherwise it will just return the pointer to the object.
  */
  static SymbolTable& instance_word();
  /**
    \brief Global point of access to SymbolTable for nonterminals.

    If this is the first call to instance() it will create a SymbolTable object and return a pointer
    to it.  Otherwise it will just return the pointer to the object.
  */
  static SymbolTable& instance_nt();

  //move constructor
  SymbolTable(SymbolTable&& o) : table(std::move(o.table)), cpt(std::move(o.cpt)) {};


  /**
     \brief Returns the number of distinct symbols in the table
  */
  unsigned int get_symbol_count() const;

  bool is_root_label(int id);

  static const std::string unknown_string; ///< the string corresponding to the unknown token

  size_t get_size() const {return table.size();}


  void load(std::string filename);
  void load(SymbolTable& st);

  std::unordered_map<int,int> build_simplification_map();

  std::vector<unsigned> get_mwe_symbols() const;
};

inline
unsigned int SymbolTable::get_symbol_count() const {return table.size();}

//TODO what if this is called as word instance?
inline
bool SymbolTable::is_root_label(int id) {return id == (int) get_label_id(LorgConstants::tree_root_name);}

inline
bool SymbolTable::token_exists(const std::string& token){return table.left.find(token) != table.left.end();}

#endif // SYMTAB_H
