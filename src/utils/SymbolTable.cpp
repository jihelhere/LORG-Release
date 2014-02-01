#include "SymbolTable.h"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#pragma clang diagnostic pop

#include <fstream>


SymbolTable * SymbolTable::NT_instancePtr = nullptr;//initialize singlton pointer
SymbolTable * SymbolTable::word_instancePtr = nullptr;//initialize singlton pointer
std::string const SymbolTable::unknown_string  = "UNK"; // TODO: give a more generic unknown string

SymbolTable& SymbolTable::instance_word()
{
  if (word_instancePtr == nullptr)  // is it the first call?
    word_instancePtr = new SymbolTable(); // create sole instance
  return *word_instancePtr; // address of sole instance
}

SymbolTable& SymbolTable::instance_nt()
{
  if (NT_instancePtr == nullptr)  // is it the first call?
    NT_instancePtr = new SymbolTable(); // create sole instance
  return *NT_instancePtr; // address of sole instance
}

SymbolTable::SymbolTable() : table(), cpt(0) {}

unsigned int SymbolTable::get(const std::string& str) const throw(Miss)
{
  return get_label_id(str);
}

unsigned int SymbolTable::get_label_id(const std::string& str) const throw(Miss)
{
  symtab::left_const_iterator find_str = table.left.find(str);
  if (find_str==table.left.end())
    throw(Miss(str));
  else
    return find_str->second;
}

int SymbolTable::insert(const std::string& str) throw()
{
  try {
    return get_label_id(str);
  }
  catch (Miss) {
    table.insert(symtab::value_type(str,cpt));
    return cpt++;
  }
}

std::string SymbolTable::translate(unsigned int i) const throw(Miss)
{
  return get_label_string(i);
}

std::string SymbolTable::get_label_string(unsigned int i) const throw(Miss)
{
  symtab::right_const_iterator find_i(table.right.find(i));
  if (find_i==table.right.end())
    throw(Miss(i));
  else
    return find_i->second;
}

void SymbolTable::load(std::string filename)
{
  std::ifstream ifs(filename);
  boost::archive::text_iarchive ia(ifs);
  ia >> *this;
}

void SymbolTable::load(SymbolTable& st)
{
  *this = st;
}



std::unordered_map<int,int> SymbolTable::build_simplification_map()
{
  std::unordered_map<int,int> result;

  typedef symtab::left_map map_type;
  map_type & m = table.left;

  static const boost::regex exp_artificial ("^\\[\\((.*)\\)>\\]$");
  static const boost::regex exp_funct ("^([A-Za-z]+\\$?)[=-].*");


  for (const auto& pair : m)
  {
    std::string name = pair.first;
    bool is_artificial = false;
    boost::cmatch matched;

    if(boost::regex_match(name.c_str(), matched, exp_artificial))
    {
      name = std::string(matched[1].first, matched[1].second);
      is_artificial = true;
    }

    if(boost::regex_match(name.c_str(),matched,exp_funct))
    {
      name = std::string(matched[1].first, matched[1].second);
    }

    if(is_artificial)
    {
      name = "[(" + name + ")>]";
    }

    result[pair.second] = this->get_label_id(name);
  }

  return result;
}

std::vector<unsigned> SymbolTable::get_mwe_symbols() const
{
  std::vector<unsigned> results;

  typedef symtab::left_map map_type;
  const map_type & m = table.left;

  for (const auto& pair : m)
  {
    //if(pair.first[0] == 'M' and pair.first[1] == 'W')
    if(pair.first[pair.first.size()-1] == '+' )
      results.push_back(pair.second);
  }

  //std::cerr << "mwe symbols.size: " << results.size() << std::endl;

  return results;
}
