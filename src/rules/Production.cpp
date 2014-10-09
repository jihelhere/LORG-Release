#include "Production.h"
#include "utils/SymbolTable.h"

#ifdef DEBUG
#include "boost/lexical_cast.hpp"
#else
#include <sstream>
#endif


Production::Production(const int llhs,const std::vector<int>& rrhs,const bool llexical)
  : MetaProduction(llhs), rhs(rrhs), lexical(llexical) {}

Production::~Production()
{
}

// format is: {lex|int} LHS RHS_1 ... RHS_n
std::ostream& operator<<(std::ostream& os, const Production& prod)
{
  if (prod.lexical)
    os << "lex ";
  else
    os << "int ";

  os << SymbolTable::instance_nt().translate(prod.lhs);

  for(const auto& p : prod.rhs)
  {
    if (prod.lexical)
      os <<" " << SymbolTable::instance_word().translate(p);
    else
      os <<" " << SymbolTable::instance_nt().translate(p);
  };

  return os;
}

#ifdef DEBUG
std::string create_debug_name(std::vector<int> rhs,bool right)
{
  std::string new_name;

  if(right) {
    for(std::vector<int>::iterator i = rhs.begin(); i != rhs.end() - 1 ; ++i)
      new_name += "(" + SymbolTable::instance_nt()->translate(*i) + ")";
  }
  else{ //left
    for(std::vector<int>::iterator i = rhs.begin()+1; i != rhs.end() ; ++i)
      new_name += "(" + SymbolTable::instance_nt()->translate(*i) + ")";
  }

  new_name = "["+new_name+"]";

  return new_name;
}
#endif

std::string create_name(std::vector<int> rhs, bool right)
{
  // to name newly created non-terminals
  static int cpt = 0; // static useful ? (yes if we call this method twice...)
  static std::map<std::vector<int>,int> record;

  std::ostringstream new_name;
  std::map<std::vector<int>, int>::iterator ans;
  std::vector<int> keyvector;

  if(right)
    keyvector.insert(keyvector.end(),rhs.begin(),rhs.end()-1);
  else // left
    keyvector.insert(keyvector.end(),rhs.begin()+1,rhs.end());

  if ((ans=record.find(keyvector)) == record.end()) {
    record[keyvector] = ++cpt;
    new_name << "[" << cpt << "]";
  }
  else
    new_name << "[" << ans->second  << "]";

  return new_name.str();

}
