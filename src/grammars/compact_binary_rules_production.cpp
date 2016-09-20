
/// template specialisations

#include "compact_binary_rules.cpp"
#include "rules/Production.h"

template<>
const Production* transform(const Production& rule)
{
  return &rule;
}


template
compact_binary_rules::vector_brules<const Production*>
compact_binary_rules::vector_brules<const Production*>::convert<Production>(std::vector<Production> const&);
