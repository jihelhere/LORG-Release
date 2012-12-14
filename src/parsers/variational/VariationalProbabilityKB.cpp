// -*- mode: c++ -*-
#include "VariationalProbabilityKB.h"

unsigned VariationalProbabilityKB::size = 0;

std::ostream & VariationalProbabilityKB::operator>> (std::ostream & out) const
{
  for(auto& cand: candidates) { out << "cand:" << cand.probability << " "; }
  return out;
}
