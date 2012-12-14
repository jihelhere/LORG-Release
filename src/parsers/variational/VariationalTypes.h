#ifndef _VARIATIONALTYPES_H_
#define _VARIATIONALTYPES_H_

#include "rules/BRuleC2f.h"
#include "rules/LexicalRuleC2f.h"
#include "rules/URuleC2f.h"



class VariationalProbabilityKB;
class VariationalEdgeDaughterProbability;

struct VariationalTypes {
  typedef VariationalProbabilityKB EdgeProbability ;
  typedef VariationalEdgeDaughterProbability EdgeDaughterProbability ;
  typedef Word ChartWord ;

  typedef PackedEdge< VariationalTypes > Edge ;
  typedef PCKYAllCell< VariationalTypes > Cell ;
  typedef ChartCKY< VariationalTypes > Chart ;

  typedef UnaryPackedEdgeDaughters<VariationalTypes> UnaryDaughter;
  typedef BinaryPackedEdgeDaughters<VariationalTypes> BinaryDaughter;
  typedef LexicalPackedEdgeDaughters<VariationalTypes> LexicalDaughter;

  typedef BRuleC2f BRule;
  typedef URuleC2f URule;
  typedef LexicalRuleC2f LRule;
};



#endif /* _VARIATIONALTYPES_H_ */
