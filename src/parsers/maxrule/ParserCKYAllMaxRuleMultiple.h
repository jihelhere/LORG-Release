// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVARMULTIPLE_H_
#define _PARSERCKYALLMAXVARMULTIPLE_H_

#include "ParserCKYAllMaxRule.h"
#include "MaxRuleProbabilityMultiple.h"


class ParserCKYAllMaxRuleMultiple : public ParserCKYAllMaxRule<MaxRuleMultipleTypes>
{
public:
  ParserCKYAllMaxRuleMultiple(ParserCKYAllFactory::MaxParsing_Calculation c,
                              std::vector<AGrammar*>& cgs,
                              const std::vector<double>& p, double b_t,
                              const std::vector< std::vector<AGrammar*> >& fgs,
                              const std::vector<annot_descendants_type>& all_annot_descendants_,
                              bool accurate_, unsigned min_beam, int stubborn, unsigned k);

  ~ParserCKYAllMaxRuleMultiple();

  /**
     \brief wraps the calculation of the best derivation
  */
  void extract_solution();

  void simple_extract_solution();

  const AGrammar& get_fine_grammar(unsigned i, unsigned j) const;

protected:
  /**
     \brief compute scores with all the fine grammars and back them up in the chart
  */
  void precompute_all_backups();

  void multiple_inside_outside_specific();


  /**
     \brief replace rules with their followers according to the defined mapping
     and reset annotations to zero (and resize thme to 1)
   */
  void change_rules_reset();

  /**
     \brief replace rules with their followers + size_grammar (to skip intermediate grammars)
     and replace current annotations with backed up ones at position backup_idx
  */
  void change_rules_load_backup(unsigned backup_idx, unsigned size_grammar) const;


  void modify_backup(unsigned backup_idx) const;

  /**
     \brief Calculates the chart specific rule probabilities of the packed edges in the chart
  */
  void calculate_maxrule_probabilities();

  /**
     \brief pick up the best derivation once the edge scores have been calculated
   */
  void calculate_best_edge();

  /**
     \brief for all edges in chart, backup current annotation
   */
  void backup_annotations() const;


protected: // attributes
  std::vector< std::vector<AGrammar*> >fine_grammars; ///< the additional grammars to be used to extract the solution
  std::vector<annot_descendants_type> all_annot_descendants; ///< all the annotations mapping for the grammars (base + fine ones)


private:
  unsigned nb_grammars;
  unsigned k;

  double log_normalisation_factor;
  std::vector<double> log_normalisation_factor_backup;

  void initialise_candidates();
  void extend_all_derivations();

  void set_log_normalisation_factor(double lnf);
  void reset_log_normalisation_factor();
  const double& get_log_normalisation_factor();
  inline const double& get_log_normalisation_factor(unsigned i);
};

inline
const ParserCKYAll::AGrammar& ParserCKYAllMaxRuleMultiple::get_fine_grammar(unsigned i, unsigned j) const
{
  return *fine_grammars[i][j];
}

const double& ParserCKYAllMaxRuleMultiple::get_log_normalisation_factor(unsigned i)
{
  //  std::cout << "size: " << log_normalisation_factor_backup.size() << std::endl;;
  return log_normalisation_factor_backup[i];
}

#endif /* _PARSERCKYALLMAXVARMULTIPLE_H_ */
