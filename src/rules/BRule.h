// -*- mode: c++ -*-
#ifndef BRULE_H
#define BRULE_H

#include <vector>
#include "AnnotatedRule.h"
#include "Types4Rules.h"

#include "utils/LorgConstants.h"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "edges/AnnotationInfo.h"

typedef std::vector< std::vector< std::vector<double> > > vector_3d;

/**
   \class BRule
   \brief a PCFG-LA binary rule
 */
class BRule : public AnnotatedRule
{

public:
  BRule() : rhs0(-1), rhs1(-1), probabilities() {};
  virtual ~BRule() {};

  /**
     \brief constructor for creating binary rules from the grammar file
     \param l the lhs of the object
     \param rhs0  the leftmost symbol of rhs
     \param rhs1  the rightmost symbol of rhs
     \param probs a vector of binary_proba_info 4-tuple (i,j,k,p)
  */
  BRule(short l, short rhs0, short rhs1, const std::vector<binary_proba_info>& probs);

  /**
     \brief constructor for the creation of a base grammar binary rule
     \param l the lhs of the object
     \param rhs0  the leftmost symbol of rhs
     \param rhs1  the rightmost symbol of rhs
     \param double probability
   */
  BRule(short l, short rhs0, short rhs1, double prob);

  /**
     \brief constructor
     \param l the lhs of the object
     \param rhs0  the leftmost symbol of rhs
     \param rhs1  the rightmost symbol of rhs
     \param probs a 3d-vector of probabilities
  */
  BRule(short l, short rhs0, short rhs1, const vector_3d& probs);



  /**
     \brief returns attribute rhs0
  */
  short get_rhs0() const;

  /**
     \brief returns attribute rhs1
  */
  short get_rhs1() const;

void set_rhs0(short r) {rhs0 = r;}
void set_rhs1(short r) {rhs1 = r;}



  /**
     \brief read access to annotated probabilities
     \param a annotation for lefthandside symbol
     \param b annotation for leftmost righthandside symbol
     \param c annotation for rightmost righthandside symbol
   */
  const double& get_probability(unsigned short a, unsigned short b, unsigned short c) const;

  /**
     \brief write access to annotated probabilities
     \param a annotation for lefthandside symbol
     \param b annotation for leftmost righthandside symbol
     \param c annotation for rightmost righthandside symbol
     \param value the new probability value
  */
  void set_probability(unsigned short a, unsigned short b, unsigned short c, const double& value);
  void set_probability(const double& value);

  /**
     \brief read/write access to the 3d-vector of probabilities
     \todo should be private ?
  */
  const vector_3d& get_probability() const;
  vector_3d& get_probability();

  /**
     \brief -> always false
   */
  inline bool is_lexical() const {return false;}
  /**
     \brief -> always false
   */
  inline bool is_unary() const {return false;}
  /**
     \brief -> always true
   */
  inline bool is_binary() const {return true;}


  /**
     \brief smooth probabilities *over all* lhs annotations
     \param alpha the small constant
  */
  void linear_smooth(const double& alpha);
  void weighted_smooth(const double& alpha, const std::vector<std::vector<double> > & weights);
  void generation_smooth(const std::vector<std::vector<std::vector<double> > > & weights);



  /**
     \brief Output operator
     \param out the ostream to write on
     \param rule the urule object to write
     \return the used ostream
  */
  friend std::ostream& operator<<(std::ostream& out, const BRule& rule);


  void update_inside_annotations(AnnotationInfo& up_annot,
                                 const AnnotationInfo& left_annot,
                                 const AnnotationInfo& right_annot) const;

  void update_inside_annotations(std::vector<double>& up,
                                 const double& left_right_precomputation) const;


  void update_outside_annotations(const AnnotationInfo& up_annot,
                                  AnnotationInfo& left_annot,
                                  AnnotationInfo& right_annot
                                  ) const;
  double update_outside_annotations_return_marginal(const AnnotationInfo& up_annot,
                                                    AnnotationInfo& left_annot,
                                                    AnnotationInfo& right_annot
                                                    ) const;
  /**
     \brief removes useless zeros from probability vector
  */
  void compact();

  /**
     \brief resize the probability vector by adding zeros
     \param rhs0_size new size for the the first right-hand side symbol
     \param rhs1_size new size for the the second right-hand side symbol
  */
  void uncompact(unsigned rhs0_size, unsigned rhs1_size);

  /**
     \brief replace annotations with 0 if they're below threshold
     \param threshold the threshold
  */
  void remove_unlikely_annotations(double threshold);


  bool is_empty() const {return probabilities.empty();}

  bool operator==(const BRule&) const;
  bool operator<(const BRule& other) const;


  bool is_invalid() const;

protected:
  short rhs0;
  short rhs1;
  vector_3d probabilities ; ///< probabilities for a CFG rule  with annotations
private:

    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/)
    {
      ar & boost::serialization::base_object<AnnotatedRule>(*this);
      ar & rhs0;
      ar & rhs1;
      ar & probabilities;
    }



};

inline
short BRule::get_rhs0() const
{
  return rhs0;
}

inline
short BRule::get_rhs1() const
{
  return rhs1;
}


inline
const double& BRule::get_probability(unsigned short a, unsigned short b, unsigned short c) const
{
  return probabilities[a][b][c];
}

inline
void BRule::set_probability(unsigned short a, unsigned short b, unsigned short c, const double& value)
{
  probabilities[a][b][c]=value;
}

inline
const vector_3d& BRule::get_probability() const
{
  return probabilities;
}

inline
vector_3d& BRule::get_probability()
{
  return probabilities;
}

inline
bool BRule::is_invalid() const
{
  return
    probabilities.empty() ||
    (probabilities.size() == 1 && probabilities[0].empty()) ||
    (probabilities.size() == 1 && probabilities[0].size() == 1 && probabilities[0][0].empty());
}



#endif //BRULE_H
