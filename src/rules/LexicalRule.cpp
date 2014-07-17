#include "LexicalRule.h"

#include "utils/SymbolTable.h"


LexicalRule::LexicalRule() : rhs0(-1), probabilities() {}

LexicalRule::~LexicalRule() {}

LexicalRule::LexicalRule(short l, int rhs0_, const std::vector<lexical_proba_info>& probs) :
  AnnotatedRule(l), rhs0(rhs0_) , probabilities()
{
  for(std::vector<lexical_proba_info>::const_iterator it = probs.begin(); it != probs.end(); ++it) {

    if(it->lhs_pos >= (int) probabilities.size())
      probabilities.resize(it->lhs_pos+1);

    probabilities[it->lhs_pos] = it->probability;
  }
}

std::ostream& operator<<(std::ostream& out, const LexicalRule& rule)
{
  out.precision( 22 ) ;

  out << "lex"
      << " " << SymbolTable::instance_nt().translate(rule.lhs)
      << " " << SymbolTable::instance_word().translate(rule.rhs0);

  for(unsigned i = 0 ; i < rule.probabilities.size(); ++i)
    if(rule.probabilities[i] != 0)
      out << " (" << i << ',' << rule.probabilities[i] << ')';
  return out;
}

#ifdef USE_MANUAL_SSE
#include <xmmintrin.h>  // Need this for SSE compiler intrinsics
#endif

void LexicalRule::update_inside_annotations(AnnotationInfo& annot) const
{
  //std::cout << "probabilities size " << probabilities.size() << std::endl;

#ifndef USE_MANUAL_SSE
  size_t size = probabilities.size();
  for(size_t i = 0 ; i < size;++i)
  {
    if(annot.invalids[i]) continue;
    annot.inside_probabilities.array[i] += probabilities[i];
  }
#else
  size_t size = probabilities.size();
  size_t N = size / 2;
  bool odd = (size % 2) != 0;

  __m128d * a = (__m128d*) annot.inside_probabilities.array.data();
  __m128d * p = (__m128d*) probabilities.data();

  for (size_t i = 0; i < N; ++i, ++a, ++p)
  {
    _mm_store_pd((double*)a, _mm_add_pd(*a,*p));
  }

  if (odd and not annot.invalids[size-1])
  {
    annot.inside_probabilities.array[size-1] += probabilities[size-1];
  }
#endif

}


#include <numeric>

//inline
void LexicalRule::update_outside_annotations(const AnnotationInfo& annot,
                                             double& left) const
{
  for(unsigned i = 0 ; i < probabilities.size();++i) {
    if(annot.invalids[i]) continue;
    //if( up[i] == 0 ||      probabilities[i] == 0) continue;
    left += annot.inside_probabilities.array[i] * probabilities[i];
  }

  //  left += std::inner_product(probabilities.begin(), probabilities.end(), up.begin(), 0.0);
}

double LexicalRule::update_outside_annotations_return_marginal(const std::vector< double >& up,
                                                               const std::vector<bool>& invalids) const
{
  double marginal = 0;
  for(unsigned i = 0 ; i < probabilities.size();++i) {
    if(invalids[i]) continue;
    //if( up[i] == 0 ||      probabilities[i] == 0) continue;
    double delta = up[i] * probabilities[i];
    marginal += delta;
  }
  return marginal;
}


void LexicalRule::remove_unlikely_annotations(const double& threshold)
{
  for(unsigned i = 0; i < probabilities.size(); ++i) {
    if(probabilities[i] < threshold)
      probabilities[i] = 0;
  }
}



bool LexicalRule::operator==(const LexicalRule& other) const
{
  return
    lhs  == other.lhs  &&
    rhs0 == other.rhs0
    //&& probabilities == other.probabilities
    ;
}

bool LexicalRule::operator<(const LexicalRule& other) const
{
  return
    lhs  < other.lhs  ||
    (lhs == other.lhs &&
     rhs0 < other.rhs0);
}

void LexicalRule::set_probability(const double& value)
{
  for(unsigned i = 0 ; i < probabilities.size(); ++i)
    probabilities[i] = value;
}

void LexicalRule::linear_smooth(const double& alpha)
{
  if(probabilities.size() == 1) return;

  double alphabar = 1 - alpha;
  double other_factor = alpha / (probabilities.size() - 1);

  std::vector<double> new_probabilities(probabilities.size(),0.0);

  for(unsigned i = 0; i < probabilities.size(); ++i) {

    for(unsigned o = 0; o < probabilities.size(); ++o) {
      new_probabilities[i] += ( i == o ? alphabar : other_factor) * probabilities[o] ;
    }
  }

  probabilities = new_probabilities;
}



void LexicalRule::weighted_smooth(const double& alpha, const std::vector<std::vector<double> > & weights)
{
  if(probabilities.size() == 1) return;

  // std::cout << lhs << " " << weights.size() << std::endl;
  // std::cout << weights[lhs].size() << std::endl;

  const std::vector<double> & myweights = weights[lhs];

  std::vector<double> new_probabilities(probabilities.size(),0.0);

  for(unsigned i = 0; i < probabilities.size(); ++i) {

    double denominator = std::accumulate(myweights.begin(),myweights.end(), - myweights[i]);
    double alphabar = (1 - alpha) / myweights[i];
    double other_factor = alpha / denominator;



    for(unsigned o = 0; o < probabilities.size(); ++o) {
      //      std::cout << "weights[lhs][o] " << weights[lhs][o] << std::endl;
      new_probabilities[i] += ( i == o ? alphabar : other_factor) * probabilities[o] * myweights[o] ;
    }
    //    std::cout << "new_probabilities[i] " << new_probabilities[i] << std::endl;
  }

  probabilities = new_probabilities;
}



void LexicalRule::generation_smooth(const std::vector<std::vector<std::vector<double> > >& weights)
{
  if(probabilities.size() == 1) return;

  const std::vector<std::vector<double> > & myweights = weights[lhs];

  std::vector<double> new_probabilities(probabilities.size(),0.0);

  // std::cout << "before" << std::endl;
  // for(unsigned i = 0; i < probabilities.size(); ++i)
  //   std::cout << probabilities[i] << " ";
  // std::cout << std::endl;

  for(unsigned i = 0; i < probabilities.size(); ++i) {
    for(unsigned o = 0; o < probabilities.size(); ++o) {
      //      std::cout << myweights[i][o] << " " << probabilities[o] << std::endl;

      new_probabilities[i] += myweights[o][i] * probabilities[o];
    }
  }
  probabilities = new_probabilities;

  // std::cout << "after" << std::endl;
  // for(unsigned i = 0; i < probabilities.size(); ++i)
  //   std::cout << probabilities[i] << " ";
  // std::cout << std::endl;

}
