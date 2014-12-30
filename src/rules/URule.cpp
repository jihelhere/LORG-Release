#include "URule.h"

#include "utils/SymbolTable.h"
#include "utils/LorgConstants.h"

#include <cassert>
#include <numeric>
#include <iostream>

URule::URule(short l, short rhs0_, const std::vector<unary_proba_info>& probs) :
  AnnotatedRule(l), rhs0(rhs0_), probabilities()
{
  // the difficult thing here is to extend the vectors to the correct sze

  //  std::cout << probs.size() << std::endl;

  for(std::vector<unary_proba_info>::const_iterator it = probs.begin(); it != probs.end(); ++it) {

    //    std::cout << it->lhs_pos << " " << it->rhs0_pos << std::endl;
    if(it->lhs_pos >= (int) probabilities.size())
      probabilities.resize(it->lhs_pos+1);

    if(it->rhs0_pos >= (int) probabilities[it->lhs_pos].size()) {

      unsigned old_size = probabilities[it->lhs_pos].size();
      probabilities[it->lhs_pos].resize(it->rhs0_pos+1);

      for(unsigned i = old_size ; i < probabilities[it->lhs_pos].size(); ++i)
	probabilities[it->lhs_pos][i] = 0 ;

    }
    probabilities[it->lhs_pos][it->rhs0_pos] = it->probability;

  }
  //  std::cout << "urule created: " << *this << std::endl;

}

URule::URule(short l, short rhs0_, double prob) :
  AnnotatedRule(l), rhs0(rhs0_),  probabilities(std::vector< std::vector<double> >(1, std::vector<double>(1, prob))) {}

URule::URule(short l, short rhs0_, const std::vector< std::vector<double> >& probs)
  : AnnotatedRule(l), rhs0(rhs0_), probabilities(probs) {}

std::ostream& operator<<(std::ostream& out, const URule& rule)
{
  out.precision( 22 ) ;

  out << "int"
      << " " << SymbolTable::instance_nt().translate(rule.lhs)
      << " " << SymbolTable::instance_nt().translate(rule.rhs0);

  for(unsigned i = 0 ; i < rule.probabilities.size(); ++i)
    for(unsigned j = 0; j < rule.probabilities[i].size(); ++j)
      if(rule.probabilities[i][j] != 0)
	out << " (" << i << ',' << j << ',' << rule.probabilities[i][j] << ')';

  return out;
}


//#include <iostream>

#ifdef USE_MANUAL_SSE
#include <xmmintrin.h>  // Need this for SSE compiler intrinsics
#endif
void URule::update_inside_annotations(AnnotationInfo& up_annot,
                                      const AnnotationInfo& left_annot) const
{

#ifdef USE_MANUAL_SSE
  size_t size = left_annot.inside_probabilities.array.size();
  size_t N = size / 2;
  bool odd = (size % 2) != 0 and not left_annot.invalids[size -1];

  const auto& llast = left_annot.inside_probabilities.array[size -1];
#endif


  for(size_t i = 0 ; i < probabilities.size();++i)
  {
    if(up_annot.invalids[i]) continue;
    //if(probabilities[i].empty()) continue;

    const std::vector<double>& rule_probs_i = probabilities[i];
    if (rule_probs_i.empty()) continue;
    auto& ui = up_annot.unary_temp.array[i];


    //   std::cout << *this << std::endl;

#ifndef USE_MANUAL_SSE
    for(size_t j = 0 ; j < rule_probs_i.size();++j)
    {
      // this test seems useless ?
      //if(left_annot.invalids[j]) continue;
      ui += left_annot.inside_probabilities.array[j] * rule_probs_i[j];
    }
#else
    __m128d * a = (__m128d*) left_annot.inside_probabilities.array.data();
    __m128d * p = (__m128d*) rule_probs_i.data();
    double tmp_double[2] = {0.0,0.0};
    __m128d * tmp = (__m128d*) tmp_double;


    for (size_t j = 0; j < N; ++j, ++a, ++p)
    {
      //std::cerr << "k " << k << std::endl;
      _mm_store_pd(tmp_double,
                   _mm_add_pd( *tmp,
                               _mm_mul_pd(*a,*p)));
    }


    ui += tmp_double[0] + tmp_double[1];
    if (odd)
    {
      ui += llast* rule_probs_i[size - 1];
    }
#endif


     // if(up[i] < 0.0 || up[i] > 1.0)
    //std::cout << *this << " " << up[i] <<std::endl;

    assert(ui >= 0.0);
    assert(ui <= 1.0);
  }
}
// {
//   std::cout << *this << std::endl;
//   std::cout << "size up: " << up.size() << std::endl;
//   std::cout << "size left: " << left.size() << std::endl;
//   for(unsigned short i = 0 ; i < probabilities.size();++i) {
//     std::cout << i << std::endl;
//     const std::vector<double>& dim_i = probabilities[i];
//     up[i] += std::inner_product(dim_i.begin(), dim_i.end(), left.begin(), 0.0);
//   }
// }




//inline
void URule::update_outside_annotations(const AnnotationInfo& up_annot,
                                       AnnotationInfo& left_annot) const
// {
//   for(unsigned short i = 0 ; i < probabilities.size();++i) {
//     //if(up[i] == 0) continue;
//     //if(probabilities[i].empty()) continue;

//     const std::vector<double>& rule_probs_i = probabilities[i];
//     for(unsigned short j = 0 ; j < rule_probs_i.size();++j) {
//       //if(rule_probs_i[j] == 0) continue;
//       left[j] += up[i] * rule_probs_i[j];
//     }
//   }
// }
{
  for(auto i = 0U; i < probabilities.size();++i) {
    if(up_annot.invalids[i]) continue;
    const auto& dim_i = probabilities[i];
    for(auto j = 0U; j < dim_i.size();++j) {
      if(left_annot.invalids[j]) continue;
      left_annot.unary_temp.array[j] += up_annot.outside_probabilities.array[i] * dim_i[j];
    }
  }
}

double URule::update_outside_annotations_return_marginal(const AnnotationInfo& up_annot,
                                                         AnnotationInfo& left_annot) const
{
  double marginal = 0.0;
  for(unsigned short i = 0 ; i < probabilities.size();++i) {
    if(up_annot.invalids[i]) continue;
    const std::vector<double>& dim_i = probabilities[i];
    for(unsigned short j = 0 ; j < dim_i.size();++j) {
      if(left_annot.invalids[j]) continue;
      double delta = up_annot.outside_probabilities.array[i] * dim_i[j] ;
      left_annot.unary_temp.array[j] += delta ;
      marginal += delta * left_annot.inside_probabilities.array[j] ;
    }
  }
  return marginal;
}

void URule::compact()
{
  //get rid of lines of zeros
  for(unsigned i = 0; i < probabilities.size(); ++i) {
    bool allzeros = true;
    for(unsigned j = 0; j < probabilities[i].size(); ++j) {
      if(probabilities[i][j] != 0) {
  	allzeros = false;
  	break;
      }
    }
    if(allzeros)
    {
      probabilities[i].clear();
    }
    probabilities[i].shrink_to_fit();
    //std::vector<double>(probabilities[i].begin(), probabilities[i].end()).swap(probabilities[i]);
  }
}


void URule::uncompact(unsigned rhs_size)
{
  for(unsigned i = 0; i < probabilities.size(); ++i) {
    probabilities[i].resize(rhs_size,0);
  }
}

void URule::remove_unlikely_annotations(const double& threshold)
{
  bool changed = false;

  for(unsigned i = 0; i < probabilities.size(); ++i)
    for(unsigned j = 0; j < probabilities[i].size(); ++j)
      if(probabilities[i][j] < threshold) {
	probabilities[i][j] = 0.0;
	changed = true;
      }

  if(changed) compact();
}


bool URule::operator==(const URule& other) const
{
  return
    lhs  == other.lhs  &&
    rhs0 == other.rhs0 &&
    probabilities == other.probabilities;
}


bool URule::operator<(const URule& other) const
{
  return
    lhs < other.lhs ||
    (lhs  == other.lhs  &&
     ( rhs0 < other.rhs0 ||
     (rhs0 == other.rhs0 &&
      probabilities < other.probabilities)));
}

void URule::set_probability(const double& value)
{
  for(unsigned i = 0 ; i < probabilities.size(); ++i)
    for(unsigned j = 0; j < probabilities[i].size(); ++j)
	probabilities[i][j] = value;
}


// p_x = p(A_x -> B_y ), y fixed !!!
// p'_x = (1 - \alpha) p_x + \alpha (average_x' (p_x'))
// assume the rule is not compacted
void URule::linear_smooth(const double& alpha)
{
  if(probabilities.size() == 1) return;

  double alphabar = 1 - alpha;
  double other_factor = alpha / (probabilities.size() - 1);

  std::vector< std::vector<double> > new_probabilities(probabilities.size(),
						       std::vector<double>(probabilities[0].size(),0.0));

  for(unsigned i = 0; i < probabilities.size(); ++i)
    for(unsigned o = 0; o < probabilities.size(); ++o)
      for(unsigned j = 0; j < probabilities[i].size(); ++j) {
	new_probabilities[i][j] += ( i == o ? alphabar : other_factor) * probabilities[o][j];
      }

  probabilities = new_probabilities;
}


// assume the rule is not compacted
void URule::weighted_smooth(const double& alpha, const std::vector<std::vector<double> > & weights)
{
  if(probabilities.size() == 1) return;

  const std::vector<double> & myweights = weights[lhs];

  std::vector< std::vector<double> > new_probabilities(probabilities.size(),
						       std::vector<double>(probabilities[0].size(),0.0));

  for(unsigned i = 0; i < probabilities.size(); ++i) {

    double denominator = std::accumulate(myweights.begin(),myweights.end(), - myweights[i]);
    double alphabar = (1 - alpha) / myweights[i];
    double other_factor = alpha / denominator;

    for(unsigned o = 0; o < probabilities.size(); ++o)
      for(unsigned j = 0; j < probabilities[i].size(); ++j) {
	new_probabilities[i][j] += ( i == o ? alphabar : other_factor) * probabilities[o][j] * myweights[o];
      }
  }

  probabilities = new_probabilities;
}

// assume the rule is not compacted
void URule::generation_smooth(const std::vector<std::vector<std::vector<double> > > & weights)
{
  if(probabilities.size() == 1) return;

  const std::vector<std::vector<double> > & myweights = weights[lhs];

  std::vector< std::vector<double> > new_probabilities(probabilities.size(),
						       std::vector<double>(probabilities[0].size(),0.0));

  for(unsigned i = 0; i < probabilities.size(); ++i) {
    for(unsigned o = 0; o < probabilities.size(); ++o)
      for(unsigned j = 0; j < probabilities[i].size(); ++j) {
	new_probabilities[i][j] +=  probabilities[o][j] * myweights[o][i];
      }
  }

  probabilities = new_probabilities;
}
