#include "BRule.h"
#include "utils/SymbolTable.h"
#include "utils/LorgConstants.h"

#include <cassert>
#include <numeric>
//#include <iostream>



BRule::BRule(short l, short rhs0_, short rhs1_, const std::vector<binary_proba_info>& probs) :
  AnnotatedRule(l), rhs0(rhs0_), rhs1(rhs1_), probabilities()
{
  // the difficult thing is to extend the vectors if needed

  for(const auto& bprob : probs)
  {

    if(bprob.lhs_pos >= (int) probabilities.size())
      probabilities.resize(bprob.lhs_pos+1);

    if(bprob.rhs0_pos >= (int) probabilities[bprob.lhs_pos].size())
      probabilities[bprob.lhs_pos].resize(bprob.rhs0_pos+1);

    if(bprob.rhs1_pos >= (int) probabilities[bprob.lhs_pos][bprob.rhs0_pos].size())
    {

      unsigned old_size = probabilities[bprob.lhs_pos][bprob.rhs0_pos].size();
      probabilities[bprob.lhs_pos][bprob.rhs0_pos].resize(bprob.rhs1_pos+1);

      for(unsigned i = old_size ; i < probabilities[bprob.lhs_pos][bprob.rhs0_pos].size(); ++i)
	probabilities[bprob.lhs_pos][bprob.rhs0_pos][bprob.rhs1_pos] = 0.0 ;

      probabilities[bprob.lhs_pos][bprob.rhs0_pos][bprob.rhs1_pos] = bprob.probability;
    }
  }
}

BRule::BRule(short l, short rhs0_, short rhs1_, double prob) :
  AnnotatedRule(l), rhs0(rhs0_), rhs1(rhs1_),
  probabilities(vector_3d(1, std::vector< std::vector<double> >(1,std::vector<double>(1,prob))))
{}

BRule::BRule(short l, short rhs0_, short rhs1_, const std::vector<std::vector<std::vector<double> > >& probs)
  : AnnotatedRule(l), rhs0(rhs0_),rhs1(rhs1_), probabilities(probs) {}


std::ostream& operator<<(std::ostream& out, const BRule& rule)
{
  out.precision( 22 ) ;

  out << "int"
      << " " << SymbolTable::instance_nt().translate(rule.lhs)
      << " " << SymbolTable::instance_nt().translate(rule.rhs0)
      << " " << SymbolTable::instance_nt().translate(rule.rhs1);


  for(unsigned i = 0 ; i < rule.probabilities.size(); ++i)
    for(unsigned j = 0; j < rule.probabilities[i].size(); ++j)
      for(unsigned k = 0; k < rule.probabilities[i][j].size(); ++k)
        if(rule.probabilities[i][j][k] != 0)
	  out << " (" <<i << ',' << j << ',' << k << ','
	      << rule.probabilities[i][j][k] << ')';

  return out;
}

#ifdef USE_MANUAL_SSE
#include <xmmintrin.h>  // Need this for SSE compiler intrinsics
#endif

void BRule::update_inside_annotations(AnnotationInfo& up_annot,
                                      const AnnotationInfo& left_annot,
                                      const AnnotationInfo& right_annot
                                      ) const
{


#ifdef USE_MANUAL_SSE
  size_t size = right_annot.inside_probabilities.array.size();
  size_t N = size / 2;
  bool odd = (size % 2) != 0  and not right_annot.invalids[size-1];
  auto& rlast = right_annot.inside_probabilities.array[size -1];
#endif //USE_MANUAL_SSE

  for(size_t i = 0 ; i < probabilities.size();++i) {
    if(up_annot.invalids[i]) continue;
    auto& ui = up_annot.inside_probabilities.array[i];
    const auto& probabilities_i = probabilities[i];

    for(size_t j = 0 ; j < probabilities_i.size(); ++j)
    {
      if(left_annot.invalids[j] || left_annot.inside_probabilities.array[j] == 0.0) continue;
      const auto& probabilities_ij = probabilities_i[j];

#ifndef USE_MANUAL_SSE
      double inner = 0.0;
      for(size_t k = 0 ; k < probabilities_ij.size();++k) {
        if(right_annot.invalids[k]) continue;
        inner += right_annot.inside_probabilities.array[k] * probabilities_ij[k];
      }
#else


      // std::cerr << *this << std::endl;
      // std::cerr << "size " << size << std::endl;

      if (probabilities[i][j].empty()) continue;

      __m128d * a = (__m128d*) right_annot.inside_probabilities.array.data();
      __m128d * p = (__m128d*) probabilities[i][j].data();
      double tmp_double[2] = {0.0,0.0};
      __m128d * tmp = (__m128d*) tmp_double;

      for (size_t k = 0; k < N; ++k, ++a, ++p)
      {
        //std::cerr << "k " << k << std::endl;
        _mm_store_pd(tmp_double,
                     _mm_add_pd( *tmp,
                                 _mm_mul_pd(*a,*p)));
      }

      double inner = tmp_double[0] + tmp_double[1];
      if (odd)
      {
        inner += rlast * probabilities[i][j][size - 1];
      }
#endif
      //std::cout << *this << " " << up[i] << " " << i<< std::endl;
      ui += left_annot.inside_probabilities.array[j] * inner;
    }
    assert(ui >= 0.0);
    assert(ui <= 1.0);
  }

  // for(size_t i = 0; i < up.size(); ++i)
  //   std::cout << i << " b : " << up[i] << " ";
  // std::cout <<std::endl;
}


// only for  forest construction : suppose only *one* annotation !
void BRule::update_inside_annotations(std::vector<double>& up,
                                      const double& left_right_precomputation) const
{
  up[0] += probabilities[0][0][0] * left_right_precomputation;

    assert(up[0] >= 0.0);
    assert(up[0] <= 1.0);
}


#ifdef USE_THREADS
#include "utils/threads.h"
#endif

void BRule::update_outside_annotations(const AnnotationInfo& up_annot,
                                       AnnotationInfo& left_annot,
                                       AnnotationInfo& right_annot
                                       ) const
{
  #ifdef USE_THREADS
  std::vector<atomic<double>> & lo =
      *reinterpret_cast<std::vector<atomic<double>> *>(&left_annot.outside_probabilities.array);
  std::vector<atomic<double>> & ro =
      *reinterpret_cast<std::vector<atomic<double>> *>(&right_annot.outside_probabilities.array);
  #else
  std::vector<double> & lo = left_annot.outside_probabilities.array ;
  std::vector<double> & ro = right_annot.outside_probabilities.array ;
  #endif
  for(unsigned short i = 0; i < probabilities.size(); ++i) {
    if(up_annot.invalids[i] || up_annot.outside_probabilities.array[i] == 0.0) continue;
    const std::vector<std::vector<double> >& dim_i = probabilities[i];
    for(unsigned short j = 0; j < dim_i.size(); ++j) {
      const std::vector<double>& dim_j = dim_i[j];
      double temp4left = 0.0;
      double factor4right = 0.0;
      if(not left_annot.invalids[j]) factor4right = up_annot.outside_probabilities.array[i] * left_annot.inside_probabilities.array[j];
      for(unsigned short k = 0; k < dim_j.size(); ++k) {
        const double& t = dim_j[k];
        // if(right_in[k] != LorgConstants::NullProba) temp4left += right_in[k] * t;
        // if(right_out[k] != LorgConstants::NullProba) right_out[k] += factor4right * t;

        // I and O are always invalid at the same time
        if(not right_annot.invalids[k]) {
          temp4left += right_annot.inside_probabilities.array[k] * t;
          ro[k] += factor4right * t;
        }
      }
      if(not left_annot.invalids[j]) lo[j] += up_annot.outside_probabilities.array[i] * temp4left;
    }
  }
}

#include "utils/threads.h" // defines += operation on tbb::atomic<double>


double
BRule::update_outside_annotations_return_marginal(const AnnotationInfo& up_annot,
                                                  AnnotationInfo& left_annot,
                                                  AnnotationInfo& right_annot
                                                  ) const
{
  double marginal = 0.;
  #ifdef USE_THREADS
  auto & lo = *reinterpret_cast<std::vector<tbb::atomic<double>> *>(&left_annot.outside_probabilities.array);
  auto & ro = *reinterpret_cast<std::vector<tbb::atomic<double>> *>(&right_annot.outside_probabilities.array);
  #else
  auto & lo = left_annot.outside_probabilities.array ;
  auto & ro = right_annot.outside_probabilities.array ;
  #endif
  for(unsigned short i = 0; i < probabilities.size(); ++i) {
    if(up_annot.invalids[i] || up_annot.outside_probabilities.array[i] == 0.0) continue;
    const std::vector<std::vector<double> >& dim_i = probabilities[i];
    for(unsigned short j = 0; j < dim_i.size(); ++j) {
      const std::vector<double>& dim_j = dim_i[j];
      double temp4left = 0.0;
      double factor4right = 0.0;
      if(not left_annot.invalids[j]) factor4right = up_annot.outside_probabilities.array[i] * left_annot.inside_probabilities.array[j];
      for(unsigned short k = 0; k < dim_j.size(); ++k) {
        const double& t = dim_j[k];
        // if(right_in[k] != LorgConstants::NullProba) temp4left += right_in[k] * t;
        // if(right_out[k] != LorgConstants::NullProba) right_out[k] += factor4right * t;

        // I and O are always Null at the same time
        if(not right_annot.invalids[k]) {
          temp4left += right_annot.inside_probabilities.array[k] * t;
          ro[k] += factor4right * t;
        }
      }
      if(not left_annot.invalids[j]) {
        double delta_left = up_annot.outside_probabilities.array[i] * temp4left;
        lo[j] += delta_left;
        marginal += delta_left * left_annot.inside_probabilities.array[j];
      }
    }
  }
  return marginal ;
}

void BRule::remove_unlikely_annotations(const double& threshold)
{
  bool changed = false;

  for(unsigned i = 0; i < probabilities.size(); ++i)
    for(unsigned j = 0; j < probabilities[i].size(); ++j)
      for(unsigned k = 0; k < probabilities[i][j].size(); ++k)
	if(probabilities[i][j][k] < threshold) {
	  changed = true;
	  probabilities[i][j][k] = 0;
	}

  if(changed) compact();

}


void BRule::compact()
{
  // int nb_elements = 0;
  // int nb_zeros = 0;

  //get rid of lines of zeros
  for(unsigned i = 0; i < probabilities.size(); ++i)
    for(unsigned j = 0; j < probabilities[i].size(); ++j) {
      bool allzeros = true;
      for(unsigned k = 0; k < probabilities[i][j].size(); ++k) {
        //        ++nb_elements;
        if(probabilities[i][j][k] != 0.0) {
          allzeros = false;
          break;
        }
        // else
        // {
        //   ++nb_zeros;
        // }
      }
      if(allzeros)
      {
        //std::cerr << "compacted something of size " << probabilities[i][j].size() << std::endl;
        probabilities[i][j].clear();
        probabilities[i][j].shrink_to_fit();

      }
      //std::vector<double>(probabilities[i][j].begin(), probabilities[i][j].end()).swap(probabilities[i][j]);
    }

  //get rid of lines of empty vectors
  for(unsigned i = 0; i < probabilities.size(); ++i) {
    bool allempty = true;
    for(unsigned j = 0; j < probabilities[i].size(); ++j) {
      if(!probabilities[i][j].empty()) {
        allempty = false;
        break;
      }
    }
    if(allempty)
    {
      //std::cerr << "compacted 2 levels" << std::endl;
      probabilities[i].clear();
    }
    probabilities[i].shrink_to_fit();
    //std::vector< std::vector<double> >(probabilities[i].begin(), probabilities[i].end()).swap(probabilities[i]);
  }

  // std::cerr << "nb zeros: " <<nb_zeros
  //           << "nb_elements: " << nb_elements
  //           << " ratio: " << double(nb_zeros) *100 / double(nb_elements)
  //           << std::endl;

}


void BRule::uncompact(unsigned rhs0_size, unsigned rhs1_size)
{

  //  std::cout << "before uncompact\n" << *this << std::endl;

  for(unsigned i = 0; i < probabilities.size(); ++i) {
    if(probabilities[i].size() != rhs0_size)
      probabilities[i].resize(rhs0_size);
    for(unsigned j = 0; j < probabilities[i].size(); ++j) {
      if(probabilities[i][j].size() != rhs1_size)
        probabilities[i][j].resize(rhs1_size,0.0);
    }
  }

  //  std::cout << "after uncompact\n" << *this << std::endl;

}



bool BRule::operator==(const BRule& other) const
{
  return
    lhs  == other.lhs  &&
    rhs0 == other.rhs0 &&
    rhs1 == other.rhs1 &&
    probabilities == other.probabilities;
}

bool BRule::operator<(const BRule& other) const
{
  return
    lhs < other.lhs ||
    (
     lhs  == other.lhs  &&
     (rhs0 < other.rhs0 ||
      (
       rhs0 == other.rhs0 &&
       (rhs1 < other.rhs1 ||
	( rhs1 == other.rhs1 &&
	  probabilities < other.probabilities)))));
}

void BRule::set_probability(const double& value)
{
  for(unsigned i = 0 ; i < probabilities.size(); ++i)
    for(unsigned j = 0; j < probabilities[i].size(); ++j)
      for(unsigned k = 0; k < probabilities[i][j].size(); ++k)
	probabilities[i][j][k] = value;
}


// p_x = p(A_x -> B_y C_z), y and z fixed !!!
// p'_x = (1 - \alpha) p_x + \alpha (average_x' (p_x'))
// assume the rule is not compacted
void BRule::linear_smooth(const double& alpha)
{
  if(probabilities.size() == 1) return;

  double alphabar = 1 - alpha;
  double other_factor = alpha / (probabilities.size() - 1);

  std::vector< std::vector< std::vector<double> > >
    new_probabilities(probabilities.size(),
		      std::vector< std::vector<double> >(probabilities[0].size(),
							 std::vector<double>(probabilities[0][0].size(),0.0)));

  for(unsigned i = 0; i < probabilities.size(); ++i)
    for(unsigned o = 0; o < probabilities.size(); ++o)
      for(unsigned j = 0; j < probabilities[i].size(); ++j)
	for(unsigned k = 0; k < probabilities[i][j].size(); ++k) {
	  new_probabilities[i][j][k] += ( i == o ? alphabar : other_factor) * probabilities[o][j][k];
	}

  probabilities = new_probabilities;
}


// assume the rule is not compacted
void BRule::weighted_smooth(const double& alpha, const std::vector<std::vector<double> > & weights)
{
  if(probabilities.size() == 1) return;

  std::vector< std::vector< std::vector<double> > >
    new_probabilities(probabilities.size(),
		      std::vector< std::vector<double> >(probabilities[0].size(),
							 std::vector<double>(probabilities[0][0].size(),0.0)));

  const std::vector<double> & myweights = weights[lhs];

  for(unsigned i = 0; i < probabilities.size(); ++i) {

    double denominator = std::accumulate(myweights.begin(),myweights.end(), - myweights[i]);
    double alphabar = (1 - alpha) / myweights[i];
    double other_factor = alpha / denominator;


    for(unsigned o = 0; o < probabilities.size(); ++o)
      for(unsigned j = 0; j < probabilities[i].size(); ++j)
	for(unsigned k = 0; k < probabilities[i][j].size(); ++k) {
	  new_probabilities[i][j][k] += ( i == o ? alphabar : other_factor ) * probabilities[o][j][k] * myweights[o];
	}
  }
  probabilities = new_probabilities;
}


// assume the rule is not compacted
void BRule::generation_smooth(const std::vector<std::vector<std::vector<double> > > & weights)
{
  if(probabilities.size() == 1) return;

  //  std::cout << *this << std::endl;


  std::vector< std::vector< std::vector<double> > >
    new_probabilities(probabilities.size(),
		      std::vector< std::vector<double> >(probabilities[0].size(),
							 std::vector<double>(probabilities[0][0].size(),0.0)));

  const std::vector<std::vector<double> > & myweights = weights[lhs];

  for(unsigned i = 0; i < probabilities.size(); ++i) {
    for(unsigned o = 0; o < probabilities.size(); ++o)  {
      //  std::cout << myweights[i][o] << std::endl;
      for(unsigned j = 0; j < probabilities[i].size(); ++j) {
	//	std::cout << "before k" << std::endl;
	for(unsigned k = 0; k < probabilities[i][j].size(); ++k) {
	  //	  std::cout << "size myweights: " << myweights.size();
	  //	  std::cout << "i: " << i << " k: " << k << std::endl;
	  new_probabilities[i][j][k] += probabilities[o][j][k] * myweights[o][i];
	}
      }
    }
  }
  probabilities = new_probabilities;


  //  std::cout << *this << std::endl;

}
