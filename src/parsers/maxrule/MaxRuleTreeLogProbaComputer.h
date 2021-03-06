// // -*- mode: c++ -*-
#ifndef _MAXRULEUPDATER_H_
#define _MAXRULEUPDATER_H_


#include "PCKYAllCell.h"

#include <numeric>
#include <math.h>

#include "parsers/ParserCKYAllFactory.h"

#include <xmmintrin.h>  // Need this for SSE compiler intrinsics

template<class ProbaModel>
class MaxRuleTreeLogProbaComputer
{
public:
  typedef typename ProbaModel::Edge Edge ;
  typedef typename Edge::LexicalDaughter LexicalDaughter;
  typedef typename Edge::UnaryDaughter UnaryDaughter;
  typedef typename Edge::BinaryDaughter BinaryDaughter;

  static  ParserCKYAllFactory::MaxParsing_Calculation calculation;

 static void set_calculation(ParserCKYAllFactory::MaxParsing_Calculation c) {calculation = c;}



  static double compute (const AnnotationInfo & up_annotations,
                         const BinaryDaughter & dtr,
                         double normalisation_factor,
                         unsigned left_idx =0, unsigned right_idx = 0)
  {

    // if (normalisation_factor == -std::numeric_limits<double>::infinity())
    //   return normalisation_factor;


    //    std::cout << left_idx << " : " << right_idx << std::endl;
    //    std::cout << *(dtr.get_rule()) << std::endl;

    const Edge & left  = dtr.left_daughter ();
    const Edge & right = dtr.right_daughter ();

    double probability = 0.0;

    const scaled_array & left_inside  = left.get_annotations ().inside_probabilities;
    const scaled_array & right_inside = right.get_annotations ().inside_probabilities;
    const scaled_array & up_outside = up_annotations.outside_probabilities;
    const auto & rule_probs = dtr.get_rule ()->get_probability ();


    // size_t size = right_inside.array.size();
    // size_t N = size / 2;
    // bool odd = (size % 2) != 0  and not right.get_annotations().invalids[size-1];
    // auto& rlast = right_inside.array[size -1];


    unsigned rsize = rule_probs.size ();
    for (unsigned i = 0; i < rsize; ++i)
    {
      if (not up_annotations.valid_prob_at (i))
        continue;
      double temp = 0;
      const auto& rule_probs_i = rule_probs[i];
      unsigned size_i = rule_probs_i.size ();
      for (unsigned j = 0; j < size_i; ++j)
      {
        if (not left.valid_prob_at (j))
          continue;
        double inner = 0;
        const auto& rule_probs_ij = rule_probs_i[j];
        unsigned size_ij = rule_probs_ij.size ();

        // if (size_ij == 0) continue;

        // __m128d * a = (__m128d*) right_inside.array.data();
        // __m128d * p = (__m128d*) rule_probs_ij.data();
        // double tmp_double[2] = {0.0,0.0};
        // __m128d * tmp = (__m128d*) tmp_double;

        // for (size_t k = 0; k < N; ++k, ++a, ++p)
        // {
        //   //std::cerr << "k " << k << std::endl;
        //   _mm_store_pd(tmp_double,
        //                _mm_add_pd( *tmp,
        //                            _mm_mul_pd(*a,*p)));
        // }

        // double inner = tmp_double[0] + tmp_double[1];
        // if (odd)
        // {
        //   inner += rlast * rule_probs_ij[size - 1];
        // }



        for (unsigned k = 0; k < size_ij; ++k)
        {
          if (right.valid_prob_at (k))
            inner += rule_probs_ij[k] * right_inside.array[k];
        }


      temp += left_inside.array[j] * inner;

      }
      probability += up_outside.array[i] * temp;
    }

   if (calculation == ParserCKYAllFactory::Product)
   {
     double logprod = (std::log (probability) - normalisation_factor);
     return  logprod
        + left.get_prob_model().get(left_idx).probability
        + right.get_prob_model().get(right_idx).probability;
   }

   if (calculation == ParserCKYAllFactory::Sum)
   {
     double logprod = (std::log (probability) - normalisation_factor);
     double sum = std::exp(logprod);
     return sum
        + left.get_prob_model().get(left_idx).probability
        + right.get_prob_model().get(right_idx).probability;
   }

   if (calculation == ParserCKYAllFactory::ProdSum)
   {
     double logprod = (std::log (probability) - normalisation_factor);
     double sum = std::exp(logprod);
     return logprod + std::pow(sum, 4.0)
        + left.get_prob_model().get(left_idx).probability
        + right.get_prob_model().get(right_idx).probability;
   }

   return -std::numeric_limits < double >::infinity ();

  }

  static void compute_best_indexes (const AnnotationInfo & up_annotations,
                                    const BinaryDaughter & dtrs,
                                    double normalisation_factor,
                                    unsigned &left_idx, unsigned &right_idx)
  {
    const Edge & left  = dtrs.left_daughter ()->get_edge (dtrs.get_rule ()->get_rhs0 ());
    const Edge & right = dtrs.right_daughter ()->get_edge (dtrs.get_rule ()->get_rhs1 ());

    const scaled_array & left_inside  = left.get_annotations ().inside_probabilities;
    const scaled_array & right_inside = right.get_annotations ().inside_probabilities;
    const scaled_array & up_outside = up_annotations.outside_probabilities;
    const std::vector < std::vector < std::vector < double >>>&rule_probs =
    dtrs.get_rule ()->get_probability ();

    double max_prob = -std::numeric_limits < double >::infinity ();
    double base = - normalisation_factor
                  +  left.get_prob_model().get (0).probability
                  + right.get_prob_model().get (0).probability;

    for (unsigned i = 0; i < rule_probs.size (); ++i)
    {
      if (not up_annotations.valid_prob_at (i))
        continue;
      const std::vector < std::vector < double >>&rule_probs_i =
      rule_probs[i];

      double contrib_i = std::log (up_outside.array[i]) + base;

      for (unsigned j = 0; j < rule_probs_i.size (); ++j)
      {
        if (!left.valid_prob_at (j))
          continue;
        const std::vector < double >&rule_probs_ij = rule_probs_i[j];

        double contrib_j = contrib_i + std::log (left_inside.array[j]);

        for (unsigned k = 0; k < rule_probs_ij.size (); ++k)
        {
          if (!right.valid_prob_at (k))
            continue;
          double contrib = std::log(right_inside.array[k]) + std::log(rule_probs_ij[k]) + contrib_j;

          if (contrib > max_prob)
          {
            left_idx = j;
            right_idx = k;
            max_prob = contrib;
          }
        }
      }
    }
  }

  static double compute (const AnnotationInfo & up_annotations,
                                            const UnaryDaughter & dtrs,
                                            double normalisation_factor,
                                            unsigned left_idx = 0)
  {
    // if (normalisation_factor == -std::numeric_limits<double>::infinity())
    //   return normalisation_factor;


    double probability = 0;

    const Edge & left = dtrs.left_daughter();

    const scaled_array & left_inside = left.get_annotations ().inside_probabilities;
    const scaled_array & up_outside = up_annotations.outside_probabilities;
    const std::vector < std::vector < double >>&rule_probs = dtrs.get_rule ()->get_probability ();

    //FIXME: this should not happen because chart is clean ???
    // only relevant in kmax parsing
    // if (left.get_prob_model().n_deriv () == 0)
    //   return -std::numeric_limits < double >::infinity ();


    for (unsigned i = 0; i < rule_probs.size (); ++i)
    {
      if (not up_annotations.valid_prob_at (i))
        continue;
      double inner (0.0);

      const auto& probs_i = rule_probs[i];

      for (unsigned j = 0; j < probs_i.size (); ++j)
      {
        if (not left.valid_prob_at (j))
          continue;
        inner += probs_i[j] * left_inside.array[j];
      }
      probability += up_outside.array[i] * inner;
    }




   if (calculation == ParserCKYAllFactory::Product)
   {
     double logprod = (std::log (probability) - normalisation_factor);
     return  logprod
         + left.get_prob_model().get(left_idx).probability;
   }

   if (calculation == ParserCKYAllFactory::Sum)
   {
     double logprod = (std::log (probability) - normalisation_factor);
     double sum = std::exp(logprod);
     return sum
         + left.get_prob_model().get(left_idx).probability;
   }

   if (calculation == ParserCKYAllFactory::ProdSum)
   {
     double logprod = (std::log (probability) - normalisation_factor);
     double sum = std::exp(logprod);
     return logprod + std::pow(sum, 4.0)
         + left.get_prob_model().get(left_idx).probability;
   }

   return -std::numeric_limits < double >::infinity ();
  }

  static void compute_best_indexes (const AnnotationInfo & up_annotations,
                                    const UnaryDaughter & dtrs,
                                    double normalisation_factor,
                                    unsigned &left_idx)
  {
    const Edge & left =
    dtrs.left_daughter ()->get_edge (dtrs.get_rule ()->get_rhs0 ());

    const scaled_array & left_inside = left.get_annotations ().inside_probabilities;
    const scaled_array & up_outside = up_annotations.outside_probabilities;
    const std::vector < std::vector < double >>&rule_probs = dtrs.get_rule ()->get_probability ();

    double max_prob = -std::numeric_limits < double >::infinity ();

    double base = - normalisation_factor
                  + left.get_prob_model().get (0).probability;


    for (unsigned i = 0; i < rule_probs.size (); ++i)
    {
      if (not up_annotations.valid_prob_at (i))
        continue;
      const std::vector < double >&rule_probs_i = rule_probs[i];
      double contrib_i = std::log (up_outside.array[i]) + base;
      for (unsigned j = 0; j < rule_probs_i.size (); ++j)
      {
        if (not left.valid_prob_at (j))
          continue;
        double contrib =
        std::log (left_inside.array[j]) + std::log (rule_probs_i[j]) +
        contrib_i;
        if (contrib > max_prob)
        {
          left_idx = j;
          max_prob = contrib;
        }
      }
    }
  }

  static double compute (const AnnotationInfo & up_annotations,
                         const LexicalRuleC2f * rule_ptr,
                         double normalisation_factor)
  {
    // if (normalisation_factor == -std::numeric_limits<double>::infinity())
    //   return normalisation_factor;


    double probability = 0.0;

    const auto& probs = rule_ptr->get_probability ();

    for (unsigned i = 0; i < probs.size (); ++i)
    {
      if (not up_annotations.valid_prob_at (i))
        continue;
      //std::cout << up_annotations.outside_probabilities.array[i] << std::endl;
        probability +=
        probs[i] *
        up_annotations.outside_probabilities.array[i];
    }


   if (calculation == ParserCKYAllFactory::Product)
   {
     double logprod = (std::log (probability) - normalisation_factor);
     return  logprod;
   }

   if (calculation == ParserCKYAllFactory::Sum)
   {
     double logprod = (std::log (probability) - normalisation_factor);
     double sum = std::exp(logprod);
     return sum;
   }

   if (calculation == ParserCKYAllFactory::ProdSum)
   {
     double logprod = (std::log (probability) - normalisation_factor);
     double sum = std::exp(logprod);
     return logprod + std::pow(sum, 4.0);
   }

   return -std::numeric_limits < double >::infinity ();
  }

  static  double compute_simple (const AnnotationInfo& up_annotations,
                                 double normalisation_factor,
                                 const AnnotationInfo& left_annotations,
                                 const AnnotationInfo& right_annotations,
                                 const std::vector<std::vector<std::vector<double>>>& rule_probs)
  {
    double probability = 0.0;

    const scaled_array & left_inside  = left_annotations.inside_probabilities;
    const scaled_array & right_inside = right_annotations.inside_probabilities;
    const scaled_array & up_outside   = up_annotations.outside_probabilities;

#ifdef USE_MANUAL_SSE
    size_t size = right_inside.array.size();
    size_t N = size / 2;
    bool odd = (size % 2) != 0  and not right_annotations.invalids[size-1];
    auto& rlast = right_inside.array[size -1];
#endif

    unsigned rsize = rule_probs.size ();
    for (unsigned i = 0; i < rsize; ++i)
    {
      if (not up_annotations.valid_prob_at (i)) continue;
      double temp = 0;
      const auto& rule_probs_i = rule_probs[i];
      unsigned size_i = rule_probs_i.size ();
      for (unsigned j = 0; j < size_i; ++j)
      {
        if (not left_annotations.valid_prob_at (j))  continue;
        double inner = 0;
        const auto& rule_probs_ij = rule_probs_i[j];
        unsigned size_ij = rule_probs_ij.size ();

#ifdef USE_MANUAL_SSE
        if (size_ij == 0) continue;

        __m128d * a = (__m128d*) right_inside.array.data();
        __m128d * p = (__m128d*) rule_probs_ij.data();
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
          inner += rlast * rule_probs_ij[size - 1];
        }
#else
        for (unsigned k = 0; k < size_ij; ++k)
        {
          if (right_annotations.valid_prob_at (k))
            inner += rule_probs_ij[k] * right_inside.array[k];
          //          std::cout << "inner " << inner << std::endl;
        }
#endif
        temp += left_inside.array[j] * inner;
      }
      probability += up_outside.array[i] * temp;
      //      std::cout << probability << std::endl;
    }

   if (calculation == ParserCKYAllFactory::Product)
   {
     double logprod = (std::log (probability) - normalisation_factor);
     return  logprod;
   }

   if (calculation == ParserCKYAllFactory::Sum)
   {
     double logprod = (std::log (probability) - normalisation_factor);
     double sum = std::exp(logprod);
     return sum;
   }

   if (calculation == ParserCKYAllFactory::ProdSum)
   {
     double logprod = (std::log (probability) - normalisation_factor);
     double sum = std::exp(logprod);
     return logprod + std::pow(sum, 4.0);
   }

   return -std::numeric_limits < double >::infinity ();
  }

  static double compute_simple (const AnnotationInfo & up_annotations,
                                double normalisation_factor,
                                const AnnotationInfo &
                                left_annotations,
                                const std::vector<std::vector<double>>&rule_probs)
  {
    double probability = 0;

    const scaled_array & left_inside = left_annotations.inside_probabilities;
    const scaled_array & up_outside = up_annotations.outside_probabilities;

    for (unsigned i = 0; i < rule_probs.size (); ++i)
    {
      if (not up_annotations.valid_prob_at (i))
        continue;
      double inner (0.0);
      const auto& rule_probs_i = rule_probs[i];
      for (unsigned j = 0; j < rule_probs_i.size (); ++j)
      {
        if (not left_annotations.valid_prob_at (j))
          continue;
        inner += rule_probs_i[j] * left_inside.array[j];
      }
      probability += up_outside.array[i] * inner;
    }

   if (calculation == ParserCKYAllFactory::Product)
   {
     double logprod = (std::log (probability) - normalisation_factor);
     return  logprod;
   }

   if (calculation == ParserCKYAllFactory::Sum)
   {
     double logprod = (std::log (probability) - normalisation_factor);
     double sum = std::exp(logprod);
     return sum;
   }

   if (calculation == ParserCKYAllFactory::ProdSum)
   {
     double logprod = (std::log (probability) - normalisation_factor);
     double sum = std::exp(logprod);
     return logprod + std::pow(sum, 4.0);
   }

   return -std::numeric_limits < double >::infinity ();
  }
};

#endif /* _MAXRULEPROBABILITY_H_ */
