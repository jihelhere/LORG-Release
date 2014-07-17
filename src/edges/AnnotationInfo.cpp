#include "AnnotationInfo.h"
#include "utils/LorgConstants.h"


scaled_array::scaled_array() : array(), scale(0) {}

scaled_array::scaled_array(unsigned i, double d) : array(i,d), scale(0) {}

scaled_array::scaled_array(const scaled_array& other) : array(other.array), scale(other.scale) {}
// {
//   scale = other.scale;
//   array.assign(other.array.begin(),other.array.end());
//   // array.resize(other.array.size());
//   // for(unsigned i = 0; i < other.array.size(); ++i)
//   //   array[i] = other.array[i];
// }

void scaled_array::scale_array(int previous_scale)
{
  //   std::cout << "scale" << std::endl;

  scale = 0;
  double scalingfactor = 1;

  double max = *(std::max_element(array.begin(),array.end()));

  while (max > SCALE) {
    //      std::cout << "loop1" << std::endl;
    max /= SCALE;
    scalingfactor *= SCALE;
    ++scale;
  }

  while (max > 0.0 && max < 1.0 / SCALE) {
    //      std::cout << "loop2" << std::endl;
    max *= SCALE;
    scalingfactor /= SCALE;
    --scale;
  }

  if(scale != 0) {
    for(unsigned i = 0; i < array.size(); ++i) {
      array[i] /= scalingfactor;
    }
  }

  scale += previous_scale;
}


double scaled_array::get_scaled_value(unsigned i) const
{
  return array[i]*std::pow(SCALE,scale);
}

double scaled_array::get_scaled_logvalue(unsigned i) const
{
  return std::log(array[i]) + LOGSCALE * scale;
}

void scaled_array::reset(double value)
{
  std::fill(array.begin(), array.end(), value);
  scale = 0;
}

void scaled_array::resize(unsigned new_size)
{
  array = std::vector<double>(new_size,0.0);
  scale = 0;
  //   std::cout << "scaled_array::resize("<<new_size<<") array at " << array.data() << std::endl; std::cout.flush();
}

//inline
double scaled_array::calculate_scalingfactor(int previous)
{
  switch (previous) {
  case 0 : return 1.0;
  case 1 : return SCALE;
  case 2 : return SCALE * SCALE;
  case -1: return 1.0/SCALE;
  case 3 : return SCALE * SCALE * SCALE; // very rare
    //   //the unseen cases in the PTB:
    // case -2:
    //   std::cout << previous << std::endl;
    //   return 1.0/SCALE/SCALE;
    // case -3:
    //   std::cout << previous << std::endl;
    //   return 1.0/SCALE/SCALE/SCALE;
  default:
    //      std::cout << previous << std::endl;
    return std::pow(SCALE,previous);

  }
}

double scaled_array::calculate_logscalingfactor(int previous)
{
  return LOGSCALE * previous;
}




////////////////////////////////


AnnotationInfo::AnnotationInfo() : inside_probabilities(), outside_probabilities(),
                                   invalids(),
                                   unary_temp()

{}

AnnotationInfo::AnnotationInfo(unsigned i, double d) : inside_probabilities(i,d), outside_probabilities(i,d),
                                                       invalids(i, false),
                                                       unary_temp(i,d)
{}

AnnotationInfo::AnnotationInfo(const AnnotationInfo& other)
  : inside_probabilities(other.inside_probabilities),
    outside_probabilities(other.outside_probabilities),
    invalids(other.invalids),
    unary_temp(other.unary_temp)
 {}


double AnnotationInfo::get_inside(int i) const
{
  return inside_probabilities.array[i];
}

double AnnotationInfo::get_outside(int i) const
{
  return outside_probabilities.array[i];
}

int AnnotationInfo::get_inside_scale() const
{
  return inside_probabilities.scale;
}

int AnnotationInfo::get_outside_scale() const
{
  return outside_probabilities.scale;
}

void AnnotationInfo::reset_inside_probabilities(double value, bool keepnullproba)
{
  inside_probabilities.reset(value);
  if (not keepnullproba)
    std::fill(invalids.begin(),invalids.end(),false);
}

void AnnotationInfo::reset_outside_probabilities(double value, bool keepnullproba)
{
  outside_probabilities.reset(value);
  if (not keepnullproba)
    std::fill(invalids.begin(),invalids.end(),false);
}

void AnnotationInfo::reset_probabilities( double value, bool keepnullproba)
{
  reset_inside_probabilities(value, keepnullproba);
  reset_outside_probabilities(value, keepnullproba);
}

unsigned AnnotationInfo::get_size() const
{
  return inside_probabilities.array.size(); //assume inside_probabilities and outside_probabilities have same size
}

void AnnotationInfo::resize(unsigned new_size)
{
  inside_probabilities.resize(new_size);
  outside_probabilities.resize(new_size);
  invalids.resize(new_size,false);
  unary_temp.resize(new_size);
}

bool AnnotationInfo::valid_prob_at(unsigned i) const
{
  return not invalids[i];
  // return (inside_probabilities.array[i] != invalid &&
  //         outside_probabilities.array[i] != invalid);
}
