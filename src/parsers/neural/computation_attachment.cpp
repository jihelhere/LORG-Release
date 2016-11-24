#include "computation_attachment.h"



std::mutex computation_attachment::cg_mutex;
dynet::ComputationGraph * computation_attachment::cg = nullptr;

computation_attachment::computation_attachment() {}

void computation_attachment::set_cg(dynet::ComputationGraph * g)
{
  cg = g;
}


dynet::ComputationGraph * computation_attachment::get_cg()
{
  return cg;
}
