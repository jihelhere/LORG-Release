// -*- mode: c++ -*-
#pragma once
#include "utils/ConfigTable.h"

#ifdef USE_THREADS
#include <tbb/task_scheduler_init.h>
#endif

class LorgApp
{
public:
  LorgApp();

  bool init(int argc, char **argv);
  virtual int run() = 0;

  virtual ~LorgApp();

protected:
  virtual bool read_config(ConfigTable& configuration);
  virtual LorgOptions get_options()const =0;

protected:
  bool verbose; // maybe use an integer to have several verbosity levels ?

  std::istream* in;
  std::ostream* out;

  #ifdef USE_THREADS
  unsigned nbthreads;
  //  tbb::task_scheduler_init tbb_task_scheduler;
  #endif

};
