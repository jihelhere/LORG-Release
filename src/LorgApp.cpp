#include "LorgApp.h"

#include "utils/SymbolTable.h"



LorgApp::LorgApp() : verbose(false), in(NULL), out(NULL), in_filename()
#ifdef USE_THREADS
                     //, tbb_task_scheduler(tbb::task_scheduler_init::deferred)
#endif
{}

LorgApp::~LorgApp()
{
  if (in != &std::cin) delete in;
  if (out != &std::cout) delete out;
}

bool LorgApp::init(int argc, char **argv)
{

  ConfigTable& configuration =  ConfigTable::create(argc,argv,get_options());

  bool res = read_config(configuration);

#ifdef USE_THREADS
  nbthreads = configuration.get_value<unsigned>("nbthreads");
  nbthreads = nbthreads ? nbthreads : tbb::task_scheduler_init::default_num_threads();
#endif

  return res;
}



bool LorgApp::read_config(ConfigTable& configuration)
{
  verbose = configuration.exists("verbose");
  // parse config file if provided
  if(configuration.exists("config-file"))
    {
      if(verbose) std::clog << "Parsing configuration file." << std::endl;
      configuration.parse_config_file(configuration.get_value<std::string>("config-file"));
    }


  if (configuration.exists("help")) {
    configuration.print_help();
    return false;
  }

  return true;
}
