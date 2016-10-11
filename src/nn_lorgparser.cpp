#include "NNLorgParseApp.h"
#include <dynet/dynet.h>

int main(int argc, char** argv)
{
  NNLorgParseApp app;


  dynet::initialize(argc, argv);
  if (!app.init(argc,argv)) return -1;

  app.run();

  return 0;
}
