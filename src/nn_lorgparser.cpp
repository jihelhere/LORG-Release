#include "NNLorgParseApp.h"
#include <cnn/cnn.h>

int main(int argc, char** argv)
{
  NNLorgParseApp app;


  cnn::initialize(argc, argv);
  if (!app.init(argc,argv)) return -1;

  app.run();

  return 0;
}
