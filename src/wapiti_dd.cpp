#include "WapitiDDApp.h"


int main(int argc, char** argv)
{
  WapitiDDApp app;

  if (!app.init(argc,argv)) {
    std::cerr << "Initialisation failed" << std::endl;
    return -1;
  }
  else
    app.run();

  return 0;
}
