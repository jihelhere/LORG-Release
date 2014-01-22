// -*- mode: c++ -*-
#pragma once

#include "LorgParseApp.h"
#include "utils/WapitiWrapper.h"



class WapitiDDApp : public LorgParseApp
{
public:
  WapitiDDApp();
  ~WapitiDDApp();
  int run();

private:
  bool read_config(ConfigTable& configuration);
  LorgOptions get_options() const;

  int find_consensus();
  std::vector<std::string> crf_tag(FILE* f, int idx);


  std::vector<wapiti_wrapper_bi> crfs;
};
