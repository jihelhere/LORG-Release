// -*- mode: c++ -*-
#pragma once

extern "C" {
#include "wapiti/src/decoder.h"
#include "wapiti/src/tools.h"
}

#include <string>
#include <vector>


struct wapiti_wrapper
{
  mdl_t * model;
  FILE * file;

  raw_t *raw;

  wapiti_wrapper(const std::string& modelfilename)
      : model(nullptr), file(nullptr), raw(nullptr)
  {
    this->model = mdl_new(rdr_new(false));

    FILE * modelfile = fopen(modelfilename.c_str(), "r");
    mdl_load(this->model, modelfile);
  }

  void set_file(const std::string& filename);
  std::vector<std::string> crf_tag();
  std::vector<std::string> redo_crf_tag();




  ~wapiti_wrapper(){}

  void clean_sentence()
  {
    rdr_freeraw(this->raw);
  }

  void clean()
  {
    mdl_free(this->model);
  }

};
