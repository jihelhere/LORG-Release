// -*- mode: c++ -*-
#pragma once

extern "C" {
#include "wapiti/src/decoder.h"
#include "wapiti/src/tools.h"
}

#include <string>
#include <vector>
#include <unordered_map>


struct wapiti_wrapper
{
  mdl_t * model;
  FILE * file;
  double coefficient;


  raw_t *raw;
  dual_t* dual;
  double score;
  std::vector<std::string> best_string_sequence;


  wapiti_wrapper(const std::string& modelfilename, int coef = 1)
      : model(nullptr), file(nullptr), coefficient(coef), raw(nullptr), dual(nullptr)
  {
    this->model = mdl_new(rdr_new(false));

    FILE * modelfile = fopen(modelfilename.c_str(), "r");
    mdl_load(this->model, modelfile);
  }

  void set_file(const std::string& filename);
  void set_coefficient(int c);

  double crf_tag();
  double crf_retag();

  void update_relaxations(const std::unordered_map<unsigned,double>& lambdas, char first, char second, int offset);

  ~wapiti_wrapper(){}

  void clean_sentence()
  {
    rdr_freeraw(this->raw);
    dual_free(this->dual);
  }

  void clean()
  {
    mdl_free(this->model);
  }

  friend std::ostream& operator<<(std::ostream& out, const wapiti_wrapper& w);


};
