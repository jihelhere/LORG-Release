// -*- mode: c++ -*-
#ifndef TWOSTAGELORGPARSEAPP_H
#define TWOSTAGELORGPARSEAPP_H

#include "LorgParseApp.h"
#include "ParseSolution.h"

#include "utils/WapitiWrapper.h"


class ParserCKYAll;

class TwoStageLorgParseApp : public LorgParseApp
{
public:
  TwoStageLorgParseApp();
  ~TwoStageLorgParseApp();
  int run();

private:
  bool read_config(ConfigTable& configuration);
  LorgOptions get_options() const;

  int find_consensus(std::vector<std::pair<PtbPsTree*,double>>&);

  std::vector<std::string> crf_tag(FILE* f, int idx);


  std::vector<ParserCKYAll *> parsers;
  std::vector<Tagger> taggers;


  std::vector<wapiti_wrapper> crfs;

  bool output_annotations;
  unsigned kbest;

  bool extract_features;
  parse_solution::parse_solution_format output_format;
};

#endif // TWOSTAGELORGPARSEAPP_H
