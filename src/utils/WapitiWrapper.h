// -*- mode: c++ -*-
#pragma once

extern "C" {
#include "wapiti/src/decoder.h"
#include "wapiti/src/tools.h"
}

struct wapiti_wrapper
{
  mdl_t * model;
  FILE * file;

  raw_t *raw;
  seq_t *seq;

  wapiti_wrapper(const std::string& model_name)
      : model(nullptr), file(nullptr), raw(nullptr), seq(nullptr)
  {};

  void set_file(const std::string& filename)
  {
    file = fopen(filename.c_str(), "r");
  }


  std::vector<std::string> crf_tag(FILE* file, int idx)
  {
    this->raw = rdr_readraw(this->mdl->reader, this->file);
    this->seq = rdr_raw2seq(this->mdl->reader, this->raw, false);

    const uint32_t T = seq->len;

    dual_t* d = dual_init(T, this->mdl->nlbl);

    uint32_t *out = (uint32_t*)xmalloc(sizeof(uint32_t) * T);
    double   *psc = (double*)xmalloc(sizeof(double) * T);
    double   *scs = (double*)xmalloc(sizeof(double));
    tag_viterbi(this->mdl, this->seq, (uint32_t*)out, scs, (double*)psc,d);

    std::vector<std::string> result(T);

    for (size_t i = 0; i < T; ++i)
    {
      result[i] = qrk_id2str(this->mdl->reader->lbl, out[i]);
    }


   free(scs);
   free(psc);
   free(out);
   rdr_freeseq(seq);
   rdr_freeraw(raw);
   return result;
  }


};
