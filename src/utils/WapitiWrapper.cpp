#include <iostream>
#include "WapitiWrapper.h"





void wapiti_wrapper::set_file(const std::string& filename)
{
  this->file = fopen(filename.c_str(), "r");
}


double wapiti_wrapper::crf_tag()
{
  this->raw = rdr_readraw(this->model->reader, this->file);
  seq_t * seq = rdr_raw2seq(this->model->reader, this->raw, false);

  const uint32_t T = seq->len;

  this->dual = dual_init(T, this->model->nlbl);

  uint32_t *out = (uint32_t*)xmalloc(sizeof(uint32_t) * T);
  double   *psc = (double*)xmalloc(sizeof(double) * T);

  tag_viterbi(this->model, seq, (uint32_t*)out, &(this->score), psc,this->dual);

  std::vector<std::string> result(T);

  best_string_sequence.clear();
  for (size_t i = 0; i < T; ++i)
  {
    best_string_sequence.push_back(qrk_id2str(this->model->reader->lbl, out[i]));
  }

  free(psc);
  free(out);
  rdr_freeseq(seq);
  //  rdr_freeraw(this->raw);
  return this->score;
}


double wapiti_wrapper::crf_retag()
{
  //this->raw = rdr_readraw(this->model->reader, this->file);
  seq_t * seq = rdr_raw2seq(this->model->reader, this->raw, false);

  const uint32_t T = seq->len;

  //this->dual = dual_init(T, this->model->nlbl);

  uint32_t *out = (uint32_t*)xmalloc(sizeof(uint32_t) * T);
  double   *psc = (double*)xmalloc(sizeof(double) * T);
  double   *scs = (double*)xmalloc(sizeof(double));
  tag_viterbi(this->model, seq, (uint32_t*)out, scs, psc,this->dual);

  best_string_sequence.clear();
  for (size_t i = 0; i < T; ++i)
  {
    best_string_sequence.push_back(qrk_id2str(this->model->reader->lbl, out[i]));
  }


  free(scs);
  free(psc);
  free(out);
  rdr_freeseq(seq);
  //rdr_freeraw(this->raw);
  return this->score;
}


void wapiti_wrapper::update_relaxations(const std::unordered_map<unsigned,double>& lambdas, char first, char second, int offset)
{
  unsigned cpt = qrk_count(this->model->reader->lbl);

  for (const auto& pair : lambdas)
  {
    for (size_t i = 0; i < cpt; ++i)
    {
      if ((first != '*') and (qrk_id2str(this->model->reader->lbl, i)[1] != first))
        continue;

      for (size_t j = 0; j < cpt; ++j)
      {
        if ((second != '*') and (qrk_id2str(this->model->reader->lbl, j)[1] != second))
          continue;

        dual_add_binary_penalty(this->dual, pair.first + offset, i, j, -pair.second);

      }
    }
  }

}
