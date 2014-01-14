#include <iostream>
#include "WapitiWrapper.h"


void wapiti_wrapper::set_file(const std::string& filename)
{
  this->file = fopen(filename.c_str(), "r");
}

void wapiti_wrapper::set_coefficient(int c)
{
  this->coefficient = c;
}


double wapiti_wrapper::crf_tag()
{
  //std::cerr << "crf_tag" << std::endl;

  this->raw = rdr_readraw(this->model->reader, this->file);
  seq_t * seq = rdr_raw2seq(this->model->reader, this->raw, false);

  const uint32_t T = seq->len;

  this->dual = dual_init(T, this->model->nlbl);

  uint32_t *out = (uint32_t*)xmalloc(sizeof(uint32_t) * T);
  double   *psc = (double*)xmalloc(sizeof(double) * T);


  //std::cerr << "before tag_viterbi" << std::endl;
  tag_viterbi(this->model, seq, (uint32_t*)out, &(this->score), psc,this->dual);
  //std::cerr << "after tag_viterbi" << std::endl;

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


  //auto save_score = this->score;

  uint32_t *out = (uint32_t*)xmalloc(sizeof(uint32_t) * T);
  double   *psc = (double*)xmalloc(sizeof(double) * T);
  //double   *scs = (double*)xmalloc(sizeof(double));
  tag_viterbi(this->model, seq, (uint32_t*)out, &(this->score), psc,this->dual);


  //auto save_string_sequence = best_string_sequence;
  best_string_sequence.clear();
  for (size_t i = 0; i < T; ++i)
  {
    best_string_sequence.push_back(qrk_id2str(this->model->reader->lbl, out[i]));
  }


  //free(scs);
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


        //std::cout << "update_relaxations1: " <<pair.second << std::endl;
        dual_add_binary_penalty(this->dual, pair.first + offset, i, j, -pair.second);
        //dual_add_unary_penalty(this->dual, pair.first, i, -pair.second);
        //dual_add_unary_penalty(this->dual, pair.first + offset, j, -pair.second);

        //        std::cout << "update_relaxations2: " <<dual_get_binary_penalty(this->dual, pair.first + offset, i,j) << std::endl;

      }
    }
  }

}


std::ostream& operator<<(std::ostream& out, const wapiti_wrapper& w)
{
  // for (size_t i = 0; i < best_string_sequence.size(); ++i)
  // {
  //   out << sentence[j].get_form() << "\t" << crfs[i].best_string_sequence[j] << std::endl;
  // }

  for(const auto& s : w.best_string_sequence)
  {
    out << s << " " ;
  }

  return out;
}
