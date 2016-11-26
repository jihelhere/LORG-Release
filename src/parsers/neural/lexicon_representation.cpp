#include "lexicon_representation.h"

std::vector<dynet::LSTMBuilder> bilstm_lexicon_representation::word_l2r_builders;
std::vector<dynet::LSTMBuilder> bilstm_lexicon_representation::word_r2l_builders;


void
simple_lexicon_representation::precompute_embeddings(const std::vector<Word>& words)
{
  embeddings.clear();

  std::lock_guard<std::mutex> guard(computation_attachment::cg_mutex);

  for (const auto& w : words)
    embeddings.push_back(wr->word_reprentation(w));
}

/////////

bilstm_lexicon_representation::bilstm_lexicon_representation(
    word_representation* word_repr,
    bool init_global,
    dynet::Model& m,
    unsigned word_embedding_size,
    unsigned lstm_hidden_size,
    unsigned l)
    : lexicon_representation(word_repr),
      computation_attachment(),
      layers(l)
{
  if (init_global)
  {
    if (layers == 1)
    {
      word_l2r_builders = std::vector<dynet::LSTMBuilder>(1,dynet::LSTMBuilder(2, word_embedding_size, lstm_hidden_size, &m));
      word_r2l_builders = std::vector<dynet::LSTMBuilder>(1,dynet::LSTMBuilder(2, word_embedding_size, lstm_hidden_size, &m));
      }
    else
    {
      word_l2r_builders = std::vector<dynet::LSTMBuilder>(1,dynet::LSTMBuilder(1, word_embedding_size, lstm_hidden_size, &m));
      word_r2l_builders = std::vector<dynet::LSTMBuilder>(1,dynet::LSTMBuilder(1, word_embedding_size, lstm_hidden_size, &m));

      for (unsigned s = 1; s < layers; ++s)
      {
        word_l2r_builders.push_back(dynet::LSTMBuilder(1, lstm_hidden_size, lstm_hidden_size, &m));
        word_r2l_builders.push_back(dynet::LSTMBuilder(1, lstm_hidden_size, lstm_hidden_size, &m));
      }
    }
  }
}


void bilstm_lexicon_representation::precompute_embeddings(
    const std::vector<Word>& words)
{
  embeddings.clear();

  std::lock_guard<std::mutex> guard(computation_attachment::cg_mutex);

  for (const auto& w : words)
    embeddings.push_back(wr->word_reprentation(w));

  for (unsigned s = 0; s < layers; ++s)
  {
    std::vector<dynet::expr::Expression> lstm_forward, lstm_backward;

    // Build forward LSTM
    word_l2r_builders[s].new_graph(*cg);
    word_l2r_builders[s].start_new_sequence();
    for (const auto& input : embeddings)
    {
      lstm_forward.push_back(word_l2r_builders[s].add_input(input));
    }

    // Build backward LSTM
    word_r2l_builders[s].new_graph(*cg);
    word_r2l_builders[s].start_new_sequence();
    for (int i = embeddings.size() -1; i >= 0; --i)
    {
      lstm_backward.push_back(word_r2l_builders[s].add_input(embeddings[i]));
    }

    for (unsigned int i = 0 ; i < lstm_forward.size() ; ++i)
    {
      auto&& e = lstm_backward[lstm_backward.size() - i - 1] +
                 lstm_forward[i];
      embeddings[i] = e;
    }
  }
}

void bilstm_lexicon_representation::set_dropout(float d)
{
  for (auto& b :  word_l2r_builders) {
    b.set_dropout(d);
  }
  for (auto& b :  word_r2l_builders) {
    b.set_dropout(d);
  }
}

void bilstm_lexicon_representation::unset_dropout()
{
  for (auto& b :  word_l2r_builders) {
    b.disable_dropout();
  }
  for (auto& b :  word_r2l_builders) {
    b.disable_dropout();
  }
}
