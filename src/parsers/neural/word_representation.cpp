#include "word_representation.h"


dynet::LookupParameter simple_word_representation::_p_word;
dynet::LookupParameter char_word_representation::_p_word;

dynet::LSTMBuilder char_word_representation::letter_l2r_builder;
dynet::LSTMBuilder char_word_representation::letter_r2l_builder;




simple_word_representation::simple_word_representation(bool init_global,
                                                       dynet::Model& m,
                                                       unsigned word_embedding_size)
    : computation_attachment()
{
  if (init_global)
  {
    unsigned s = SymbolTable::instance_word().get_symbol_count() + 1;
    _p_word = m.add_lookup_parameters(s, {word_embedding_size});
  }
}

dynet::expr::Expression simple_word_representation::word_reprentation(const Word& w)
{
  return dynet::expr::lookup(*cg,_p_word, w.get_id());
}

dynet::expr::Expression simple_word_representation::get_pad_expression()
{
  static int pad = SymbolTable::instance_word().get_symbol_count();
  return dynet::expr::lookup(*cg, _p_word, pad);
}

///////////

char_word_representation::char_word_representation(bool init_global,
                                                   dynet::Model& m,
                                                   unsigned word_embedding_size)
    : computation_attachment()
{
  if (init_global)
  {
    unsigned s = SymbolTable::instance_word().get_symbol_count() + 1;
    _p_word = m.add_lookup_parameters(s, {word_embedding_size});

    letter_l2r_builder = dynet::LSTMBuilder(1, 1, word_embedding_size, &m);
    letter_r2l_builder = dynet::LSTMBuilder(1, 1, word_embedding_size, &m);
  }
}

dynet::expr::Expression char_word_representation::word_reprentation(const Word& w)
{
  // Build forward LSTM
  letter_l2r_builder.new_graph(*cg);
  letter_l2r_builder.start_new_sequence();
  // initialize with <w>
  letter_l2r_builder.add_input(dynet::expr::input(*cg, 0.0));
  auto&& form = w.get_form();
  for (unsigned c = 0; c < form.size(); ++c)
  {
    letter_l2r_builder.add_input(dynet::expr::input(*cg, float(form[c])));
  }

  // Build backward LSTM
  letter_r2l_builder.new_graph(*cg);
  letter_r2l_builder.start_new_sequence();
  // initialize with <w> (encoded as zero)
  letter_r2l_builder.add_input(dynet::expr::input(*cg,0.0));
  for (int i = form.size() - 1; i >= 0; --i)
  {
    letter_r2l_builder.add_input(dynet::expr::input(*cg,float(form[i])));
  }

  // finalize with </w> and concatenate (encoded as one)
  return
      dynet::expr::lookup(*cg,_p_word, w.get_id()) +
      letter_l2r_builder.add_input(dynet::expr::input(*cg,1.0)) +
      letter_r2l_builder.add_input(dynet::expr::input(*cg,1.0))
      ;
}


dynet::expr::Expression char_word_representation::get_pad_expression()
{
  static int pad = SymbolTable::instance_word().get_symbol_count();
  return dynet::expr::lookup(*cg, _p_word, pad);
}
