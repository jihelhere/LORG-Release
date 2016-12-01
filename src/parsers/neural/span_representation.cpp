#include "span_representation.h"


dynet::Parameter all_span_representation::_p_W_span_distance;
dynet::Parameter all_span_representation::_p_W_span_extra;

dynet::Parameter all_span_representation::_p_b_span_bin;
dynet::Parameter all_span_representation::_p_o_span_bin;

dynet::Parameter all_span_representation::_p_b_span_un;
dynet::Parameter all_span_representation::_p_o_span_un;

dynet::Parameter all_span_representation::_p_W_span_mid;

dynet::Parameter all_span_representation::_p_b_span_end;
dynet::Parameter all_span_representation::_p_o_span_end;
dynet::Parameter all_span_representation::_p_b_span_end_un;
dynet::Parameter all_span_representation::_p_o_span_end_un;

//dynet::Parameter all_span_representation::_p_W_span_init;
dynet::Parameter all_span_representation::_p_W_span_init_n;
dynet::Parameter all_span_representation::_p_W_span_init_w;

dynet::Parameter all_span_representation::_p_b_span_init;
dynet::Parameter all_span_representation::_p_o_span_init;
dynet::Parameter all_span_representation::_p_b_span_init_un;
dynet::Parameter all_span_representation::_p_o_span_init_un;

dynet::Parameter all_span_representation::_p_W_span_left;
dynet::Parameter all_span_representation::_p_W_span_right;

dynet::Parameter all_span_representation::_p_b_span_split;
dynet::Parameter all_span_representation::_p_o_span_split;

//dynet::Parameter all_span_representation::_p_W_span0_init;
dynet::Parameter all_span_representation::_p_W_span0_init_n;
dynet::Parameter all_span_representation::_p_W_span0_init_w;
dynet::Parameter all_span_representation::_p_b_span0_init;
dynet::Parameter all_span_representation::_p_o_span0_init;
dynet::Parameter all_span_representation::_p_b_span0_init_un;
dynet::Parameter all_span_representation::_p_o_span0_init_un;
dynet::Parameter all_span_representation::_p_b_span0_end;
dynet::Parameter all_span_representation::_p_o_span0_end;
dynet::Parameter all_span_representation::_p_b_span0_end_un;
dynet::Parameter all_span_representation::_p_o_span0_end_un;
dynet::Parameter all_span_representation::_p_b_span0_split;
dynet::Parameter all_span_representation::_p_o_span0_split;


//dynet::Parameter all_span_representation::_p_W_span1_init;
dynet::Parameter all_span_representation::_p_W_span1_init_n;
dynet::Parameter all_span_representation::_p_W_span1_init_w;
dynet::Parameter all_span_representation::_p_b_span1_init;
dynet::Parameter all_span_representation::_p_o_span1_init;
dynet::Parameter all_span_representation::_p_b_span1_end;
dynet::Parameter all_span_representation::_p_o_span1_end;
dynet::Parameter all_span_representation::_p_b_span1_split;
dynet::Parameter all_span_representation::_p_o_span1_split;

all_span_representation::all_span_representation(bool init_global,
                                                 dynet::Model& m,
                                                 unsigned span_leve,
                                                 unsigned input_size,
                                                 unsigned nt_embedding_size,
                                                 unsigned hidden_size,
                                                 bool span_midpoints,
                                                 lexicon_representation * l,
                                                 cfg_rule_representation * c)
: span_representation(), computation_attachment(),
  lr(l), cfg(c),
  span_level(span_leve), use_span_midpoints(span_midpoints)
{
  if (init_global)
  {
    //_p_W_span_init = m.add_parameters({hidden_size, input_size + nt_embedding_size});
    _p_W_span_init_n = m.add_parameters({hidden_size, nt_embedding_size});
    _p_W_span_init_w = m.add_parameters({hidden_size, input_size});
    _p_b_span_init = m.add_parameters({hidden_size});
    _p_o_span_init = m.add_parameters({1,hidden_size});
    _p_b_span_init_un = m.add_parameters({hidden_size});
    _p_o_span_init_un = m.add_parameters({1,hidden_size});

    // _p_W_span_end =  m.add_parameters({hidden_size, span_input_dim_base + nt_embedding_size});
    _p_b_span_end = m.add_parameters({hidden_size});
    _p_o_span_end = m.add_parameters({1,hidden_size});

    _p_b_span_end_un = m.add_parameters({hidden_size});
    _p_o_span_end_un = m.add_parameters({1,hidden_size});

    // _p_W_span_split =  m.add_parameters({hidden_size, span_input_dim_base + nt_embedding_size});
    _p_b_span_split = m.add_parameters({hidden_size});
    _p_o_span_split = m.add_parameters({1,hidden_size});


    //_p_W_span0_init = m.add_parameters({hidden_size, input_size + nt_embedding_size});
    _p_W_span0_init_n = m.add_parameters({hidden_size, nt_embedding_size});
    _p_W_span0_init_w = m.add_parameters({hidden_size, input_size});
    _p_b_span0_init = m.add_parameters({hidden_size});
    _p_o_span0_init = m.add_parameters({1,hidden_size});
    _p_b_span0_end = m.add_parameters({hidden_size});
    _p_o_span0_end = m.add_parameters({1,hidden_size});
    _p_b_span0_init_un = m.add_parameters({hidden_size});
    _p_o_span0_init_un = m.add_parameters({1,hidden_size});
    _p_b_span0_end_un = m.add_parameters({hidden_size});
    _p_o_span0_end_un = m.add_parameters({1,hidden_size});
    _p_b_span0_split = m.add_parameters({hidden_size});
    _p_o_span0_split = m.add_parameters({1,hidden_size});

    //_p_W_span1_init = m.add_parameters({hidden_size, input_size + nt_embedding_size});
    _p_W_span1_init_n = m.add_parameters({hidden_size, nt_embedding_size});
    _p_W_span1_init_w = m.add_parameters({hidden_size, input_size});

    _p_b_span1_init = m.add_parameters({hidden_size});
    _p_o_span1_init = m.add_parameters({1,hidden_size});
    _p_b_span1_end = m.add_parameters({hidden_size});
    _p_o_span1_end = m.add_parameters({1,hidden_size});
    _p_b_span1_split = m.add_parameters({hidden_size});
    _p_o_span1_split = m.add_parameters({1,hidden_size});

    _p_W_span_left = m.add_parameters({hidden_size, input_size });
    _p_W_span_right = m.add_parameters({hidden_size, input_size });
    //if (use_span_midpoints) // TODO templatize code to be able to remove this
    _p_W_span_mid = m.add_parameters({hidden_size, input_size });
    _p_W_span_distance = m.add_parameters({hidden_size,1});

    switch (span_level) {
      case 1: { // only the span with distance
        break;
      }
      case 2: { // the span + distance + natural/artificial nt
        _p_W_span_extra = m.add_parameters({hidden_size,1});
        break;
      }
      default: { // the span + distance + nt embeddding
        _p_W_span_extra = m.add_parameters({hidden_size,nt_embedding_size});
        break;
      }
    }
    _p_b_span_bin = m.add_parameters({hidden_size});
    _p_o_span_bin = m.add_parameters({1,hidden_size});

    _p_b_span_un = m.add_parameters({hidden_size});
    _p_o_span_un = m.add_parameters({1,hidden_size});
  }
}


void all_span_representation::clear()
{
  span_scores_bin.clear();
  span_scores_un.clear();

  span_expressions_bin.clear();
  span_expressions_un.clear();
}

void all_span_representation::precompute_span_expressions(const std::vector<int>& lhs_int,
                                                          const std::vector<int>& rhs0_int,
                                                          const std::vector<int>& rhs1_int,
                                                          const std::vector<Word>& words,
                                                          bool train_mode)
{
  std::lock_guard<std::mutex> guard(computation_attachment::cg_mutex);

  // span bigrams
  auto&& Wl = dynet::expr::parameter(*cg, _p_W_span_left);
  auto&& Wr = dynet::expr::parameter(*cg, _p_W_span_right);
  auto&& Wm = dynet::expr::parameter(*cg, _p_W_span_mid);

  auto&& Wd = dynet::expr::parameter(*cg, _p_W_span_distance);
  auto&& bbin = dynet::expr::parameter(*cg, _p_b_span_bin);
  auto&& obin = dynet::expr::parameter(*cg, _p_o_span_bin);
  auto&& bun = dynet::expr::parameter(*cg, _p_b_span_un);
  auto&& oun = dynet::expr::parameter(*cg, _p_o_span_un);


  // lhs + span init
  //auto&& Wi = dynet::expr::parameter(*cg,_p_W_span_init);
  auto&& Win = dynet::expr::parameter(*cg,_p_W_span_init_n);
  auto&& Wiw = dynet::expr::parameter(*cg,_p_W_span_init_w);
  auto&& bi = dynet::expr::parameter(*cg,_p_b_span_init);
  auto&& oi = dynet::expr::parameter(*cg,_p_o_span_init);
  auto&& biu = dynet::expr::parameter(*cg,_p_b_span_init_un);
  auto&& oiu = dynet::expr::parameter(*cg,_p_o_span_init_un);


  // lhs + span end
  //auto&& We = dynet::expr::parameter(*cg,_p_W_span_end);
  auto&& be = dynet::expr::parameter(*cg,_p_b_span_end);
  auto&& oe = dynet::expr::parameter(*cg,_p_o_span_end);
  auto&& beu = dynet::expr::parameter(*cg,_p_b_span_end_un);
  auto&& oeu = dynet::expr::parameter(*cg,_p_o_span_end_un);

  //lhs + span split
  //auto&& Ws = dynet::expr::parameter(*cg,_p_W_span_split);
  auto&& bs = dynet::expr::parameter(*cg,_p_b_span_split);
  auto&& os = dynet::expr::parameter(*cg,_p_o_span_split);


  //   // rhs0 + span init
  //auto&& W0i = dynet::expr::parameter(*cg,_p_W_span0_init);
  auto&& W0in = dynet::expr::parameter(*cg,_p_W_span0_init_n);
  auto&& W0iw = dynet::expr::parameter(*cg,_p_W_span0_init_w);
  auto&& b0i = dynet::expr::parameter(*cg,_p_b_span0_init);
  auto&& o0i = dynet::expr::parameter(*cg,_p_o_span0_init);
  auto&& b0iu = dynet::expr::parameter(*cg,_p_b_span0_init_un);
  auto&& o0iu = dynet::expr::parameter(*cg,_p_o_span0_init_un);

  // rhs0 + span end
  auto&& b0e = dynet::expr::parameter(*cg,_p_b_span0_end);
  auto&& o0e = dynet::expr::parameter(*cg,_p_o_span0_end);
  auto&& b0eu = dynet::expr::parameter(*cg,_p_b_span0_end_un);
  auto&& o0eu = dynet::expr::parameter(*cg,_p_o_span0_end_un);

  //rhs0 + span split
  auto&& b0s = dynet::expr::parameter(*cg,_p_b_span0_split);
  auto&& o0s = dynet::expr::parameter(*cg,_p_o_span0_split);


    //   // rhs1 + span init
  //auto&& W1i = dynet::expr::parameter(*cg,_p_W_span1_init);
  auto&& W1in = dynet::expr::parameter(*cg,_p_W_span1_init_n);
  auto&& W1iw = dynet::expr::parameter(*cg,_p_W_span1_init_w);
  auto&& b1i = dynet::expr::parameter(*cg,_p_b_span1_init);
  auto&& o1i = dynet::expr::parameter(*cg,_p_o_span1_init);

  // rhs0 + span end
  auto&& b1e = dynet::expr::parameter(*cg,_p_b_span1_end);
  auto&& o1e = dynet::expr::parameter(*cg,_p_o_span1_end);

  //rhs0 + span split
  auto&& b1s = dynet::expr::parameter(*cg,_p_b_span1_split);
  auto&& o1s = dynet::expr::parameter(*cg,_p_o_span1_split);


  std::vector<dynet::expr::Expression> lefts,rights,mids,distances,extras;
  std::vector<dynet::expr::Expression> lhss, rhs0s, rhs1s;


  auto wl = words.size();
  auto nl = SymbolTable::instance_nt().get_symbol_count();



  span_expressions_init.resize(wl);
  span_scores_init.resize(wl);
  span_expressions_init_un.resize(wl);
  span_scores_init_un.resize(wl);

  span_expressions_end.resize(wl);
  span_scores_end.resize(wl);
  span_expressions_end_un.resize(wl);
  span_scores_end_un.resize(wl);

  span_expressions_split.resize(wl);
  span_scores_split.resize(wl);

  span_expressions_rhs0_init.resize(wl);
  span_scores_rhs0_init.resize(wl);
  span_expressions_rhs0_init_un.resize(wl);
  span_scores_rhs0_init_un.resize(wl);

  span_expressions_rhs0_end.resize(wl);
  span_scores_rhs0_end.resize(wl);
  span_expressions_rhs0_end_un.resize(wl);
  span_scores_rhs0_end_un.resize(wl);

  span_expressions_rhs0_split.resize(wl);
  span_scores_rhs0_split.resize(wl);

  span_expressions_rhs1_init.resize(wl);
  span_scores_rhs1_init.resize(wl);

  span_expressions_rhs1_end.resize(wl);
  span_scores_rhs1_end.resize(wl);

  span_expressions_rhs1_split.resize(wl);
  span_scores_rhs1_split.resize(wl);


  lhss.resize(nl);
  for (auto l : lhs_int)
  {
    lhss[l] = Win * cfg->get_nt_expr(l);
  }

  rhs0s.resize(nl);
  for (auto r : rhs0_int)
  {
    rhs0s[r] = W0in * cfg->get_nt_expr(r);
  }

  rhs1s.resize(nl);
  for (auto r : rhs1_int)
  {
    rhs1s[r] = W1in * cfg->get_nt_expr(r);
  }




  auto&& embeddings = lr->get_embeddings();
  // TODO: distinguish symbols based on un/bin
  for (unsigned i = 0; i < wl; ++i)
  {
    lefts.push_back(Wl * embeddings[i]);
    rights.push_back(Wr * embeddings[i]);
    if (use_span_midpoints) mids.push_back(Wm * embeddings[i]);
    distances.push_back(Wd * dynet::expr::input(*cg, i));


    //if (i != (wl-1)) // no binary rhs0 at last position
    if (true)          // TODO: distinguish binaries/unaries
    {
      span_expressions_rhs0_init[i].resize(nl);
      span_scores_rhs0_init[i].resize(nl);
      span_expressions_rhs0_end[i].resize(nl);
      span_scores_rhs0_end[i].resize(nl);
      span_expressions_rhs0_split[i].resize(nl);
      span_scores_rhs0_split[i].resize(nl);

      span_expressions_rhs0_init_un[i].resize(nl);
      span_scores_rhs0_init_un[i].resize(nl);
      span_expressions_rhs0_end_un[i].resize(nl);
      span_scores_rhs0_end_un[i].resize(nl);


      auto&& h0l = W0iw * embeddings[i];
      for (auto r0 : rhs0_int)
      {
        auto&& h0p = rhs0s[r0] + h0l;

        auto&& e0i = o0i * dynet::expr::rectify(h0p + b0i);
        auto&& e0e = o0e * dynet::expr::rectify(h0p + b0e);
        auto&& e0s = o0s * dynet::expr::rectify(h0p + b0s);

        auto&& e0iu = o0iu * dynet::expr::rectify(h0p + b0iu);
        auto&& e0eu = o0eu * dynet::expr::rectify(h0p + b0eu);

        if (train_mode)
        {
          span_expressions_rhs0_init[i][r0] = e0i;
          span_expressions_rhs0_end[i][r0] = e0e;
          span_expressions_rhs0_split[i][r0] = e0s;

          span_expressions_rhs0_init_un[i][r0] = e0iu;
          span_expressions_rhs0_end_un[i][r0] = e0eu;
        }
        span_scores_rhs0_init[i][r0]  = as_scalar(e0i.value());
        span_scores_rhs0_end[i][r0]   = as_scalar(e0e.value());
        span_scores_rhs0_split[i][r0] = as_scalar(e0s.value());

        span_scores_rhs0_init_un[i][r0]  = as_scalar(e0iu.value());
        span_scores_rhs0_end_un[i][r0]   = as_scalar(e0eu.value());
      }
    }


    //////////////////

    //    if (i != 0) // no rhs1 at first
    if (true)
    {
      span_expressions_rhs1_init[i].resize(nl);
      span_scores_rhs1_init[i].resize(nl);
      span_expressions_rhs1_end[i].resize(nl);
      span_scores_rhs1_end[i].resize(nl);
      span_expressions_rhs1_split[i].resize(nl);
      span_scores_rhs1_split[i].resize(nl);

      auto&& h1l = W1iw * embeddings[i];
      for (auto r1 : rhs1_int)
      {
        auto&& h1p = rhs1s[r1] + h1l;

        auto&& e1i = o1i * dynet::expr::rectify(h1p + b1i);
        auto&& e1e = o1e * dynet::expr::rectify(h1p + b1e);
        auto&& e1s = o1s * dynet::expr::rectify(h1p + b1s);

        if (train_mode)
        {
          span_expressions_rhs1_init[i][r1] = e1i;
          span_expressions_rhs1_end[i][r1]   = e1e;
          span_expressions_rhs1_split[i][r1] = e1s;

        }
        span_scores_rhs1_init[i][r1]  = as_scalar(e1i.value());
        span_scores_rhs1_end[i][r1]   = as_scalar(e1e.value());
        span_scores_rhs1_split[i][r1] = as_scalar(e1s.value());
      }
    }

    /////////////////





    span_expressions_init[i].resize(nl);
    span_scores_init[i].resize(nl);


    span_expressions_end[i].resize(nl);
    span_scores_end[i].resize(nl);

    span_expressions_split[i].resize(nl);
    span_scores_split[i].resize(nl);

    span_expressions_init_un[i].resize(nl);
    span_scores_init_un[i].resize(nl);
    span_expressions_end_un[i].resize(nl);
    span_scores_end_un[i].resize(nl);


    auto&& hl = Wiw * embeddings[i];
    for (auto l : lhs_int)
    {
      auto&& hp = lhss[l] + hl;
      // auto&& ei = oi * dynet::expr::rectify(Wi * inp + bi);
      // auto&& ee = oe * dynet::expr::rectify(We * inp + be);
      // auto&& es = os * dynet::expr::rectify(Ws * inp + bs);

      auto&& ei = oi * dynet::expr::rectify(hp + bi);
      auto&& ee = oe * dynet::expr::rectify(hp + be);
      auto&& es = os * dynet::expr::rectify(hp + bs);

      auto&& eiu = oiu * dynet::expr::rectify(hp + biu);
      auto&& eeu = oeu * dynet::expr::rectify(hp + beu);

      if (train_mode)
      {
        span_expressions_init[i][l] = ei;
        span_expressions_end[i][l] = ee;
        span_expressions_split[i][l] = es;

        span_expressions_init_un[i][l] = eiu;
        span_expressions_end_un[i][l] = eeu;
      }
      span_scores_init[i][l]  = as_scalar(ei.value());
      span_scores_end[i][l]   = as_scalar(ee.value());
      span_scores_split[i][l] = as_scalar(es.value());

      span_scores_init_un[i][l]  = as_scalar(eiu.value());
      span_scores_end_un[i][l]   = as_scalar(eeu.value());
    }
  }

  switch (span_level) {
    case 1: {
      for (unsigned i = 0; i < wl; ++i)
      {
        for (unsigned j = i; j < wl; ++j)
        {
          auto&& e = lefts[i] + rights[j] + distances[j-i];
          auto&& g = oun * dynet::expr::rectify(e + bun);
          auto&& t = std::make_tuple(i,j,0);
          if (train_mode) span_expressions_un[t] = g;
          span_scores_un[t] = as_scalar(g.value());

          if (use_span_midpoints)
            for (unsigned k = i+1; k <= j; ++k)
            {
              auto&& eb = e + mids[k];
              auto&& f = obin * dynet::expr::rectify(eb + bbin);
              auto&& t = std::make_tuple(i,j,k,0);
              if (train_mode) span_expressions_bin[t] = f;
              span_scores_bin[t] = as_scalar(f.value());
            }
          else
          {
            auto&& f = obin * dynet::expr::rectify(e + bbin);
            auto&& t = std::make_tuple(i,j,0,0);
            if (train_mode) span_expressions_bin[t] = f;
            span_scores_bin[t] = as_scalar(f.value());
          }
        }
      }
      break;
    }
    case 2: {
      auto&& We = dynet::expr::parameter(*cg, _p_W_span_extra);
      std::vector<dynet::expr::Expression> extras;
      for (unsigned l = 0; l < 2; ++l)
        extras.push_back(We * dynet::expr::input(*cg, l));

      for (unsigned l = 0; l < 2; ++l)
      {
        auto&& e1 = extras[l];
        for (unsigned i = 0; i < wl; ++i)
        {
          auto&& e2 = e1 + lefts[i];
          for (unsigned j = i; j < wl; ++j)
          {
            auto&& e = e2 + rights[j] + distances[j-i];
            auto&& g = oun * dynet::expr::rectify(e + bun);
            if (l==1) // artificial for unaries makes no sense
            {
              auto&& t = std::make_tuple(i,j,l);
              if (train_mode) span_expressions_un[t] = g;
              span_scores_un[t] = as_scalar(g.value());
            }
            if (use_span_midpoints)
              for (unsigned k = i+1; k <= j; ++k)
              {
                auto&& eb = e + mids[k];
                auto&& f = obin * dynet::expr::rectify(eb + bbin);
                auto&& t = std::make_tuple(i,j,k,l);
                if (train_mode) span_expressions_bin[t] = f;
                span_scores_bin[t] = as_scalar(f.value());
              }
            else
            {
              auto&& f = obin * dynet::expr::rectify(e + bbin);
              auto&& t = std::make_tuple(i,j,0,l);
              if (train_mode) span_expressions_bin[t] = f;
              span_scores_bin[t] = as_scalar(f.value());
            }
          }
        }
      }
      break;
    }
    default:
      auto&& We = dynet::expr::parameter(*cg, _p_W_span_extra);
      std::vector<dynet::expr::Expression> extras;
      for (unsigned l =0;l < lhs_int.size(); ++l)
        extras.push_back(We * cfg->get_nt_expr(lhs_int[l]));

      for (unsigned l = 0;l < lhs_int.size(); ++l)
      {
        auto&& e1 = extras[l];
        for (unsigned i = 0; i < wl; ++i)
        {
          auto&& e2 = e1 + lefts[i];
          for (unsigned j = i; j < wl; ++j)
          {
            auto&& e = e2 + rights[j] + distances[j-i];
            auto&& g = oun * dynet::expr::rectify(e + bun);

            auto&& t = std::make_tuple(i,j,lhs_int[l]);
            if (train_mode) span_expressions_un[t] = g;
            span_scores_un[t] = as_scalar(g.value());

            if (use_span_midpoints)
              for (unsigned k = i+1; k <= j; ++k)
              {
                auto&& eb = e + mids[k];
                auto&& f = obin * dynet::expr::rectify(eb + bbin);

                auto&& t = std::make_tuple(i,j,k,lhs_int[l]);
                if (train_mode) span_expressions_bin[t] = f;
                span_scores_bin[t] = as_scalar(f.value());

              }
            else
            {
                auto&& f = obin * dynet::expr::rectify(e + bbin);

                auto&& t = std::make_tuple(i,j,0,lhs_int[l]);
                if (train_mode) span_expressions_bin[t] = f;
                span_scores_bin[t] = as_scalar(f.value());
            }
          }
        }
      }
      break;
  }
}


double all_span_representation::get_span_score_lhs_begin(int lhs, int begin) {return span_scores_init[begin][lhs];}
double all_span_representation::get_span_score_lhs_end  (int lhs, int end)   {return span_scores_end[end][lhs];}

double all_span_representation::get_span_score_lhs_begin_unary(int lhs, int begin) {return span_scores_init_un[begin][lhs];}
double all_span_representation::get_span_score_lhs_end_unary  (int lhs, int end)   {return span_scores_end_un[end][lhs];}


double all_span_representation::get_span_score_lhs_split(int lhs, int split) {return span_scores_split[split][lhs];}

double all_span_representation::get_span_score_rhs0_begin(int rhs,int begin) {return span_scores_rhs0_init[begin][rhs];}
double all_span_representation::get_span_score_rhs0_end(int rhs,int end) {return span_scores_rhs0_end[end][rhs];}
double all_span_representation::get_span_score_rhs0_begin_unary(int rhs,int begin) {return span_scores_rhs0_init_un[begin][rhs];}
double all_span_representation::get_span_score_rhs0_end_unary(int rhs,int end) {return span_scores_rhs0_end_un[end][rhs];}
double all_span_representation::get_span_score_rhs0_split(int rhs,int split) {return span_scores_rhs0_split[split][rhs];}



double all_span_representation::get_span_score_rhs1_begin(int rhs,int begin) {return span_scores_rhs1_init[begin][rhs];}
double all_span_representation::get_span_score_rhs1_end(int rhs,int end) {return span_scores_rhs1_end[end][rhs];}
double all_span_representation::get_span_score_rhs1_split(int rhs,int split) {return span_scores_rhs1_split[split][rhs];}


double all_span_representation::get_span_score_bin_info(int begin, int end, int split, int root_info)
{
  return span_scores_bin[std::make_tuple(begin,end, split, root_info)];
}

double all_span_representation::get_span_score_una_info(int begin, int end, int root_info)
{
  return span_scores_un[std::make_tuple(begin,end,root_info)];
}

dynet::expr::Expression& all_span_representation::get_span_expr_lhs_init(int lhs, int begin)
{
  return span_expressions_init[begin][lhs];
}

dynet::expr::Expression& all_span_representation::get_span_expr_lhs_init_unary(int lhs, int begin)
{
  return span_expressions_init_un[begin][lhs];
}

dynet::expr::Expression& all_span_representation::get_span_expr_lhs_end (int lhs, int end)
{
  return span_expressions_end[end][lhs];
}

dynet::expr::Expression& all_span_representation::get_span_expr_lhs_end_unary (int lhs, int end)
{
  return span_expressions_end_un[end][lhs];
}

dynet::expr::Expression& all_span_representation::get_span_expr_lhs_split(int lhs, int split)
{
  return span_expressions_split[split][lhs];
}


dynet::expr::Expression& all_span_representation::get_span_expr_rhs0_init(int rhs0, int begin)
{
  return span_expressions_rhs0_init[begin][rhs0];
}


dynet::expr::Expression& all_span_representation::get_span_expr_rhs0_end(int rhs0, int end)
{
  return span_expressions_rhs0_end[end][rhs0];
}


dynet::expr::Expression& all_span_representation::get_span_expr_rhs0_init_unary(int rhs0, int begin)
{
  return span_expressions_rhs0_init_un[begin][rhs0];
}


dynet::expr::Expression& all_span_representation::get_span_expr_rhs0_end_unary(int rhs0, int end)
{
  return span_expressions_rhs0_end_un[end][rhs0];
}


dynet::expr::Expression& all_span_representation::get_span_expr_rhs0_split(int rhs0, int split)
{
  return span_expressions_rhs0_split[split][rhs0];
}


dynet::expr::Expression& all_span_representation::get_span_expr_rhs1_init(int rhs1, int begin)
{
  return span_expressions_rhs1_init[begin][rhs1];
}


dynet::expr::Expression& all_span_representation::get_span_expr_rhs1_end(int rhs1, int end)
{
  return span_expressions_rhs1_end[end][rhs1];
}


dynet::expr::Expression& all_span_representation::get_span_expr_rhs1_split(int rhs1, int split)
{
  return span_expressions_rhs1_split[split][rhs1];
}





dynet::expr::Expression& all_span_representation::get_span_expr_lhs_info(int lhs, int begin, int end, int medium)
{
  int lhs_code = 0;

  switch (span_level)
  {
    case 1:
      {
        lhs_code = 0;
        break;
      }
    case 2:
      {
        lhs_code = is_artificial(lhs) ? 0 : 1;
        break;
      }
    default:
      lhs_code = lhs;
      break;
  }

  if (medium >=0 && not use_span_midpoints) medium = 0;

  return
      medium >= 0 ?
      span_expressions_bin[std::make_tuple(begin,end, medium, lhs_code)]
      :
      span_expressions_un[std::make_tuple(begin,end, lhs_code)];
}
