// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVAR_H_
#define _PARSERCKYALLMAXVAR_H_

#include "parsers/ParserCKYAll.h"
#include "parsers/ParserCKYAllFactory.h"
//#include "utils/lorg_functional.h"

template<class Types>
class ParserCKYAllMaxRule : public ParserCKYAll_Impl<Types>
{
 public:
  typedef typename ParserCKYAll_Impl<Types>::AGrammar AGrammar;
  typedef typename Types::Chart Chart;
  typedef typename Types::Cell Cell;
  typedef typename Types::Edge Edge;
  typedef typename Types::EdgeProbability ProbaModel;

    ParserCKYAllMaxRule(ParserCKYAllFactory::MaxParsing_Calculation c,
                        const std::vector<AGrammar*>& cgs,
                        const std::vector<double>& priors,
                        double beam_threshold,
                        const annot_descendants_type& annot_descendants_,
                        bool accurate_, unsigned min_beam_length, int stubborn)
    : ParserCKYAll_Impl<Types>(cgs,
                               priors,
                               beam_threshold,
                               annot_descendants_,
                               accurate_, min_beam_length, stubborn)
    {
     ProbaModel::set_calculation(c);
    }

   /**
   * @brief calculate the chart specific rule probabilities for all packed edges in all cells
   *
   * @return void
   **/
   void calculate_maxrule_probabilities(double log_norms)
   {
     this->chart->opencells_apply_bottom_up(
      [&log_norms](Cell & cell)
      {
        // std::cout << "filling cell " << &cell
        //           << " : ======================================================"
        //           << cell << std::endl;


        for(unsigned i=0; i<cell.get_max_size(); ++i) {/*std::cout << "edge " << i << std::endl ;*/
          if(cell.exists_edge(i))
          {
            auto& e =  cell.get_edge(i);

            auto& p = e.get_prob_model();
            p.init();
            for(auto& d: e.get_lexical_daughters())
              p.update_lexical(e, d, log_norms);
            for(auto& d: e.get_binary_daughters())
              p.update_binary(e, d, log_norms);
          }
        }

        for(unsigned i=0; i<cell.get_max_size(); ++i) {/*std::cout << "edge " << i << std::endl ;*/
          if(cell.exists_edge(i))
          {
            auto& e =  cell.get_edge(i);

            auto& p = e.get_prob_model();
            for(auto& d: e.get_unary_daughters())
              p.update_unary(e, d, log_norms);

            p.finalize();
          }
        }




        // WHY NOT ??
        // cell.apply_on_edges (
        //     toFunc(&ProbaModel::init),
        //     std::make_pair(toFunc(&ProbaModel::update_lexical), log_normalisation_factor),
        //     std::make_pair(toFunc (&ProbaModel::update_binary), log_normalisation_factor));

        // cell.apply_on_edges (
        //     std::make_pair(toFunc(&ProbaModel::update_unary), log_normalisation_factor),
        //     toFunc (&ProbaModel::finalize));

        // std::cout << "best filled for cell " << &cell << " : " << cell << std::endl;
      }
    );
  }
};

#endif /* _PARSERCKYALLMAXVAR_H_ */
