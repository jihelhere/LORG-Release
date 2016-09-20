// -*- mode: c++ -*-
#pragma once

#include <vector>
#include <string>

#ifdef USE_THREADS
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_scheduler_init.h>
#endif

#include <functional>
#include <iostream>

#include "./Bracketing.h"
#include "utils/PtbPsTree.h"
#include "utils/LorgConstants.h"
#include "utils/SymbolTable.h"


/**
  \class ChartCKY
  \brief represents a chart of cells
*/
template<class Types>
class ChartCKY
{
 public:
  typedef typename Types::Cell Cell;
  typedef typename Types::Edge Edge;
  typedef typename Types::ChartWord MyWord;

 private:
  Cell * the_cells; ///< the chart itself
  Edge * the_edges; ///< the edges of the chart
  unsigned size;     ///< the size of the chart (width)
  unsigned nb_cells; ///< number of cells in the chart
  const std::vector< MyWord >& sentence;
  const std::vector<bracketing>& brackets;

  // prevents unwanted conversions
  ChartCKY(const ChartCKY&);
  ChartCKY& operator=(const ChartCKY&);

 public:
#ifdef USE_THREADS
  static unsigned nbthreads;
#endif

  ~ChartCKY();

  /**
     \brief constructor with initialisation
     \param sentence the sentence to create the chart
     \param grammar_size  the number non-terminals in the grammar
     \param brackets chunks
  */
  ChartCKY(const std::vector< MyWord >& sentence, unsigned grammar_size, const std::vector<bracketing>& brackets);

  /**
     \brief get the size of the chart
     \return the size of the chart
  */
  unsigned get_size() const;


  /**
     \brief access a cell of the chart by its coordinates
     \param start starting point of the cell's span
     \param end   ending point of the cell's span
     \return a cell (may segfault if coordinates are out of bounds)
  */

  const Cell& access(unsigned start, unsigned end) const;
  Cell& access(unsigned start, unsigned end);

  inline const Cell& get_root() const;
  Cell& get_root();

  PtbPsTree* get_best_tree(int start_symbol, unsigned k) const;

  SET<const Edge*> get_rules_best_solution(int start_symbol) const;

  double get_score(int start_symbol, unsigned k) const;

  void init(const std::vector< MyWord >& sentence);

  void reset_probabilities();

  bool has_solution(int symb, unsigned i) const;

  void clear();

  void prepare_retry();

  bool is_valid(int start_symbol) const;

  std::ostream & to_stream(std::ostream & s) const;
  std::string toString() const;

  void opencells_apply(const std::function<void(Cell &)>& f);
  void opencells_apply_nothread(const std::function<void(Cell &)>& f);
  void opencells_apply_bottom_up(const std::function<void(Cell &)>& f, unsigned min_span = 0);
  void opencells_apply_bottom_up_nothread(const std::function<void(Cell &)>& f, unsigned min_span = 0);
  void opencells_apply_top_down(const std::function<void(Cell &)>& f);
  void opencells_apply_top_down_nothread(const std::function<void(Cell &)>& f);
  void opencells_apply_top_down_nothread(const std::function<void(const Cell &)>& f) const;

  std::ostream & operator<<(std::ostream & out) { opencells_apply_bottom_up([&out](Cell & cell){return out << cell << std::endl; }); return out; }


  void update_relaxations(bool, const MAP<int, MAP<int, MAP<int, double>>>&,
                          const std::unordered_map<int, int>&);
};
