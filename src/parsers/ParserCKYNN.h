// -*- mode: c++ -*-
#pragma once


// #include "grammars/Grammar.h"
#include "SimpleChartCKY.h"
#include "edges/Edge.h"
#include "PCKYBestCell.h"

#include "ParserCKY.h"

#include "NNScorer.h"



struct PNNTypes
{
  typedef nn_scorer scorer;
  typedef ::Edge Edge;
  typedef Word ChartWord ;
  typedef PCKYBestCell<Edge, scorer> Cell ;
  typedef ChartCKY< PNNTypes > Chart ;
};


class ParserCKYNN :
    public ParserCKY<Grammar<Rule,Rule,Rule>>,
    public PNNTypes
{
  typedef Grammar<Rule,Rule,Rule> MyGrammar;

  //typedef compact_binary_rules::vector_brules<const Rule*> vector_brules;

 public:
  //  ParserCKYNN();
  ~ParserCKYNN() {};

  ParserCKYNN(MyGrammar& grammar) : ParserCKY(grammar), word_signature(nullptr)
                                    {};


  void set_word_signature(WordSignature* ws) {word_signature = ws;};
  const WordSignature* get_word_signature() const {return word_signature;}


  /**
      \brief parses the sentence using the grammar
      \param chart the chart to fill
  */
  void parse(Chart& chart, scorer& scorer) const;

  /** \brief Add unary rules at this position in the chart
      (only consider non-terminals created from pos)
      \param cell the cell to fill
      \param isroot true if cell is root
  */
  void add_unary_init(Cell&cell, bool isroot, scorer& scorer) const;


  /** \brief Add unary rules at this position in the chart
      (only consider non-terminals created from binary rules)
      \param cell the cell to fill
      \param isroot true if cell is root
  */
  void add_unary(Cell& cell, bool isroot, scorer& scorer) const;


    /**
     \brief finds all extensions of an edge using unary rules
     \param cell the cell to fill
     \param edge_ptr the address of the edge to extend
     \param isroot true if cell is root
   */
  void follow_unary_chain(Cell& cell, const Edge * edge_ptr, bool isroot, scorer& scorer) const;

  /**
     \brief processes the internal rules, adding edges to the chart where appropriate
     \param chart the chart to fill
  */
  void process_internal_rules(Chart& chart, scorer& scorer) const;

  /**
     \brief fill the result cell with the most probable edges for each lhs,
     created from edges contained in left_cell and right_cell
     \param left_cell the  leftmost cell to combine
     \param right_cell the rightmost cell to combine
     \param result_cell the cell to store new edges
  */
  void get_candidates(const Cell& left_cell, const Cell& right_cell, Cell& result_cell, scorer& scorer) const;


 private:

  WordSignature * word_signature;



};
