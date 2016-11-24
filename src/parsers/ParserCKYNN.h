// -*- mode: c++ -*-
#pragma once


// #include "grammars/Grammar.h"
#include "SimpleChartCKY.h"
#include "edges/Edge.h"
#include "PCKYBestCell.h"

#include "ParserCKY.h"

#include "neural/NNScorer.h"


struct PNNTypes
{

  typedef nn_scorer scorer;
  typedef ::Edge Edge;
  typedef Word ChartWord ;
  typedef PCKYBestCell<Edge, scorer> Cell ;
  typedef ChartCKY<PNNTypes> Chart ;
};



class ParserCKYNN :
    public ParserCKY<Grammar<Rule,Rule,Rule>>,
    public PNNTypes
{
  typedef Grammar<Rule,Rule,Rule> MyGrammar;
  typedef PNNTypes MyTypes;

 public:
  ~ParserCKYNN() {};

  ParserCKYNN(MyGrammar& grammar) : ParserCKY(grammar), word_signature(nullptr)
                                    {};


  void set_word_signature(WordSignature* ws) {word_signature = ws;};
  const WordSignature* get_word_signature() const {return word_signature;}


  /**
      \brief parses the sentence using the grammar
      \param chart the chart to fill
  */
  void parse(typename MyTypes::Chart& chart,
             typename MyTypes::scorer& scorer) const;

  /** \brief Add unary rules at this position in the chart
      (only consider non-terminals created from pos)
      \param cell the cell to fill
      \param isroot true if cell is root
  */
  void add_unary_init(typename MyTypes::Cell& cell,
                      bool isroot,
                      typename MyTypes::scorer& scorer) const;


  /** \brief Add unary rules at this position in the chart
      (only consider non-terminals created from binary rules)
      \param cell the cell to fill
      \param isroot true if cell is root
  */
  void add_unary(typename MyTypes::Cell& cell, bool isroot,
                 typename MyTypes::scorer& scorer) const;


    /**
     \brief finds all extensions of an edge using unary rules
     \param cell the cell to fill
     \param edge_ptr the address of the edge to extend
     \param isroot true if cell is root
   */
  void follow_unary_chain(typename MyTypes::Cell& cell,
                          const Edge * edge_ptr, bool isroot,
                          typename MyTypes::scorer& scorer) const;

  /**
     \brief processes the internal rules, adding edges to the chart where appropriate
     \param chart the chart to fill
  */
  void process_internal_rules(typename MyTypes::Chart& chart,
                              typename MyTypes::scorer& scorer) const;

  /**
     \brief fill the result cell with the most probable edges for each lhs,
     created from edges contained in left_cell and right_cell
     \param left_cell the  leftmost cell to combine
     \param right_cell the rightmost cell to combine
     \param result_cell the cell to store new edges
  */
  void get_candidates(const typename MyTypes::Cell& left_cell,
                      const typename MyTypes::Cell& right_cell,
                      typename MyTypes::Cell& result_cell,
                      typename MyTypes::scorer& scorer) const;


  void set_nbthreads(unsigned u) {nbthreads = u;}
  unsigned get_nbthreads() const {return nbthreads;}

 private:

  WordSignature * word_signature;

  unsigned nbthreads;

};
