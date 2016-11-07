// -*- mode: c++ -*-
#pragma once

#include <cassert>
#include <cstring>
#include "edges/Edge.h"
#include "Word.h"
#include "rules/Rule.h"

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winconsistent-missing-override"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weffc++"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-local-typedef"
#else
//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#if defined(__clang__)
#pragma clang diagnostic pop
#pragma clang diagnostic pop
#pragma clang diagnostic pop
#pragma clang diagnostic pop
#else
//#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#endif


/**
  \class PCKYBestCell
  \brief represents a cell in a chart
  that only accepts new or more probable edges
*/
#define UNARY_LENGTH 2

template<class MyEdge, typename probability_extractor>
class PCKYBestCell {
 private:
  MyEdge ** edges;
  bool closed;
  const MyEdge * word_edge;
  unsigned begin;
  unsigned end;
  bool top;

  static unsigned max_size;

public:
  static const unsigned unary_length = UNARY_LENGTH;

  typedef MyEdge CellEdge;

  /**
     \brief Simple constructor
     \note this cell is not initialised: trying
     to use it will result in segfault !
     You have to call init first
   */
  PCKYBestCell() : edges(nullptr), closed(true) {};

  /**
     \brief Constructor
     \param cl true is closed
   */
  PCKYBestCell(bool cl);

  /**
     \brief destructor
   */
  ~PCKYBestCell();


  /**
     \brief initialise the cell
     \param cl true is closed
   */
  void init(bool cl, unsigned begin, unsigned end, bool top);
  void reinit(bool cl);
  inline bool get_top() const { return top; }
  inline unsigned get_begin() const { return begin; }
  inline unsigned get_end() const { return end; }

  /**
     \brief insert a candidate edge in the cell
     \param e a pointer to the candidate edge
     \return a pointer to the newly created edge or nullptr if no insertion took place
  */
  const MyEdge * process_candidate(int s, const MyEdge &e);


  /**
     \brief test if there's an edge of given lhs in the cell
     \param label a lhs label
     \return true if label is in the cell
  */
  bool exists_edge(int s, int label) const;
  bool exists_edge_all_height(int label) const;


  /**
     \brief add an edge  to the cell
     \param edge the edge to add
     \return a pointer to the newly created edge
  */
  const MyEdge * add_edge(int s, const MyEdge& edge);


  /**
     \brief access an edge by its lhs
     \param i the label of the edge
     \return the edge with i as lhs
  */
  const MyEdge& at(int s, int i ) const;


  /**
     \brief access the mostprobable edge by its lhs
     \param i the label of the edge
     \return the best edge with i as lhs
  */
  const MyEdge& get_edge(int s, int i) const;
  const MyEdge& get_best_edge(int i) const;

  /**
     \brief access
     \return true if the cell is closed
  */
  bool is_closed() const;



  /**
     \brief removes edges to far from the most probable one
   */

  void apply_beam();


  void set_word_edge(const MyEdge * we);
  const MyEdge * get_word_edge() const;

  void add_word(const Word& word, probability_extractor& scorer);


  /**
     \brief Output operator
     \param out the ostream to write on
     \param cell the cell object to write
    \return the used ostream
  */
  template<class E, typename p>
  friend std::ostream& operator<<(std::ostream& out, const PCKYBestCell<E,p>& cell);


  static void set_max_size(unsigned size);

};

template<class MyEdge, typename probability_extractor>
unsigned PCKYBestCell<MyEdge,probability_extractor>::max_size = 0;

template<class MyEdge, typename probability_extractor>
PCKYBestCell<MyEdge,probability_extractor>::~PCKYBestCell()
{
  if(!closed) {
    size_t s = top ? UNARY_LENGTH + 1 : UNARY_LENGTH;
    for(unsigned i = 0; i < max_size * s;++i) delete edges[i];
    delete edges;
    delete word_edge;
  }
}

template<class MyEdge, typename probability_extractor>
inline
bool PCKYBestCell<MyEdge,probability_extractor>::exists_edge(int s, int label) const
{
  return (edges[s * max_size + label] != nullptr);
}

template<class MyEdge, typename probability_extractor>
inline
bool PCKYBestCell<MyEdge,probability_extractor>::exists_edge_all_height(int label) const
{
  size_t s = top ? unary_length +1 : unary_length;
  for (unsigned h = 0; h < s; ++h)
  {
    if (exists_edge(h,label)) return true;
  }
  return false;
}


template<class MyEdge, typename probability_extractor>
inline
const MyEdge * PCKYBestCell<MyEdge,probability_extractor>::process_candidate(int s, const MyEdge& candidate)
{
  MyEdge ** current = &edges[s * max_size + candidate.get_lhs()];

  if(*current)
  {
    if (candidate.get_probability() > (*current)->get_probability())
    {
        (*current)->replace(candidate);
    }
    else
      return nullptr;
  }
  else
  {
    if (s > 0)
    {
      int sp = 0;
      MyEdge ** previous = nullptr;
      while (sp != s)
      {
        previous = &edges[sp * max_size + candidate.get_lhs()];
        if (*previous) break;
        ++sp;
      }

      if(*previous)
      {
        if (candidate.get_probability() > (*previous)->get_probability())
        {
          *current = new MyEdge(candidate);
        }
        else
          return nullptr;
      }
      else
      {
        *current = new MyEdge(candidate);
      }
    }
    else
      *current = new MyEdge(candidate);
  }
  return *current;
}

template<class MyEdge, typename probability_extractor>
inline
const MyEdge * PCKYBestCell<MyEdge,probability_extractor>::add_edge(int s, const MyEdge& edge)
{
  return edges[s * max_size + edge.get_lhs()] = new MyEdge(edge);
}



template<class MyEdge, typename probability_extractor>
inline
const MyEdge& PCKYBestCell<MyEdge,probability_extractor>::at(int s, int i) const
{
  //assert(i>=0);
  //assert( i < (int) max_size);
  assert(i>=0 && i < (int) max_size);

  return *edges[s * max_size + i];
}

template<class MyEdge, typename probability_extractor>
inline
const MyEdge& PCKYBestCell<MyEdge,probability_extractor>::get_edge(int s, int i) const
{
  return at(s, i);
}


// assume the edge is in the cell
template<class MyEdge, typename probability_extractor>
inline
const MyEdge& PCKYBestCell<MyEdge,probability_extractor>::get_best_edge(int i) const
{
  int bh = -1;
  double bs = - std::numeric_limits<double>::infinity();

  size_t s = top ? unary_length + 1 : unary_length;

  for (unsigned h = 0; h < s; ++h)
  {
    const auto edgep = edges[h * max_size + i];
    if (edgep)
    {
      auto p = edgep->get_probability();
      if (p > bs)
      {
        bh = h;
        bs = p;
      }
    }
  }

  return at(bh, i);
}



template<class MyEdge, typename probability_extractor>
inline
void PCKYBestCell<MyEdge,probability_extractor>::init(bool cl, unsigned int b, unsigned int e, bool t)
{
  begin = b; end = e; top = t;
  reinit(cl);
}

template<class MyEdge, typename probability_extractor>
inline
void PCKYBestCell<MyEdge,probability_extractor>::reinit(bool cl)
{
  if(!(closed = cl)) {
    size_t s = top ? unary_length + 1: unary_length;
    word_edge = nullptr;
    edges =  new MyEdge * [s * max_size];
    memset(edges, 0, s * max_size * sizeof(MyEdge*));
    //   for(unsigned i = 0; i < max_size;++i)
    //     edges[i]=nullptr;
  }
}


template<class MyEdge, typename probability_extractor>
inline
bool PCKYBestCell<MyEdge,probability_extractor>::is_closed() const
{ return closed; }


template<class MyEdge, typename probability_extractor>
inline
void PCKYBestCell<MyEdge,probability_extractor>::set_word_edge(const MyEdge * we)
{
  word_edge = we;
}

template<class MyEdge, typename probability_extractor>
inline
const MyEdge * PCKYBestCell<MyEdge,probability_extractor>::get_word_edge() const
{
  return word_edge;
}


template<class MyEdge, typename probability_extractor>
inline
void PCKYBestCell<MyEdge,probability_extractor>::add_word(const Word & word,
                                                          probability_extractor& scorer)
{
  if(!word_edge)
  {
    set_word_edge(new MyEdge(word.get_id(),0,true, &word.get_form()));
  }

  for(const auto& r : word.get_rules())
  {
    //std::cerr << r << std::endl;
    MyEdge e(r->get_lhs(), word_edge, scorer.compute_lexical_score(get_begin(), r));
    //MyEdge e(r->get_lhs(), word_edge, static_cast<const Rule*>(r)->get_probability());
    e.set_pruning_probability(static_cast<const Rule*>(r)->get_probability());
    (void) add_edge(0, e);
  }
}

template<class MyEdge, typename probability_extractor>
inline
void PCKYBestCell<MyEdge,probability_extractor>::apply_beam()
{
  // TODO: remove the define
#define BEAM_CONST 60

  double max_prob = - std::numeric_limits<double>::infinity();

  size_t s = top ? unary_length + 1 : unary_length;


  for(unsigned i = 0 ; i < max_size * s; ++i)
    if(edges[i] && max_prob < edges[i]->get_probability())
      max_prob = edges[i]->get_probability();


  for(unsigned i = 0 ; i < max_size * s; ++i)
    if(edges[i] && max_prob > edges[i]->get_probability() + BEAM_CONST) {
      delete edges[i];
      edges[i] =  nullptr;
    }
}

template<class MyEdge, typename probability_extractor>
inline
void PCKYBestCell<MyEdge,probability_extractor>::set_max_size(unsigned size)
{
  max_size = size;
}
