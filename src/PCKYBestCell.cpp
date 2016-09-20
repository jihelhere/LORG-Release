#include "PCKYBestCell.h"
#include "utils/SymbolTable.h"
#include<cmath>

template<class MyEdge, typename probability_extractor>
PCKYBestCell<MyEdge,probability_extractor>::PCKYBestCell(bool cl):
  closed(cl), word_edge(NULL)
{
  if(closed)
    edges = NULL;
  else {
    size_t s = max_size * (top ? unary_length +1 : unary_length);
    edges =  new MyEdge * [s];
    memset(edges, 0, s * sizeof(MyEdge*));
    //for(unsigned i = 0; i < max_size;++i) edges[i]=NULL;
  }
}

template<class MyEdge, typename probability_extractor>
std::ostream& operator<<(std::ostream& out, const PCKYBestCell<MyEdge,probability_extractor>& cell)
{
  if(cell.is_closed())
    out << "cell closed\n";
  else {
    size_t s = cell.max_size * (cell.get_top() ?
                                PCKYBestCell<MyEdge,probability_extractor>::unary_length + 1 : PCKYBestCell<MyEdge,probability_extractor>::unary_length);
    for(unsigned i = 0; i < s ; ++i) {
      MyEdge * edge = cell.edges[i];
      if(edge) {
        out << SymbolTable::instance_nt().translate(edge->get_lhs()) << ": "
        << std::exp(edge->get_probability()) << "\n";
      }
    }
  }

  return out;
}
