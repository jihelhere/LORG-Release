#ifndef _PACKEDEDGE_HPP_
#define _PACKEDEDGE_HPP_

#include "utils/ConfigTable.h"

#include "PackedEdge.h"

template<class Types>
inline
void PackedEdge<Types>::add_daughters(Edge & left, const UnaryRule* rule)
{
  if (not this->open) {
    this->open = true;
    this->local_resize_annotations(1);
  }
  unary_daughters.push_back(UnaryDaughter(left,rule));
}

template<class Types>
inline
void PackedEdge<Types>::add_daughters(Edge & left, Edge & right, const BinaryRule* rule)
{
  //               BLOCKTIMING("PackedEdge add_daughters(binary)");
  if (not this->open) {
    this->open = true;
    this->local_resize_annotations(1);
  }
  binary_daughters.push_back(BinaryDaughter(left,right,rule));
}

template<class Types>
inline
void PackedEdge<Types>::add_daughters(const LexicalRule* rule, const Word* w)
{
  if (not this->open) {
    this->open = true;
    this->local_resize_annotations(1);
  }
  lexical_daughters.push_back(LexicalDaughter(rule, w));
}


template<class Types>
inline
AnnotationInfo& PackedEdge<Types>::get_annotations() {return annotations;}

template<class Types>
inline
const AnnotationInfo& PackedEdge<Types>::get_annotations() const {return annotations;}

// template<class Types>
// inline
// std::vector<AnnotationInfo>& PackedEdge<Types>::get_annotations_backup() {return best.get_annotations_backup();}

// template<class Types>
// inline
// const std::vector<AnnotationInfo>& PackedEdge<Types>::get_annotations_backup() const {return best.get_annotations_backup();}

template<class Types>
inline
const typename PackedEdge<Types>::bvector& PackedEdge<Types>::get_binary_daughters() const {return binary_daughters;}


template<class Types>
inline
typename PackedEdge<Types>::bvector& PackedEdge<Types>::get_binary_daughters() {return binary_daughters;}

template<class Types>
inline
const typename PackedEdge<Types>::uvector& PackedEdge<Types>::get_unary_daughters() const {return unary_daughters;}

template<class Types>
inline
typename PackedEdge<Types>::uvector& PackedEdge<Types>::get_unary_daughters() {return unary_daughters;}

template<class Types>
inline
typename PackedEdge<Types>::BinaryDaughter& PackedEdge<Types>::get_binary_daughter(unsigned i) {return binary_daughters[i];}

template<class Types>
inline
const typename PackedEdge<Types>::BinaryDaughter& PackedEdge<Types>::get_binary_daughter(unsigned i) const {return binary_daughters[i];}


template<class Types>
inline
typename PackedEdge<Types>::UnaryDaughter& PackedEdge<Types>::get_unary_daughter(unsigned i) {return unary_daughters[i];}

template<class Types>
inline
const typename PackedEdge<Types>::UnaryDaughter& PackedEdge<Types>::get_unary_daughter(unsigned i) const {return unary_daughters[i];}

template<class Types>
inline
bool PackedEdge<Types>::get_lex() const {return !(lexical_daughters.empty());}

template<class Types>
inline
const typename PackedEdge<Types>::lvector& PackedEdge<Types>::get_lexical_daughters() const {return lexical_daughters;}

template<class Types>
inline
typename PackedEdge<Types>::lvector& PackedEdge<Types>::get_lexical_daughters() {return lexical_daughters;}

template<class Types>
inline
typename PackedEdge<Types>::LexicalDaughter& PackedEdge<Types>::get_lexical_daughter(unsigned i) {return lexical_daughters[i];}

template<class Types>
inline
const typename PackedEdge<Types>::LexicalDaughter& PackedEdge<Types>::get_lexical_daughter(unsigned i) const {return lexical_daughters[i];}

template<class Types>
inline
void PackedEdge<Types>::local_resize_annotations(unsigned size) {annotations.resize(size);}

/////////////////////////////////////////////////////
// parsing
////////////////////////////////////////////////

template<class Types>
inline
const typename Types::EdgeProbability & PackedEdge<Types>::get_prob_model() const {return best;}

template<class Types>
inline
typename Types::EdgeProbability& PackedEdge<Types>::get_prob_model() {return best;}

template<class Types>
inline void PackedEdge<Types>::extend_derivation(unsigned i, bool licence_unaries, const std::vector<double>& log_norms)
{
  best.extend_derivation(this,i, licence_unaries, log_norms);
}



template<class Types>
inline
bool PackedEdge<Types>::valid_prob_at(unsigned i) const
{
  return get_annotations().valid_prob_at(i, LorgConstants::NullProba);
}


//////////

template <class Types>
PathMatrix PackedEdge<Types>::unary_chains = PathMatrix();








// should be renamed (Why Viterbi?) and moved
template <class Types>
void PackedEdge<Types>::set_unary_chains(const PathMatrix& pathmatrix)
{
  unary_chains = pathmatrix;
}

template <class Types>
const PathMatrix& PackedEdge<Types>::get_unary_chains()
{
  return unary_chains;
}


struct c2f_replace_struct_helper
{
  unsigned idx;
  c2f_replace_struct_helper(unsigned i) : idx(i) {};

  //reset the rule to the 'finer_id'th rule in the vector
  //TODO implement the Charniak-style coarse-to-fine method - this will involve replacing this rule with several
  template <typename T>
  void operator ()(T& daughter) const
  {
    daughter.set_rule(static_cast<const typename T::Rule *>(daughter.get_rule()->get(idx)));
  }
};


template <class Types>
void PackedEdge<Types>::replace_rule_probabilities(unsigned i)
{
  // for all possible daughters
  c2f_replace_struct_helper c2f_replacer(i);
  std::for_each(this->binary_daughters.begin(),this->binary_daughters.end(), c2f_replacer);
  std::for_each(this->unary_daughters.begin(),this->unary_daughters.end(),   c2f_replacer);
  std::for_each(this->lexical_daughters.begin(), this->lexical_daughters.end(),   c2f_replacer);
}




template <class Types>
void PackedEdge<Types>::prepare_inside_probability()
{
  this->get_annotations().inside_probabilities_unary_temp.array = this->get_annotations().inside_probabilities.array ;
  for(auto & prob: this->get_annotations().inside_probabilities_unary_temp.array) {
    if (prob != LorgConstants::NullProba) prob = 0;
  }
}

template <class Types>
void PackedEdge<Types>::adjust_inside_probability()
{
  auto& annot = this->get_annotations();

  for (unsigned i = 0; i < annot.inside_probabilities.array.size(); ++i)
    {
      if(annot.inside_probabilities.array[i] != LorgConstants::NullProba)
        annot.inside_probabilities.array[i] += annot.inside_probabilities_unary_temp.array[i];

      if(annot.inside_probabilities_unary_temp.array[i] == LorgConstants::NullProba)
        annot.inside_probabilities.array[i] = LorgConstants::NullProba;
    }
}


template <class Types>
void PackedEdge<Types>::prepare_outside_probability()
{
  this->get_annotations().outside_probabilities_unary_temp.array = this->get_annotations().outside_probabilities.array ;
  for(auto & prob: this->get_annotations().outside_probabilities_unary_temp.array) {
    if (prob != LorgConstants::NullProba) prob = 0;
  }
}

template <class Types>
void PackedEdge<Types>::adjust_outside_probability()
{
  for (unsigned i = 0; i < this->get_annotations().outside_probabilities.array.size(); ++i)
    {
      if(this->get_annotations().outside_probabilities.array[i] != LorgConstants::NullProba)
        this->get_annotations().outside_probabilities.array[i] += this->get_annotations().outside_probabilities_unary_temp.array[i];
    }
}

/////////////
/// clean
////////////

template <class Types>
void PackedEdge<Types>::clean_invalidated_binaries()
{
//   BLOCKTIMING("PackedEdge<Types>::clean_invalidated_binaries()");

//   auto it = binary_daughters.begin();
//   auto removed_begin = binary_daughters.end();
//   while (it != removed_begin) {
//     if (it->points_towards_invalid_edges())
//       *it = *(--removed_begin);
//     else
//       ++it;
//   }

  auto removed_begin = std::remove_if(binary_daughters.begin(),
                         binary_daughters.end(),
                         toFunc(& BinaryDaughter::points_towards_invalid_edges));

  binary_daughters . erase(removed_begin, binary_daughters.end());

  // Reclaim memory !
  binary_daughters . shrink_to_fit() ;
  //   decltype(binary_daughters)(binary_daughters).swap(binary_daughters);

  assert(binary_daughters.capacity() == binary_daughters.size());
}




///////////////////////////////
////// printing
//////////////////////////////

#include <sstream>


// decode the path-matrix to read a unary chain between start and stop
// TODO: create a class for the matrix and encapsulate
// (because matrix for max-rule has a different type)
inline
unsigned decode_path(PtbPsTree& tree,
                     PtbPsTree::depth_first_iterator& pos,
                     const PathMatrix& paths,
                     const asymb& start,
                     const asymb& end)
{
  bool final = false;
  asymb candidate = start;
  unsigned path_length = 0;


  while(!final) {

    std::unordered_map<asymb, std::unordered_map< asymb, asymb> >::const_iterator it1 =
      paths.find(candidate);

    if(it1 == paths.end())
      final = true;

    else {
      std::unordered_map< asymb, asymb>::const_iterator it2 = it1->second.find(end);

      if(it2 == it1->second.end())
        final = true;
      else {
        //  std::cout << "adding node " << SymbolTable::instance_nt()->translate(it2->second.first) << std::endl;

        std::ostringstream node_content;
        node_content << SymbolTable::instance_nt().translate(it2->second.first);
        //        if(append_annot) node_content << "_" << it2->second.second;

        pos=tree.add_last_daughter(pos,node_content.str());
        candidate = it2->second;
        ++path_length;
      }
    }
  }

  return path_length;

}


template <class Types>
void PackedEdge<Types>::to_set(SET<const PackedEdge<Types>*>& results ) const
{
  results.insert(this);

  if(best.get(0).dtrs->is_binary())
  {
    const BinaryDaughter * daughters =  static_cast<const BinaryDaughter*>(best.get(0).dtrs);
      daughters->left_daughter().to_set(results);
      daughters->right_daughter().to_set(results);
  }
  else
    if(best.get(0).dtrs->is_lexical())
    {
      // const LexicalDaughter * daughters =  static_cast<const LexicalDaughter*>(best.get(0).dtrs);
      // results.insert(daughters->left_daughter());
    }
    else
    {
      const UnaryDaughter * daughters =  static_cast<const UnaryDaughter*>(best.get(0).dtrs);
      daughters->left_daughter().to_set(results);
      //std::cout << *(daughters->get_rule()) << std::endl;
    }

}






template <class Types>
PtbPsTree * PackedEdge<Types>::to_ptbpstree(int lhs, unsigned ith_deriv) const
{
  PtbPsTree * tree = NULL;

  std::ostringstream node_content;
  node_content << SymbolTable::instance_nt().translate(lhs) ;

  //  std::cout << SymbolTable::instance_nt()->translate(get_lhs()) << std::endl;
  //  std::cout << "size daughters: " << best_dtrs_vector.size() << std::endl;
  //  if(best_dtrs_vector[0] == NULL) std::cout << "daughter is NULL" << std::endl;


  //  std::cout << "log prob deriv: " << best.get(ith_deriv).probability << std::endl;


  // TODO raise exception
  // This is only to prevent wtf cases
  if(best.get(ith_deriv).probability == - std::numeric_limits<double>::infinity()) {
    //std::cout << "-inf" << std::endl;
    return NULL;
  }

  if(best.get(ith_deriv).dtrs == NULL) {
    //std::cout << "no daughters" << std::endl;
    return NULL;
  }
  if(std::isnan(best.get(ith_deriv).probability)) {
    //std::cout << "invalid prob" << std::endl ;
    return NULL;
  }
  // if(my_isinvalid(best.get(ith_deriv).probability)) {
  //   //g    std::cerr << "invalid prob" << std::endl;
  //   return NULL;
  // }


  //  std::cout << best.get(ith_deriv).probability << std::endl;


    tree = new PtbPsTree(node_content.str());
    PtbPsTree::depth_first_iterator pos = tree->dfbegin();

    //    std::cout << "allocation done" << std::endl;

    if(best.get(ith_deriv).dtrs->is_binary()) {

      //       std::cout
      //           << "b rule: "
      //     << *(static_cast<const
      //     BinaryDaughter*>(best.get(ith_deriv).dtrs)->get_rule())
      //           << std::endl;

      // std::cout << "deriv " << ith_deriv
      //           << "\t" << best.get(ith_deriv).get_left_index()
      //           << "\t" << best.get(ith_deriv).get_right_index() << std::endl;

      const BinaryDaughter * daughters =  static_cast<const BinaryDaughter*>(best.get(ith_deriv).dtrs);
      int rhs0 = daughters->get_rule()->get_rhs0();
      daughters->left_daughter().to_ptbpstree(*tree, pos, rhs0, best.get(ith_deriv).get_left_index());
      int rhs1 = daughters->get_rule()->get_rhs1();
      daughters->right_daughter().to_ptbpstree(*tree, pos, rhs1, best.get(ith_deriv).get_right_index());
    }
    else {

       // std::cout
       //     << "u rule: "
       //     << *(static_cast<const UnaryDaughter*>(best.get(ith_deriv).dtrs)->get_rule())
       //     << std::endl;


       // std::cout << "deriv " << ith_deriv << std::endl;
       // std::cout << "\t" << best.get(ith_deriv).get_left_index() << std::endl;

      const UnaryDaughter * daughters =  static_cast<const UnaryDaughter*>(best.get(ith_deriv).dtrs);
      decode_path(*tree,pos,
                  PackedEdge<Types>::get_unary_chains(),
                  std::make_pair(lhs,0),
                  std::make_pair(daughters->get_rule()->get_rhs0(), best.get(ith_deriv).get_left_index()));
      int rhs0 = daughters->get_rule()->get_rhs0();
      daughters->left_daughter().to_ptbpstree(*tree, pos, rhs0, best.get(ith_deriv).get_left_index());
    }

    return tree;
}

template <class Types>
void PackedEdge<Types>::to_ptbpstree(PtbPsTree& tree,
                                   PtbPsTree::depth_first_iterator& pos, int lhs, unsigned index) const
{
  std::ostringstream node_content;

  // default height for unary chains
  unsigned added_height = 1;

  // either a valid internal node or a lexical one
  assert(best.get(index).dtrs || get_lex());

  // if we can't find a branching (lexical node)
  // the nwe exit
  if(!best.get(index).dtrs)
  {
    //std::cout << "dtrs is null" << std::endl;
    return;
  }


  // if the *branching* leads to a lexical node
  // the we write the preterminal *and* the word
  if(best.get(index).dtrs->is_lexical()) {

    // std::cout << "lex: " << std::endl;
    // std::cout << "lex: " << get_lhs() << std::endl;
    // std::cout << SymbolTable::instance_nt()->translate(get_lhs()) << std::endl;
    // pos = tree.add_last_daughter(pos, SymbolTable::instance_word()->translate(get_lhs()));

    const LexicalDaughter * daughters =  static_cast<const LexicalDaughter*>(best.get(index).dtrs);

     // std::cout
     //     << "l rule: "
    //     << *(static_cast<const
    //     LexicalDaughter*>(best.get(index).dtrs)->get_rule())
         // << std::endl;


    node_content << SymbolTable::instance_nt().translate(lhs);
    pos = tree.add_last_daughter(pos, node_content.str());

    const Word& w = *(daughters->get_word());

    static bool output_forms = ConfigTable::access().exists("always-output-forms");

    std::string s = output_forms ? w.get_form() :
                    (w.get_id() != -1 ? SymbolTable::instance_word().get_label_string(w.get_id()) : LorgConstants::token_unknown);

    pos = tree.add_last_daughter(pos, s);

    // we added two nodes here
    added_height = 2;

  } // if(best.get(index).dtrs->is_lexical())
  else {
    // std::cout << "int: " << lhs << std::endl;
    // std::cout << SymbolTable::instance_nt().translate(lhs) << std::endl;

    node_content << SymbolTable::instance_nt().translate(lhs);

    pos = tree.add_last_daughter(pos, node_content.str());

    //assert(best);
    assert(best.get(index).dtrs);

    if(best.get(index).dtrs->is_binary()) {

       // std::cout
       //     << "b rule: "
      //     << *(static_cast<const
      //     BinaryDaughter*>(best.get(index).dtrs)->get_rule())
           // << std::endl;


      // // std::cout << "binary" << std::endl;
      // std::cout << "deriv " << index
      //           << "\t" << best.get(index).get_left_index()
      //           << "\t" << best.get(index).get_right_index() << std::endl;

      const BinaryDaughter * daughters =  static_cast<const BinaryDaughter*>(best.get(index).dtrs);

      int rhs0 = daughters->get_rule()->get_rhs0();
      daughters->left_daughter().to_ptbpstree(tree, pos, rhs0, best.get(index).get_left_index());
      int rhs1 = daughters->get_rule()->get_rhs1();
      daughters->right_daughter().to_ptbpstree(tree, pos, rhs1, best.get(index).get_right_index());
    } // if(best.get(index).dtrs->is_binary())
    else { //unary branching

       // std::cout
       //     << "u rule: "
       //     << *(static_cast<const UnaryDaughter*>(best.get(index).dtrs)->get_rule())
       //     << std::endl;

      // std::cout << "unary" << std::endl;
      // std::cout << "deriv " << index
      //           << "\t" << best.get(index).get_left_index()
      //           << "\t" << best.get(index).get_right_index() << std::endl;

      const UnaryDaughter * daughters =  static_cast<const UnaryDaughter*>(best.get(index).dtrs);

      std::pair<int,int> fro = std::make_pair(lhs, index);
      std::pair<int,int> to =  std::make_pair(daughters->get_rule()->get_rhs0(),
                                              best.get(index).get_left_index());

      added_height += decode_path(tree,pos,
                                  PackedEdge<Types>::get_unary_chains(),
                                  fro,
                                  to);

      // std::cout << index << "\t" << best_left_indices.size() << std::endl;
      // std::cout << best_left_indices[index] << std::endl;

      int rhs0 = daughters->get_rule()->get_rhs0();
      //      std::cout << *daughters->get_rule() << std::endl;
      daughters->left_daughter().to_ptbpstree(tree, pos, rhs0, best.get(index).get_left_index());
    }

  }

  // go up because you processed the last daughter
  while(added_height--)
    pos.up();
}


// TODO refile this method above


template<class OPEP>
std::ostream& operator<<(std::ostream& out, const PackedEdge<OPEP>& edge)
{
  out << "(edge: " << &edge << ": inside(0): " << edge.get_annotations().get_inside(0) << " best:" << edge.best;
  edge.apply(std::function<void(const typename OPEP::LexicalDaughter&)>([&out](const typename OPEP::LexicalDaughter & dtr){out<<"(dtr: " << *dtr.get_rule() << ") ";}));
  edge.apply(std::function<void(const typename OPEP::BinaryDaughter&)>([&out](const typename OPEP::BinaryDaughter & dtr){out<<"(dtr: " << *dtr.get_rule() << ") ";}));
  edge.apply(std::function<void(const typename OPEP::UnaryDaughter&)>([&out](const typename OPEP::UnaryDaughter & dtr){out<<"(dtr: " << *dtr.get_rule() << ") ";}));
  return out << ") ";
}


template <class Types>
inline
bool PackedEdge<Types>::has_solution(unsigned i) const
{
  return best.has_solution(i);
}


template <typename Types>
double PackedEdge<Types>::marginalise() const
{
  // for(auto & d: this->binary_daughters) {
  //   std::cout << *(d.get_rule()) << std::endl;
  // }
  // for(auto & d: this->unary_daughters) {
  //   std::cout << *(d.get_rule()) << std::endl;
  // }
  // for(auto & d: this->lexical_daughters) {
  //   std::cout << *(d.get_rule()) << std::endl;
  // }


  double normalisation_factor = 0;
  const AnnotationInfo & a = this->get_annotations();
  for (unsigned i = 0; i < a.inside_probabilities.array.size(); ++i)
  {
    if(a.inside_probabilities.array[i] == LorgConstants::NullProba
       || a.outside_probabilities.array[i] == LorgConstants::NullProba
       )
      continue;

    //    std::cout << a.inside_probabilities.array[i] << " " << a.outside_probabilities.array[i] << std::endl;
    normalisation_factor += a.inside_probabilities.array[i] * a.outside_probabilities.array[i];
  }

  //  std::cout << std::log(normalisation_factor) << std::endl;
  return std::log(normalisation_factor);

}




template <typename Types>
void PackedEdge<Types>::update_relaxations(const double & u)
{
  for(auto& le : lexical_daughters)
  {
    le.update_relaxation(-u);
  }
  for(auto& un : unary_daughters)
  {
    un.update_relaxation(-u);
  }
  for(auto& bi : binary_daughters)
  {
    bi.update_relaxation(-u);
  }
}

#endif
