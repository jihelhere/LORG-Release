#include "PtbPsTree.h"

#include <list>
#include <iostream>
#include <queue>
#include <vector>
#include <cassert>
#include <sstream>

#include "SymbolTable.h"
#include "LorgConstants.h"
#include "RandomGenerator.h"

namespace {

    inline
    bool is_artificial( const std::string& str )
    {
        return str[0] == '['; // hardcoded == Bad -> true
    }

    void replace_number( std::string& s, const boost::regex& num_exp)
    {
        boost::cmatch matched;
        if(boost::regex_match(s.c_str(),matched,num_exp))
            s = LorgConstants::token_number;
    }


void remove_function_from_nt_string( std::string& orig )
  {
    static const boost::regex exp ("^([A-Za-z]+\\$?)[=-].*");
    boost::cmatch matched;

    if(boost::regex_match(orig.c_str(),matched,exp))
      orig = std::string(matched[1].first, matched[1].second);
  }


  void remove_trailing_numbers_from_nt_string( std::string& orig )
  {
    static const boost::regex exp ("^([A-Za-z=-]+\\$?)[=-][0-9]*");
    boost::cmatch matched;

    if(boost::regex_match(orig.c_str(),matched,exp))
      orig = std::string(matched[1].first, matched[1].second);
  }


    // void remove_function_from_nt_string( std::string& orig )
    // {
    //     static const boost::regex exp ("^([A-Za-z]+\\$?)([=-])(.*)");
    //     boost::cmatch matched;

    //     if(boost::regex_match(orig.c_str(),matched,exp))
    //     {
    //       orig = std::string(matched[1].first, matched[1].second);

    //       // -LRB- => need to capture final dash
    //       if( std::string(matched[3].first,matched[3].second) == "")
    //       {
    //         orig += std::string(matched[2].first, matched[2].second);
    //       }

    //     }
    // }


    // void remove_trailing_numbers_from_nt_string( std::string& orig )
    // {
    //     static const boost::regex exp ("^([A-Za-z=-]+\\$?)([=-])([0-9-]*)");
    //     boost::cmatch matched;

    //     if(boost::regex_match(orig.c_str(),matched,exp))
    //     {
    //         orig = std::string(matched[1].first, matched[1].second);

    //       // -LRB- => need to capture final dash
    //       if( std::string(matched[3].first,matched[3].second) == "")
    //       {
    //         orig += std::string(matched[2].first, matched[2].second);
    //       }

    //     }
    // }

}

PtbPsTree::PtbPsTree() : PsTree<Content>() {}

PtbPsTree::PtbPsTree(const Content& content) :PsTree<Content>(content) {}


PtbPsTree::PtbPsTree(const Content& content, const std::vector<PtbPsTree>& daughters)
    : PsTree<Content>(content)
{
  auto current = dfbegin();
  for(auto d_iter = daughters.begin(); d_iter != daughters.end();++d_iter)
    (void) add_last_daughter(current,*d_iter);
}


PtbPsTree::~PtbPsTree(){}

//remove node with unwanted labels
void PtbPsTree::clean(const std::unordered_set<std::string>& labels_to_remove)
{
    for(auto i(dfbegin()); i != dfend(); ++i)
    {
      if(labels_to_remove.count(*i))
      {
        while(!i->has_left_sister() && !i->has_right_sister())
        {
          i.up();
        }
        erase(i);
      }
    }
}

//remove unary chains X -> X
void PtbPsTree::remove_useless_unary_chains()
{
  for(auto i(dfbegin()); i != dfend(); ++i)
  {
    // look for configuration  X -> Y
    if(!i->has_left_sister() && !i->has_right_sister() && i->has_mother() && !i->leaf())
    {
      auto save(i); save.up();
      if(*save != *i) continue;

      auto j(i); j.down_first();
      while(j != dfend())
      {
        add_last_daughter(save,subtree(j));
        j.right();
      }
      erase(i);
    }
  }
}


void PtbPsTree::remove_function()
{
    for(auto i = dfbegin(); i != dfend(); ++i)
      if(!i->leaf()) //words don't have functions
        remove_function_from_nt_string(*i);
}


void PtbPsTree::remove_trailing_numbers()
{
  for(auto i = dfbegin(); i != dfend(); ++i)
    if(!i->leaf()) // don't touch words
      remove_trailing_numbers_from_nt_string(*i);
}


void PtbPsTree::remove_numbers(const boost::regex& num_regex)
{
    for(auto i = lbegin(); i != lend(); ++i)
      replace_number(*i,num_regex);
}




void PtbPsTree::unbinarise()
{
  for(auto i = dfbegin(); i != dfend(); ++i)
  {
    if(is_artificial(*i) && subtree(i).height() > 1)
    {
      auto dter = i;
      for(dter.down_last(); dter != dfend(); dter.left())
        add_right_sister(i, subtree(dter));
      erase(i);
    }
  }
  add_root(subtree(++dfbegin()));
}


namespace {

typedef PtbPsTree::depth_first_iterator iterator;

inline
std::string daughters_to_string(const std::list<iterator>& daughters, const std::string&, int)
{
  std::ostringstream treename;
  treename << '[';
  for(auto iter = daughters.rbegin(); iter != daughters.rend();++iter)
    treename << '(' << *(*iter) <<')';
  treename << ']';

  return treename.str();
}

inline
std::string
n_daughters_left_to_string(const std::list<iterator>& daughters, const std::string& ancestor_name, int markov)
{
  std::ostringstream treename;
  treename << "[|"  + ancestor_name + "|>";
  //  treename << "@"+ ancestor_name;

  auto iter = daughters.rbegin();
  while(markov > 0 && iter != daughters.rend())
  {
    treename << '(' << *(*iter) << ')';
    --markov;
    ++iter;
  }
  treename << ']';

  return treename.str();
}


inline
std::string
n_daughters_right_to_string(const std::list<iterator>& daughters, const std::string& ancestor_name, int markov)
{
  std::ostringstream treename;
  treename << "[|"  + ancestor_name + "|>";
  // treename << "@"+ ancestor_name;

  auto iter = daughters.begin();
  while(markov > 0 && iter != daughters.end())
  {
    treename << '(' << *(*iter) << ')';
    --markov;
    ++iter;
  }
  treename << ']';

  return treename.str();
}


inline
PtbPsTree
create_bintree(const std::string& current_name,
               std::list<iterator>& daughters, Bin_Direction dir,
               const std::string& ancestor_name,
               HorizMarkov mark,
               std::string (*give_name)(const std::list<iterator>&, const std::string&, HorizMarkov))
{
  iterator child;
  if(dir == LEFT)
  {
    child = daughters.back();
    daughters.pop_back();
  }
  else
  {
    child = daughters.front();
    daughters.pop_front();
  }

  PtbPsTree bin(give_name(daughters, ancestor_name,mark));
  iterator bin_iter = bin.dfbegin();
  for(const auto& daughter : daughters)
  {
    bin.add_first_daughter(bin_iter,bin.subtree(daughter));
  }

  PtbPsTree new_tree(current_name);
  iterator new_iter = new_tree.dfbegin();
  if (dir ==  LEFT)
  {
    new_tree.add_first_daughter(new_iter,bin);
    new_tree.add_first_daughter(new_iter,new_tree.subtree(child));
  }
  else
  {
    new_tree.add_first_daughter(new_iter,new_tree.subtree(child));
    new_tree.add_first_daughter(new_iter,bin);
  }
  return new_tree;
}

} // namespace


void PtbPsTree::binarise(Bin_Direction direction, HorizMarkov mark)
{
  assert(direction != NONE );

  std::string (*give_name)(const std::list<iterator>&, const std::string&, HorizMarkov) = NULL;


  if(mark < 0)
    give_name = daughters_to_string;
  else if(direction == LEFT)
    give_name = n_daughters_left_to_string;
  else
    give_name = n_daughters_right_to_string;

  assert(give_name != NULL);

  std::string ancestor_name;

  // visit all nodes, starting from root (depth-first)
  for(iterator pos = dfbegin(); pos != dfend(); ++pos)
  {

    std::list<iterator> daughters;

    //get daughters (in reverse order)
    iterator dter = pos;
    for(dter.down_last(); dter != dfend(); dter.left())
      daughters.push_back(dter);

    if(daughters.size() > 2 ) {

      //should be a const_iterator !
      iterator ancestor_iter = pos;
      while(is_artificial(*ancestor_iter))
        ancestor_iter.up();
      std::string ancestor_name = *ancestor_iter;



      Bin_Direction real_direction;
      if (direction == RAND)
      {
        if (RandomGenerator::instance_binarization()->next() < 0.5)
        {
          real_direction = LEFT;
          give_name = n_daughters_left_to_string;
        }
        else
        {
          real_direction = RIGHT;
          give_name = n_daughters_right_to_string;
        }
      }
      else
      {
        real_direction = direction;
      }
      PtbPsTree new_tree = create_bintree(*pos,daughters, real_direction, ancestor_name, mark, give_name);

      //replace content with binarized content
      add_left_sister(pos,new_tree);
      iterator new_pos = pos;
      new_pos.left();
      erase(pos);
      pos=new_pos;
    }
  }
}


struct Higher : public std::binary_function<const PtbPsTree::depth_first_iterator&,
const PtbPsTree::depth_first_iterator&,
bool>
{
    bool operator()(const PtbPsTree::depth_first_iterator& lhs,
                    const PtbPsTree::depth_first_iterator& rhs) const
    {
        return lhs->height() < rhs->height();
    }
};

void PtbPsTree::parent_annotate(unsigned level, bool annotate_pos)
{
  // initialise priority queue, sorts elements according to height in tree
  std::priority_queue<PtbPsTree::depth_first_iterator, std::vector<PtbPsTree::depth_first_iterator>, Higher> pqueue;
  for(iterator i = dfbegin(); i != dfend(); ++i)
    pqueue.push(i);

  // annotate tree nodes bottom up
  while(!pqueue.empty())
  {
    if(!pqueue.top()->leaf())
    {
      //      std::cout << "height " << pqueue.top()->height() << std::endl;

      iterator copy = pqueue.top();
      if (!copy.down_first()->leaf() || annotate_pos)
      {
        unsigned i = 0;
        iterator m = pqueue.top();
        while(i++ < level && m->has_mother()) {
          std::string mother_name = *(m.up());
          *(pqueue.top()) = *(pqueue.top()) + '^' + mother_name ;
        }
      }
    }
    pqueue.pop();
  }
}


void PtbPsTree::productions(std::vector<Production>& internals, std::vector<Production>& lexicals) const
{
  typedef PtbPsTree::const_depth_first_iterator const_iterator;
  static SymbolTable& sym_tab_nt = SymbolTable::instance_nt();
  static SymbolTable& sym_tab_word = SymbolTable::instance_word();

  for(const_iterator i(this->dfbegin()); i != this->dfend(); ++i)
  {
    if(!i->leaf()) {
      int lhs = sym_tab_nt.insert(*i);
      std::vector<int> rhs;

      const_iterator j = i;
      j.down_first();
      if(j->leaf())
      {  //lexical rule
        rhs.push_back(sym_tab_word.insert(*j));
        lexicals.push_back(Production(lhs,rhs,true));
      }
      else
      { // internal rule
        for(; j != this->dfend(); j.right())
          rhs.push_back(sym_tab_nt.insert(*j));

        internals.push_back(Production(lhs,rhs,false));
      }
    }
  }
}

void PtbPsTree::collect_internal_counts(std::map<Production, double> & binary_counts,
                                        std::map<Production, double> & unary_counts,
                                        std::map< int, double> & LHS_counts) const
{
    //static RandomGenerator random_gen(0,1);
    static SymbolTable& sym_tab_nt = SymbolTable::instance_nt();

    for(auto i = this->dfbegin(); i != this->dfend(); ++i)
    {
      if(!i->leaf())
      {
        auto j = i;
        j.down_first();

        double update = 1 /*+ random_gen.next() / 100*/ ;
        if(!j->leaf()) {//internal rule
          int lhs = sym_tab_nt.insert(*i);

          LHS_counts[lhs] += update;
          std::vector<int> rhs;

          for(; j != this->dfend(); j.right()){
            rhs.push_back(sym_tab_nt.insert(*j));
          }

          assert(rhs.size() >0 && rhs.size() <= 2);

          if (rhs.size()==1)//unary
            unary_counts[Production(lhs,rhs,false)] += update;
          else
          {
            binary_counts[Production(lhs,rhs,false)] += update;
          }
        }
        else
        {
          int lhs = sym_tab_nt.insert(*i);
          LHS_counts[lhs] += update;
        }
      }
    }
}





void PtbPsTree::anchored_productions(std::vector<std::tuple<int,int,int,Production>>& anchored_binaries,
                                     std::vector<std::tuple<int,int,Production>>& anchored_unaries,
                                     std::vector<std::tuple<int,Production>>& anchored_lexicals) const
{
  return anchored_productions(anchored_binaries, anchored_unaries, anchored_lexicals,
                              0, number_of_leaves());

}
void PtbPsTree::anchored_productions(std::vector<std::tuple<int,int,int,Production>>& anchored_binaries,
                                     std::vector<std::tuple<int,int,Production>>& anchored_unaries,
                                     std::vector<std::tuple<int,Production>>& anchored_lexicals,
                                     int begin, int end) const
{
  typedef PtbPsTree::const_depth_first_iterator const_iterator;
  static SymbolTable& sym_tab_nt = SymbolTable::instance_nt();
  static SymbolTable& sym_tab_word = SymbolTable::instance_word();

  auto i = dfbegin();

  if(!i->leaf()) {
    int lhs = sym_tab_nt.insert(*i);
    std::vector<int> rhs;

    const_iterator j = i;
    j.down_first();
    if(j->leaf())
    {  //lexical rule
      rhs.push_back(sym_tab_word.insert(*j));
      anchored_lexicals.push_back(std::make_tuple(begin,Production(lhs,rhs,true)));
    }
    else
    { // internal rule

      std::vector<PtbPsTree> daughters;

      for(; j != this->dfend(); j.right())
      {
        rhs.push_back(sym_tab_nt.insert(*j));
        auto t = subtree(j);
        daughters.push_back(*static_cast<PtbPsTree*>(&t));
      }
      if (rhs.size() == 1) //unary rule
      {
        anchored_unaries.push_back(std::make_tuple(begin,end,Production(lhs,rhs,false)));
        daughters[0].anchored_productions(anchored_binaries, anchored_unaries, anchored_lexicals,
                                          begin,end);
      }
      else //binaries
      {
        assert(rhs.size() == 2);

        int m = begin + daughters[0].number_of_leaves();


        anchored_binaries.push_back(std::make_tuple(begin,m,end,Production(lhs,rhs,false)));
        daughters[0].anchored_productions(anchored_binaries, anchored_unaries, anchored_lexicals,
                                          begin, m);
        daughters[1].anchored_productions(anchored_binaries, anchored_unaries, anchored_lexicals,
                                          m, end);
      }




    }
  }
}
