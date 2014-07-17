// -*- mode: c++ -*-
#ifndef PACKEDEDGEDAUGHTERS_H
#define PACKEDEDGEDAUGHTERS_H

#ifdef USE_THREADS
#include "utils/tick_count.h"
#endif

#include "rules/BRuleC2f.h"
#include "rules/URuleC2f.h"
#include "rules/LexicalRuleC2f.h"
#include "AnnotationInfo.h"

class Word;
//class PCKYAllCell;

/**
   \class PackedEdgeDaughters
   \brief represents a branching of possible daughters
*/

class PackedEdgeDaughters
{
protected:
  const AnnotatedRule * rule;
  bool calculated;
  double q_score;

public:
  PackedEdgeDaughters(const AnnotatedRule * r = 0)
      : rule(r), calculated(false), q_score(-std::numeric_limits<double>::infinity())
  {};
  ~PackedEdgeDaughters() {};

  inline bool is_unary() const {return get_rule()->is_unary();}
  inline bool is_binary() const {return get_rule()->is_binary();}
  inline bool is_lexical() const {return get_rule()->is_lexical();}

  inline bool is_calculated() const {return calculated;}
  inline void set_calculated(bool v) {calculated = v;}

  inline double get_q_score() const {return q_score;}
  inline void set_q_score(double v) {q_score = v;}

  // protected:
  inline const AnnotatedRule * get_rule() const {return rule;}
  inline void set_rule(const AnnotatedRule * r) {rule = r;}

};



template<class Rule>
class TypedRulePackedEdgeDaughters : public PackedEdgeDaughters
{
public:
  TypedRulePackedEdgeDaughters (const Rule * r) : PackedEdgeDaughters(r) {}
  inline const Rule * get_rule() const {return (Rule*) rule;}
};


/**
   \class BinaryPackedEdgeDaughters
   \brief represents a binary branching of possible daughters + a binary rule
*/
template<class Types>
class BinaryPackedEdgeDaughters : public TypedRulePackedEdgeDaughters<typename Types::BRule>
{
public:
  typedef TypedRulePackedEdgeDaughters<typename Types::BRule> Parent;
  typedef typename Types::Cell Cell;
  typedef typename Types::Edge Edge;
  typedef typename Types::Edge * EdgePtr ;
  typedef typename Types::Edge & EdgeRef ;
  typedef typename Types::BRule Rule;

protected:
  EdgePtr left;
  EdgePtr right;
  // EdgeRef left;
  // EdgeRef right;

  //   friend class Types::Edge ;
public:
  BinaryPackedEdgeDaughters(Edge& le, Edge& ri, const typename Types::BRule * ru)
      : Parent(ru), left(&le),right(&ri)
        //: Parent(ru), left(le),right(ri)
  {};

  // inline BinaryPackedEdgeDaughters& operator=(const BinaryPackedEdgeDaughters& other)
  // {
  //   if (this == &other)
  //   {
  //     return *this;
  //   }
  //   this->rule = other.rule;
  //   this->left = other.left;
  //   this->right = other.right;

  //   return *this;
  // }


  ~BinaryPackedEdgeDaughters() {};

  inline Edge& left_daughter() const {return *left;}
  inline Edge& right_daughter() const {return * right;}
  // inline Edge& left_daughter() const {return left;}
  // inline Edge& right_daughter() const  {return right;}

  inline bool operator==(const BinaryPackedEdgeDaughters& other)
  {
    return Parent::rule == other.rule && left == other.left && right == other.right;
  }
  inline bool points_towards_invalid_edges() const
  {
    return left->is_closed() or right->is_closed() ;
    //return left.is_closed() or right.is_closed() ;
  }

  inline void update_inside_annotations(AnnotationInfo & annotations) const {
    assert(Parent::rule != NULL);

    const auto& lannot = left->get_annotations();
    const auto& rannot = right->get_annotations();

    Parent::get_rule()->update_inside_annotations(annotations,
                                                  lannot,
                                                  rannot
                                                  );
  }

  inline void update_outside_annotations(AnnotationInfo & annotations) const
  {
    auto& lannot = left->get_annotations();
    auto& rannot = right->get_annotations();



    Parent::get_rule()->update_outside_annotations(annotations,
                                                   lannot,
                                                   rannot
                                                   );
  }
};


/**
   \class UnaryPackedEdgeDaughters
   \brief represents a unary branching (!) + a unary rule
*/
template<class Types>
class UnaryPackedEdgeDaughters : public TypedRulePackedEdgeDaughters<typename Types::URule>
{
public:
  typedef TypedRulePackedEdgeDaughters<typename Types::URule> Parent;
  typedef typename Types::Cell Cell;
  typedef typename Types::Edge Edge;
  typedef typename Types::Edge * EdgePtr ;
  typedef typename Types::URule Rule;

protected:
  EdgePtr left;

public:
  UnaryPackedEdgeDaughters(Edge & le, const Rule * ru)
      : Parent(ru), left(&le)
  {};

  ~UnaryPackedEdgeDaughters() {};

  inline bool is_binary() const {return false;}
  inline bool is_lexical() const {return false;}
  //   inline const Edge& left_daughter() const  {return *left;}
  inline Edge& left_daughter() const {return *left;}

  inline bool points_towards_invalid_edges() const
  {
    return left->is_closed();
  }
  inline void update_inside_annotations(AnnotationInfo & annotations) const {
    assert(Parent::rule != NULL);

    const auto& lannot = left->get_annotations();

    Parent::get_rule()->update_inside_annotations(annotations, lannot);
  }
  inline void update_outside_annotations(AnnotationInfo & annotations) const
  {
    auto& lannot = left->get_annotations();

    Parent::get_rule()->update_outside_annotations(annotations, lannot);
  }
};


/**
   \class LexicalPackedEdgeDaughters
   \brief represents a zero-ary branching (!) + a lexical rule
*/
template<class Types>
class LexicalPackedEdgeDaughters : public TypedRulePackedEdgeDaughters<typename Types::LRule>
{
public:
  typedef typename Types::LRule Rule;
  typedef TypedRulePackedEdgeDaughters<typename Types::LRule> Parent;

  protected:
  const Word* word;

public:
  LexicalPackedEdgeDaughters(const Rule * ru, const Word* w) :
  Parent(ru), word(w)
  {};

  ~LexicalPackedEdgeDaughters() {};

  inline bool is_binary() const {return false;}
  inline bool is_lexical() const {return true;}
  inline const Word* get_word() const {return word;}

  inline void update_inside_annotations(AnnotationInfo & annotations) const {
//     BLOCKTIMING("LexicalPackedEdgeDaughters - update_inside_annotations");
    assert(Parent::get_rule() != NULL);
    // Parent::get_rule()->update_inside_annotations(annotations.inside_probabilities.array,
    //                                               annotations.invalids);
    Parent::get_rule()->update_inside_annotations(annotations);
  }
};


#endif //PACKEDEDGEDAUGHTERS
