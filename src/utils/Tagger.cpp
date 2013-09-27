#include "Tagger.h"
#include "rules/MetaProduction.h"
#include "utils/ConfigTable.h"


#include <algorithm>


Tagger::Tagger(const std::vector< std::vector<const MetaProduction*> >* word_2_rule) :
    word_2_rule_(word_2_rule)
{
  switch(ConfigTable::access().get_value<std::string>("filter-tagging")[0])
  {
    case 's' :
      filter_tagging = Sancl;
      break;
    case 'b' :
      filter_tagging = Bikel;
      break;
    case 'n' :
      filter_tagging = None;
      break;
    default:
      throw std::out_of_range("filter tagging unknown");
  }
}


void Tagger::set_word_rules(const std::vector< std::vector<const MetaProduction*> >* word_2_rule)
{
  word_2_rule_ = word_2_rule;
}


void Tagger::tag(Word& word, const WordSignature& ws) const
{
  // only do this once ...
  static int unknown_id = SymbolTable::instance_word().insert(LorgConstants::token_unknown);
  static const std::vector<const MetaProduction*>& unknown_tags =  (*word_2_rule_)[unknown_id];

  static
  std::function<bool(Word&,int,const std::vector<const MetaProduction*>&)>
      find_tags_in_w2r_and_tag
      = [](Word& w, int tag, const std::vector<const MetaProduction*>& w2r)
      {
        bool result = false;

        auto r_iter =  std::find_if(w2r.begin(), w2r.end(),
                                    [&](const MetaProduction* mp)
                                       {return mp->get_lhs() == tag;}
                                    );

        if (r_iter != w2r.end())
        {
          w.rules.push_back(*r_iter);
          result = true;
        }

        return result;
      };


  word.initialize_id(ws);




  if(word.is_tagged()) {
    if((this->filter_tagging == Bikel && word.id == -1) or this->filter_tagging == None or this->filter_tagging == Sancl)
    {

      // for all the predicted tags for this word
      for (unsigned i = 0; i < word.tags.size(); ++i)
      {
        int given_tag = word.get_given_tag(i);
        bool found = false;
        if(word.id != -1)
        {
          // we try to find the first rule with the given tag in lhs
          // position
          found = find_tags_in_w2r_and_tag(word, given_tag, (*word_2_rule_)[word.id]);

          if(this->filter_tagging == Sancl)
            continue;
        }


        // if word is known but tag -> word rule is not found : look for
        // tag -> sig rule
        if(!found && word.id != -1 && word.sigid != word.id && word.sigid != -1)
        {
          found = find_tags_in_w2r_and_tag(word, given_tag, (*word_2_rule_)[word.sigid]);
        }

        //if we couldn't find TAG ->  signature
        //try to look for TAG -> UNKNOWN
        if(!found)
        {
          found = find_tags_in_w2r_and_tag(word, given_tag, (*word_2_rule_)[unknown_id]);
        }

        if(!found)
        {
          std::cerr << "Could not find assigned tag "
                    <<  SymbolTable::instance_nt().translate(word.tags[i])
                    << " for word " <<  word.form << std::endl;
        }
      }
    }
  }

  if (not word.is_tagged() or word.rules.empty())
  {
    //read tags from grammar
    //find all lexical rules with this word as its rhs
    //    std::cout << word.id << " : "  << word.form << std::endl;
    //    std::cout << word_2_rule_.size() << std::endl;
    if(word.id != -1)
      word.rules= (*word_2_rule_)[word.id];
    else
      word.rules= unknown_tags;

    //    std::cout << "nb tags :" << word.rules.size() << std::endl;
  }

  if(word.rules.empty()) {
    std::cerr << word.form << std::endl;
    throw std::logic_error("Can't tag word. Have you set unknown cutoff to zero ?");
  }
}


void Tagger::tag( std::vector< Word >& sentence, const WordSignature& ws ) const
{
  std::for_each(sentence.begin(), sentence.end(),
                [&](Word& w){this->tag(w,ws);}
                );
}
