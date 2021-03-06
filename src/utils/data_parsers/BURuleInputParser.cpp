#include "BURuleInputParser.h"
#include "utils/SymbolTable.h"

#include <fstream>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/phoenix/stl.hpp>
#include <boost/phoenix/object.hpp>
#include <boost/phoenix/bind/bind_member_function.hpp>
#include <boost/phoenix/bind/bind_function.hpp>

#include <boost/phoenix/fusion.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/variant/recursive_variant.hpp>

#include <boost/config/warning_disable.hpp>

#include <boost/fusion/include/std_pair.hpp>
#pragma clang diagnostic pop


namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;
namespace qi = boost::spirit::qi;
//namespace ascii = boost::spirit::ascii;
//namespace unicode = boost::spirit::unicode;
//namespace mychar = boost::spirit::ascii;
namespace mychar = boost::spirit::standard;

template <typename Iterator>
struct burule_parser : boost::spirit::qi::grammar<Iterator, std::vector<AnnotatedRule*>(), mychar::space_type>
{

burule_parser() : burule_parser::base_type(rules)
{
  using namespace qi::labels;

  symbol %= qi::lexeme[+(qi::char_ - ' ' - '\n')];

  nt_symbol= symbol [_val = phoenix::bind(&SymbolTable::insert,
					  phoenix::ref(SymbolTable::instance_nt()),
					  qi::labels::_1)];

  t_symbol = symbol [_val = boost::phoenix::bind(&SymbolTable::insert,
						 phoenix::ref(SymbolTable::instance_word()),
						 qi::labels::_1)];


  binary_proba_info_ = ('(' >> qi::short_ >> ','
			>> qi::short_ >> ','
			>> qi::short_ >> ','
			>> qi::long_double >> ')')
    [_val = phoenix::construct<binary_proba_info>(qi::labels::_1,
						  qi::labels::_2,
						  qi::labels::_3,
						  qi::labels::_4)];


  unary_proba_info_ = ('(' >> qi::int_ >> ','
		       >> qi::int_ >> ','
		       >> qi::long_double >> ')')
    [_val = phoenix::construct<unary_proba_info>(qi::labels::_1,
						 qi::labels::_2,
						 qi::labels::_3)];

  lexical_proba_info_ = ('(' >> qi::int_ >> ','
			 >> qi::long_double >> ')')
    [_val = phoenix::construct<lexical_proba_info>(qi::labels::_1,
						   qi::labels::_2)];



  binary = (qi::lit("int") >> nt_symbol >> nt_symbol >> nt_symbol
	    >> +binary_proba_info_  )
    [ _val = phoenix::new_<BRule>(qi::labels::_1,
				  qi::labels::_2,
				  qi::labels::_3,
				  qi::labels::_4)];

  unary_int = (qi::lit("int") >> nt_symbol >> nt_symbol
	       >> +unary_proba_info_ )
    [ _val = phoenix::new_<URule>(qi::labels::_1,
				  qi::labels::_2,
				  qi::labels::_3)];

  unary_lex = (qi::lit("lex") >> nt_symbol >> t_symbol
	       >> *lexical_proba_info_ )
    [ _val = phoenix::new_<LexicalRule>(qi::labels::_1,
					qi::labels::_2,
					qi::labels::_3)];

  rule %= unary_int | binary | unary_lex;

  rules %= +rule;

  rules.name("rules");
  rule.name("rule");
  unary_lex.name("unary_lex");
  unary_int.name("unary_int");
  binary.name("binary");
  nt_symbol.name("nt_symbol");
  symbol.name("symbol");
  t_symbol.name("t_symbol");
  binary_proba_info_.name("binary_proba_info");
  unary_proba_info_.name("unary_proba_info");
  lexical_proba_info_.name("lexical_proba_info");

  qi::on_error<qi::fail>
    (
     rule
     , std::clog
     << phoenix::val("Error! Expecting ")
     <<  qi::labels::_4                           // what failed?
     << phoenix::val(" here: \"")
     << phoenix::construct<std::string>( qi::labels::_3,
				qi::labels::_2)
     // iterators to error-pos, end
     << phoenix::val("\"")
     << std::endl
     );
}

  boost::spirit::qi::rule<Iterator, std::vector<AnnotatedRule*>(), mychar::space_type> rules;
  boost::spirit::qi::rule<Iterator, AnnotatedRule*(), mychar::space_type> rule;
  boost::spirit::qi::rule<Iterator, BRule*(), mychar::space_type> binary;
  boost::spirit::qi::rule<Iterator, URule*(), mychar::space_type> unary_int;
  boost::spirit::qi::rule<Iterator, LexicalRule*(), mychar::space_type> unary_lex;
  boost::spirit::qi::rule<Iterator, int(), mychar::space_type> nt_symbol, t_symbol;
  boost::spirit::qi::rule<Iterator, std::string(), mychar::space_type> symbol;
  boost::spirit::qi::rule<Iterator, binary_proba_info(), mychar::space_type> binary_proba_info_;
  boost::spirit::qi::rule<Iterator, unary_proba_info(), mychar::space_type> unary_proba_info_;
  boost::spirit::qi::rule<Iterator, lexical_proba_info(), mychar::space_type> lexical_proba_info_;
};

template <typename Iterator>
struct annotation_map_parser : boost::spirit::qi::grammar<Iterator, std::pair<int, Tree<unsigned> >(), mychar::space_type>
{
  annotation_map_parser() : annotation_map_parser::base_type(kvpair)
  {
    using namespace qi::labels;

    kvpair %= qi::lit("ainfos")  >> nt_symbol >> ptbtree;
    nt_symbol = symbol [_val = phoenix::bind(&SymbolTable::insert,
					     phoenix::ref(SymbolTable::instance_nt()),
					     qi::labels::_1)];
    symbol %= qi::lexeme[+(qi::char_ - ' ' - '\n')];


    ptbtree = '(' >> tree [_val = qi::labels::_1] >> ')';

    // tree --> ( nonterm [ term | tree+ ] )
    // if the parser finds a nonterm or term node, it will add a daughter to the current node
    // if the parser finds a closing bracket, it will make the mother of the current node the current node
    // this will happen after every terminal node as well

    tree
      =
      (id  [_val = phoenix::construct< Tree<unsigned> >(qi::labels::_1)] )
      |
      ( qi::char_('(') >>  id [_a = qi::labels::_1] >> +tree [phoenix::push_back(_b,qi::labels::_1)] >> qi::char_(')') [ _val = phoenix::construct< Tree<unsigned> >(_a,_b)])
      ;


    id %= qi::uint_;





    kvpair.name("kvpair");
    nt_symbol.name("nt_symbol");
    symbol.name("symbol");

    ptbtree.name("ptbtree");
    tree.name("tree");
    id.name("id");

    qi::on_error<qi::fail>(kvpair, std::clog << phoenix::val("Error!") << std::endl);
  }

  qi::rule<Iterator, std::pair<int, Tree<unsigned> >(),mychar::space_type> kvpair;
  qi::rule<Iterator, int(),mychar::space_type> nt_symbol;
  qi::rule<Iterator, std::string(),mychar::space_type> symbol;

  qi::rule<Iterator,  Tree<unsigned> (), mychar::space_type> ptbtree;
  qi::rule<Iterator,  Tree<unsigned> (), boost::spirit::locals< unsigned, std::vector< Tree<unsigned> > >, mychar::space_type> tree;
  qi::rule<Iterator, unsigned(), mychar::space_type> id;


};

//////////////////////

template <typename Iterator>
struct conf_parser : boost::spirit::qi::grammar<Iterator, std::pair<std::string,std::string>(), mychar::space_type>
{
  conf_parser() : conf_parser::base_type(kvpair)
  {
    using namespace qi::labels;

    kvpair %= qi::lit("// conf:")  >> symbol >> symbol ;
    symbol %= qi::lexeme[+(qi::char_ - ' ' - '\n')];

    kvpair.name("kvpair");
    symbol.name("symbol");

    qi::on_error<qi::fail>(kvpair, std::clog << phoenix::val("Error!") << std::endl);
  }

  qi::rule<Iterator, std::pair<std::string, std::string>(),mychar::space_type> kvpair;
  qi::rule<Iterator, std::string(),mychar::space_type> symbol;
};

/////////////////////



void BURuleInputParser::read_rulestring(const std::string& str, AnnotatedRule** rule_ptr) throw(ParseError)
{
  std::string::const_iterator iter = str.begin();
  std::string::const_iterator end = str.end();

  typedef std::string::const_iterator iterator_type;
  typedef burule_parser<iterator_type> parser;

  parser p;

  std::vector<AnnotatedRule*> v;

  bool r = phrase_parse(iter, end, p, mychar::space, v);

  if(!r) {
    //std::clog << "not read: " << std::string(iter,end) << std::endl;
    throw(ParseError());
  }

  *rule_ptr = v[0];

}

void BURuleInputParser::read_rulefile(const std::string& filename,
				      std::vector<LexicalRule>& lexicals,
				      std::vector<URule>& unaries,
				      std::vector<BRule>& n_aries,
                                      std::map<short, Tree<unsigned> >& history_tree_map,
                                      std::map<std::string, std::string>& conf
				      ) throw(ParseError)
{
  std::ifstream in_file(filename.c_str(),std::ios::in);
  if (!in_file)
  {
    std::cerr << "can't open file: " + filename  << std::endl;
    exit(1);
  }

  std::string str;

  in_file.unsetf(std::ios::skipws); // No white space skipping!
  std::copy(
	    std::istream_iterator<char>(in_file),
	    std::istream_iterator<char>(),
	    std::back_inserter(str));


  std::string::const_iterator iter = str.begin();
  std::string::const_iterator end = str.end();

  typedef std::string::const_iterator iterator_type;
  typedef burule_parser<iterator_type> parser;
  typedef annotation_map_parser<iterator_type> am_parser;
  typedef conf_parser<iterator_type> c_parser;

  parser p;
  am_parser a;

  c_parser c;

  bool res;


  do
  {
    std::pair<std::string, std::string> conf_pair;
    res = phrase_parse(iter, end, c, mychar::space, conf_pair);
    conf.insert(conf_pair);
  } while (res && iter != end);

  // Read annotations Map
  // lines like: Symbol History_Tree
  do {
    std::pair<int, Tree<unsigned> > m;
    res = phrase_parse(iter, end, a, mychar::space, m);
    history_tree_map.insert(m);
  } while(res && iter != end);


  std::vector<AnnotatedRule *> v;
  // Read rules
  do {
    //std::cout << "in rule parser" << std::endl;

    res = phrase_parse(iter, end, p, mychar::space, v);

    if(!res) {
      //std::cout << "not read: " << std::string(iter,end) << std::endl;
      throw(ParseError(filename));
    }

    // if(r->is_lexical()) {
    //   if(static_cast<LexicalRule*>(r)->get_probability().size()>0)
    //     lexicals.push_back(*(static_cast<LexicalRule*>(r)));
    // }
    // else
    //   if(r->is_unary()) {
    //     if(static_cast<URule*>(r)->get_probability().size()>0)
    //       unaries.push_back(*static_cast<URule*>(r));
    //   }
    //   else
    //     if(static_cast<BRule*>(r)->get_probability().size()>0)
    //       n_aries.push_back(*static_cast<BRule*>(r));

    // delete(r);
  } while(res && iter != end);


  lexicals.reserve(v.size());
  unaries.reserve(v.size());
  n_aries.reserve(v.size());
  for(auto& r : v)
  {
    if(r->is_lexical())
    {
      if(static_cast<LexicalRule*>(r)->get_probability().size()>0)
        lexicals.push_back(*(static_cast<LexicalRule*>(r)));
    }
    else
    {
      if(r->is_unary())
      {
        if(static_cast<URule*>(r)->get_probability().size()>0)
          unaries.push_back(*static_cast<URule*>(r));
      }
      else
        if(static_cast<BRule*>(r)->get_probability().size()>0)
          n_aries.push_back(*static_cast<BRule*>(r));
    }
    delete(r);
  }
  lexicals.shrink_to_fit();
  unaries.shrink_to_fit();
  n_aries.shrink_to_fit();

  //std::cout << "got here " << std::endl;
  in_file.close();
}
