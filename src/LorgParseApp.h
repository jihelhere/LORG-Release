// -*- mode: c++ -*-
#ifndef _LORGPARSEAPP_H_
#define _LORGPARSEAPP_H_

#include "LorgApp.h"


//should not be included here: defines TokMode
 #include "utils/TokeniserFactory.h"
#include "utils/Tagger.h"


class LorgParseApp : public LorgApp
{
public:
  LorgParseApp() : LorgApp(){} ;
  ~LorgParseApp() {};

protected:
  bool read_config(ConfigTable& configuration);
  LorgOptions get_options() const;

protected:
  std::auto_ptr<Tokeniser> tokeniser;

  unsigned max_length;

  std::string number_regex;
  bool replace_numbers;


  bool always_output_forms;
};

#endif /* _LORGPARSEAPP_H_ */
