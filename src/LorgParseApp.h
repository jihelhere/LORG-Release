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
  std::unique_ptr<Tokeniser> tokeniser;

  unsigned max_length;
};

#endif /* _LORGPARSEAPP_H_ */
