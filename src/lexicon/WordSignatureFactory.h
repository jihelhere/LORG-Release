// -*- -mode: c++ *-
#ifndef _WORDSIGNATUREFACTORY_H_
#define _WORDSIGNATUREFACTORY_H_

#include "utils/ConfigTable.h"

#include "WordSignature.h"

namespace WordSignatureFactory {

  WordSignature *  create_wordsignature(unknown_word_map unknown_map, bool verbose);
  WordSignature *  create_wordsignature(ConfigTable& config);

}






#endif /* _WORDSIGNATUREFACTORY_H_ */
