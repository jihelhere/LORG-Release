#include "TokeniserFactory.h"
#include "utils/EnglishSpec.h"

#include "utils/ConfigTable.h"


Tokeniser * TokeniserFactory::create_tokeniser()

{
  const ConfigTable& configuration = ConfigTable::access();
  bool verbose = configuration.exists("verbose");

  // get rm_punctuation
  bool rm_punctuation = configuration.exists("remove-punctuation");
  if(rm_punctuation) {
    if(verbose) std::clog << "Removing punctuation from input.\n";
  }
  else {
    if(verbose) std::clog << "Not removing punctuation from input.\n";
  }



  // get input mode
  TokMode input_mode;
  try {
    input_mode = Tokeniser::string_to_tokmode(configuration.get_value<std::string>("input-mode"));
  }
  catch(std::out_of_range& e) { // should be a proper exception here
    std::clog << "Unknown input mode: " <<  configuration.get_value<std::string>("input-mode") << std::endl;
    throw e;
  }

  if(verbose)
    std::clog << "Input mode set to: " << Tokeniser::tokmode_to_string(input_mode) << std::endl;


  char comment_char = configuration.get_value<char>("comment-char");

  std::string tokens_to_remove = "";
  std::string delim = " ";

  return new Tokeniser(new EnglishSpec(), rm_punctuation, input_mode,
                       tokens_to_remove, delim, comment_char);

}
