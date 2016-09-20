#include "LorgParseApp.h"

#include "lexicon/WordSignatureFactory.h"

LorgOptions LorgParseApp::get_options() const
{
  LorgOptions options;
  options.add_parser_options();
  options.add_unknown_word_mapping_options();
  //options.add_parser_positionals();
  return options;
}

bool LorgParseApp::read_config(ConfigTable& configuration)
{
  if(LorgApp::read_config(configuration) == false)
    return false;

 // get output stream or write to stdout
 if(configuration.exists("output")) {
   const std::string& output_filename = configuration.get_value<std::string>("output");
   std::ofstream * outstream = new std::ofstream(output_filename.c_str());
   if(*outstream) {
     out = outstream;
     if(verbose) std::clog << "Setting ouput file to " << output_filename << std::endl;
   }
   else {
     std::clog << "Could not open file " << output_filename << std::endl;
     return false;
   }
 }
 else {
   if(verbose) std::clog << "Writing output to stdout." << std::endl;
   out = &std::cout;
 }

  // get input stream or read from stdin
  if(configuration.exists("test-data")) {
    const std::string& input_filename = configuration.get_value<std::string>("test-data");
    std::ifstream * instream = new std::ifstream(input_filename.c_str());
    if(*instream) {
      in = instream;
      if(verbose) std::clog << "Setting test file to " << input_filename << std::endl;
    }
    else {
      std::clog << "Could not open file " << input_filename << std::endl;
      return false;
    }
  }
  else {
    if(verbose)	std::clog << "Reading input from stdin." << std::endl;
    in = &std::cin;
  }

  //get max_length
  max_length = configuration.get_value<unsigned>("max-length");


  //always_output_forms = configuration.get_value<bool>("always-output-forms");





  //creating tokeniser for english: TODO write tokeniser for other languages ?
  tokeniser = std::unique_ptr<Tokeniser>(
      TokeniserFactory::create_tokeniser());

  return true;
}
