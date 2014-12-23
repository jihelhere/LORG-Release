#include "ExtractApp.h"
#include "utils/data_parsers/PTBInputParser.h"


ExtractApp::ExtractApp() : LorgApp(),  extractor()
{
  in = &std::cin;
  out= &std::cout;

  extractor = std::make_shared<Extract>
              (Extract({
                  RuleFeature(),
                      RuleParentFeature(),
                      RuleGrandParentFeature(),
                      BiGramNodeFeature(),
                      HeavyFeature(),
                      NeighboursFeature(),
                      NeighboursExtFeature(),
                      WordFeature2(),
                      WordFeature3(),
                      WordFeatureGen2(),
                      WordFeatureGen3()}));
  }

LorgOptions ExtractApp::get_options() const
{
  LorgOptions options;
  return options;
}


int ExtractApp::run()
{
  std::string raw_string;

  while(std::getline((*in),raw_string)) {

    if(raw_string == "" || raw_string == "no parse found" || raw_string == "(())")
      (*out) << std::endl;
    else {
      std::string result;
      PtbPsTree tree = PTBInputParser::from_string(raw_string)[0]; // assume one tree per line
      tree.remove_function(); // FIXME : make it a command-line switch
      extractor->extract(tree, result);
      (*out) << result << std::endl;
    }
  }

  (*out) << std::endl;
  return 0;
}

int main(int, char **)
{
  ExtractApp ea;
  ea.run();

  return 0;
}
