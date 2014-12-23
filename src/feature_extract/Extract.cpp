#include "Extract.h"

//#include <functional>
//#include <algorithm>
#include <iostream>
#include <sstream>


void Extract::extract(const PtbPsTree& tree, std::string& result)
{
  std::map<std::string, unsigned> counts;


  for(auto d = tree.dfbegin(); d != tree.dfend(); ++d) {

    for (size_t i = 0; i < features.size(); ++i)
      {
        std::string s;
        if((s = features[i].extract_with_name(d,tree)) != "")
        counts[s]++;
      }

  }

  std::stringstream ss(std::stringstream::in |std::stringstream::out);
  for (auto i = counts.begin(); i != counts.end(); ++i)
    {
       ss << i->first << ":" << i->second << " ";
    }

  result = ss.str();

}
