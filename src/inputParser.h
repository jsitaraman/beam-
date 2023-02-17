#include<unordered_map>
#include<vector>
#include <algorithm>
#include <unordered_map>
// Adapted from
// https://www.walletfox.com/course/parseconfigfile.php
std::unordered_map<std::string,std::vector<double>> readInputs(std::string fname)
{
  // std::ifstream is RAII, i.e. no need to call close
  std::unordered_map<std::string,std::vector<double>> inputMap;
  std::ifstream cFile (fname);
  if (cFile.is_open())
    {
      std::string line;
      while(getline(cFile, line)){
	line.erase(std::remove_if(line.begin(), line.end(), isspace),
		   line.end());
	if(line[0] == '#' || line.empty())
	  continue;
	auto delimiterPos = line.find("=");
	auto name = line.substr(0, delimiterPos);
	bool endReached=false;
	auto value = line.substr(delimiterPos + 1);
	while(!endReached) {
	  auto commaPos=value.find(",");
	  if ( commaPos!=-1) {
	    inputMap[name].push_back(std::stod(value.substr(0,commaPos)));
	    value=value.substr(commaPos+1);
	  } else {
	    inputMap[name].push_back(std::stod(value));
            endReached=true;
          }
	}
      }
    }
  else {
    std::cerr << "Couldn't open file" << fname << "for reading.\n";
  }
  return inputMap;
}
