#include <iostream>
#include <fstream>
#include "inputParser.h"
#include "beam.h"
using namespace beam_pp;
int main(int argc, char *argv[])
{
  auto inputMap=readInputs(std::string(argv[1]));
  for(auto &val : inputMap) {
    std::cout << val.first  ;
    for(auto &datum : val.second) {
      std::cout << " " << datum;
    }
    std::cout << "\n";
  }
  auto beam1=EulerBernoulliBeam(inputMap);
}
