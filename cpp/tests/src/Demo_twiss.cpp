#include <ibs>
#include <iostream>
#include <map>
#include <stdio.h>
#include <string>
#include <vector>

int main() {
  string twissfilename = "../src/b2_design_lattice_1996.twiss";
  map<string, double> twissheadermap;
  twissheadermap = GetTwissHeader(twissfilename);

  printf("Test %12.6e\n", twissheadermap["GAMMA"]);
  return 0;
}