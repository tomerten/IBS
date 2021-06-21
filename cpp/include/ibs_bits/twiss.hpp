#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <sstream>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

map<string, double> GetTwissHeader(string filename);
vector<vector<double>> GetTable(string filename, vector<string> columns);
map<string, vector<double>> GetTwissTableAsMap(string filename);