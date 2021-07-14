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

#ifndef TWISS_H
#define TWISS_H

/**
 * Load Twiss header / summary  from file.
 *
 * @param filename Path to the Twiss file.
 * @return A map of twiss header parameters and their values.
 *
 * @note Twiss file needs to be produced with Madx.
 */
map<string, double> GetTwissHeader(string filename);

vector<vector<double>> GetTable(string filename, vector<string> columns);
map<string, vector<double>> GetTwissTableAsMap(string filename);

#endif