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

/**
 * Load the Twiss table as vector of vectors.
 *
 * @param filename Path to the Twiss file.
 * @param columns Columns to extract from the Twiss file and store in the table.
 *
 * @return Table with the values of the selected Twiss columns.
 *
 * @note Twiss file needs to be produced with Madx.
 */
vector<vector<double>> GetTable(string filename, vector<string> columns);

/**
 * Loat the Twiss table as a map (columnNames -> values).
 *
 * @param filename Path to the Twiss file.
 *
 * @return Map with column names as keys and values are the respective Twiss
 values for that column.
 *
 * @note Twiss file needs to be produced with Madx.
 * @warning Pre-selection of columns is taken : {"L",    "BETX", "ALFX",
 "BETY","ALFY", "DX",   "DPX",  "DY", "DPY", "K1L", "K1SL", "ANGLE"}
 *
 */
map<string, vector<double>> GetTwissTableAsMap(string filename);

#endif