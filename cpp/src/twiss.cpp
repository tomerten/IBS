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
map<string, double> GetTwissHeader(string filename) {
  vector<string> TWISSHEADERKEYS /* */ {
      "MASS",     "CHARGE",  "ENERGY",  "PC",      "GAMMA",   "KBUNCH",
      "BCURRENT", "SIGE",    "SIGT",    "NPART",   "EX",      "EY",
      "ET",       "BV_FLAG", "LENGTH",  "ALFA",    "ORBIT5",  "GAMMATR",
      "Q1",       "Q2",      "DQ1",     "DQ2",     "DXMAX",   "DYMAX",
      "XCOMAX",   "YCOMAX",  "BETXMAX", "BETYMAX", "XCORMS",  "YCORMS",
      "DXRMS",    "DYRMS",   "DELTAP",  "SYNCH_1", "SYNCH_2", "SYNCH_3",
      "SYNCH_4",  "SYNCH_5"};
  map<string, double> out;
  string line;
  ifstream file(filename);

  vector<int> labels;
  map<int, int> colmap;

  int counter = 0;
  int numberoflines = 1;

  getline(file, line);
  // check if file is open
  if (file.is_open()) {
    printf("File is open. Reading Twiss header.\n");
    // read lines until eof
    int counter = 0;
    string key;
    string at;
    double value;
    while (!file.eof() && counter < 41) {
      counter++;
      getline(file, line);
      istringstream iss(line);
      iss >> at >> key >> at >> value;
      vector<string>::iterator it =
          find(TWISSHEADERKEYS.begin(), TWISSHEADERKEYS.end(), key);
      // cout << key << " " << value << " " << endl;
      if (it != TWISSHEADERKEYS.end()) {
        out[key] = value;
      }
    }
  }
  file.close();
  printf("File is closed. Done reading Twiss header.\n");
  return out;
}

vector<vector<double>> GetTable(string filename, vector<string> columns) {
  string line;
  ifstream file(filename);

  vector<int> labels;
  map<int, int> colmap;

  int counter = 0;
  int numberoflines = 1;

  // check if file is open
  if (file.is_open()) {
    // printf("File is open\n");

    // read number of lines
    while (!file.eof()) {
      getline(file, line);
      ++numberoflines;
    }
    // cout << numberoflines << endl;
    // printf("Number of lines %i\n", numberoflines);
    file.clear();
    file.seekg(0, ios::beg);
    // getline(file, line);
    // cout << line << endl;
    vector<vector<double>> output(numberoflines - 50,
                                  vector<double>(columns.size()));

    // read lines until eof
    while (!file.eof()) {
      // increase line counter
      ++counter;
      // cout << counter << endl;

      // read a line
      getline(file, line);

      // check if eof
      if (file.eof())
        break;

      // if line is 47 read the column names
      if (counter == 47) {
        // load the current line as stream
        istringstream iss(line);

        // split the line and save in vector
        // vector<string> labels;
        int colcounter = 0;
        // cout << "Col idx: ";
        do {
          string sub;
          iss >> sub;
          vector<string>::iterator it =
              find(columns.begin(), columns.end(), sub);
          if (it != columns.end()) {
            colmap[colcounter - 1] = it - columns.begin();
            labels.push_back(colcounter - 1);
            // cout << colcounter << " ";
          }
          ++colcounter;
        } while (iss);
        // cout << endl;
      }

      if (counter > 48) {
        istringstream iss(line);
        int colcounter = 0;
        do {
          string sub;
          iss >> sub;
          vector<int>::iterator iti =
              find(labels.begin(), labels.end(), colcounter);
          if (iti != labels.end()) {
            // cout << counter << " " << counter - 49 << " " << sub;
            // output[counter - 49].push_back(stod(sub));
            output[counter - 49][colmap[colcounter]] = stod(sub);
          }
          ++colcounter;
        } while (iss);
        // cout << endl;
      }
    }
    return output;
  }
  vector<vector<double>> output(1, vector<double>(1));
  output[0][0] = 0.0;
  return output;
}

map<string, vector<double>> GetTwissTableAsMap(string filename) {
  vector<string> TWISSCOLS /* */ {"L",    "BETX", "ALFX", "BETY",
                                  "ALFY", "DX",   "DPX",  "DY",
                                  "DPY",  "K1L",  "K1SL", "ANGLE"};
  map<string, vector<double>> out;
  map<int, string> columnnames;

  string line;
  ifstream file(filename);

  if (file.is_open()) {
    printf("File is open\n");

    vector<double> row;
    int rowcounter = 0;
    while (!file.eof()) {
      // update row counter
      rowcounter++;

      // read a line
      getline(file, line);

      // if line is 47 read the column names
      if (rowcounter == 47) {
        // load the current line as stream
        istringstream iss(line);

        // split the line and save in vector
        // vector<string> labels;
        int colcounter = 0;
        // cout << "Col idx: ";
        do {
          string sub;
          iss >> sub;
          vector<string>::iterator it =
              find(TWISSCOLS.begin(), TWISSCOLS.end(), sub);
          if (it != TWISSCOLS.end()) {
            columnnames[colcounter - 1] = sub;
            // cout << colcounter - 1 << " " << sub << endl;
          }
          ++colcounter;
        } while (iss);
      }
      if (rowcounter > 48) {
        istringstream iss(line);
        int colcounter = 0;
        do {
          string sub;
          string key;
          double value;
          iss >> sub;
          if (columnnames.count(colcounter) > 0) {
            key = columnnames[colcounter];
            value = stod(sub);
            // cout << key << " " << value << endl;
            out[key].push_back(value);
          }
          ++colcounter;
        } while (iss);
        // cout << endl;
      }
    }
  }
  file.close();
  return out;
}
