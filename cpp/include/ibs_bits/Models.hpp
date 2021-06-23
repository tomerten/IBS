#include "CoulombLogFunctions.hpp"
#include "Integrators.hpp"
#include "NumericFunctions.hpp"
#include <iostream>
#include <map>
#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>

void printouts(double output[3]);

double *PiwinskiSmooth(double pnumber, double ex, double ey, double sigs,
                       double dponp, map<string, double> &twiss, double r0);

double *PiwinskiLattice(double pnumber, double ex, double ey, double sigs,
                        double dponp, map<string, double> &twissheader,
                        map<string, vector<double>> &twissdata, double r0);

double *PiwinskiLatticeModified(double pnumber, double ex, double ey,
                                double sigs, double dponp,
                                map<string, double> &twissheader,
                                map<string, vector<double>> &twissdata,
                                double r0);

double *Nagaitsev(double pnumber, double ex, double ey, double sigs,
                  double dponp, map<string, double> &twissheader,
                  map<string, vector<double>> &twissdata, double r0);

double *Nagaitsevtailcut(double pnumber, double ex, double ey, double sigs,
                         double dponp, map<string, double> &twissheader,
                         map<string, vector<double>> &twissdata, double r0,
                         double aatom);
/*
================================================================================

IBS MODELS USING DEDICATED INTEGRATORS (SEE INTEGRATORS)

================================================================================
*/
double *ibsmadx(double pnumber, double ex, double ey, double sigs, double sige,
                map<string, double> &twissheader,
                map<string, vector<double>> &twissdata, double r0,
                bool printout);

double *ibsmadxtailcut(double pnumber, double ex, double ey, double sigs,
                       double sige, map<string, double> &twissheader,
                       map<string, vector<double>> &twissdata, double r0,
                       double aatom);

double *BjorkenMtingwa2(double pnumber, double ex, double ey, double sigs,
                        double dponp, map<string, double> &twissheader,
                        map<string, vector<double>> &twissdata, double r0);

double *BjorkenMtingwa(double pnumber, double ex, double ey, double sigs,
                       double dponp, map<string, double> &twissheader,
                       map<string, vector<double>> &twissdata, double r0);

double *BjorkenMtingwatailcut(double pnumber, double ex, double ey, double sigs,
                              double sige, map<string, double> &twissheader,
                              map<string, vector<double>> &twissdata, double r0,
                              double aatom);

double *ConteMartini(double pnumber, double ex, double ey, double sigs,
                     double dponp, map<string, double> &twissheader,
                     map<string, vector<double>> &twissdata, double r0);

double *ConteMartinitailcut(double pnumber, double ex, double ey, double sigs,
                            double dponp, map<string, double> &twissheader,
                            map<string, vector<double>> &twissdata, double r0,
                            double aatom);

double *MadxIBS(double pnumber, double ex, double ey, double sigs, double dponp,
                map<string, double> &twissheader,
                map<string, vector<double>> &twissdata, double r0);