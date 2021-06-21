#include "NumericFunctions.hpp"
#include <iostream>
#include <map>
#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

void printradint(double out[7]);

double *RadiationDampingApprox(double latticeLength, double gamma,
                               double gammaTransition,
                               double dipoleBendingRadius, double betax,
                               double betay);

double *RadiationDampingLattice(map<string, vector<double>> &table);

double *RadiationDampingLifeTimesAndEquilibriumEmittancesWithPartitionNumbers(
    map<string, double> &twissheadermap, double radiationIntegrals[7],
    double aatom, double qs);

double RadiationLossesPerTurn(map<string, double> &twiss, double I2,
                              double aatom);

double *RadiationCriticalEnergy(double rho, double gamma, double omega);

double *RadiationDampingElement(double l, double bx, double by, double dx,
                                double dpx, double dy, double dpy, double ax,
                                double ay, double angle, double k1l,
                                double k1sl);