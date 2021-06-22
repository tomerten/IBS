#include "NumericFunctions.hpp"
#include "RadiationDamping.hpp"
#include <map>
#include <math.h>
#include <stdio.h>
#include <string>

void twclog(double pnumber, double bx, double by, double dx, double dy,
            double ex, double ey, double r0, double gamma, double charge,
            double en0, double amass, double sige, double sigt, double *clog);

void twclogtail(double pnumber, double l, double bx, double by, double dx,
                double dpx, double dy, double dpy, double ax, double ay,
                double angle, double k1l, double k1sl, double ex, double ey,
                double r0, double aatom, double gamma, double en0, double len,
                double amass, double charge, double sige, double sigt,
                double *clog);

void CoulombLog(double pnumber, double ex, double ey,
                map<string, double> &twissheader, double sige, double sigt,
                double r0, bool printout, double *clog);

void TailCutCoulombLog(double pnumber, double ex, double ey,
                       map<string, double> &twissheader, double sige,
                       double sigt, double tauradx, double taurady,
                       double taurads, double r0, bool printout, double *clog);