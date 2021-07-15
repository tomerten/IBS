#include "CoulombLogFunctions.hpp"
#include "Integrators.hpp"
#include "NumericFunctions.hpp"
#include <iostream>
#include <map>
#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>

/**
 * Method that prints out IBS growth rates and times from given input array.
 *
 * @param output IBS growth rates array.
 */
void printouts(double output[3]);

/**
 * Piwinski smooth approximation model.
 *
 * @param pnumber numper of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param dponp @f$ \frac{dp}{p} @f$
 * @param twiss Twiss Header Map
 * @param r0 Classical particle radius
 *
 * @return IBS amplitude growth rates (longitudinal, horizontal, vertical)
 *
 * @note HANDBOOK FOR ACCELERATOR PHYSICISTS AND ENGINEERS P.126..
 */
double *PiwinskiSmooth(double pnumber, double ex, double ey, double sigs,
                       double dponp, map<string, double> &twiss, double r0);

/**
 * Piwinski approximation element weighted.
 *
 * @param pnumber numper of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param dponp @f$ \frac{dp}{p} @f$
 * @param twissheader Twiss Header Map
 * @param twissdata Twiss Table Map
 * @param r0 Classical particle radius
 *
 *
 * @return IBS amplitude growth rates (longitudinal, horizontal, vertical)
 * @note HANDBOOK FOR ACCELERATOR PHYSICISTS AND ENGINEERS P.126..
 */
double *PiwinskiLattice(double pnumber, double ex, double ey, double sigs,
                        double dponp, map<string, double> &twissheader,
                        map<string, vector<double>> &twissdata, double r0);
/**
 * Piwinski approximation element weighted taking also dispersion derivatives
 * into account.
 *
 * @param pnumber numper of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param dponp @f$ \frac{dp}{p} @f$
 * @param twissheader Twiss Header Map
 * @param twissdata Twiss Table Map
 * @param r0 Classical particle radius
 *
 * @return IBS amplitude growth rates (longitudinal, horizontal, vertical)
 * @note HANDBOOK FOR ACCELERATOR PHYSICISTS AND ENGINEERS P.126..
 *
 */
double *PiwinskiLatticeModified(double pnumber, double ex, double ey,
                                double sigs, double dponp,
                                map<string, double> &twissheader,
                                map<string, vector<double>> &twissdata,
                                double r0);
/**
 * Nagaitsev's high energy approximation element weighted with Coulomb logs
 * calculated element by element.
 *
 * @param pnumber numper of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param dponp @f$ \frac{dp}{p} @f$
 * @param twissheader Twiss Header Map
 * @param twissdata Twiss Table Map
 * @param r0 Classical particle radius
 *
 * @return IBS amplitude growth rates (longitudinal, horizontal, vertical)
 * @note PRSTAB 8, 064403 (2005
 *
 */
double *Nagaitsev(double pnumber, double ex, double ey, double sigs,
                  double dponp, map<string, double> &twissheader,
                  map<string, vector<double>> &twissdata, double r0);

/**
 * Nagaitsev's high energy approximation element weighted with Coulomb logs
 * calculated element by element and taking tailcut into account.
 *
 * @param pnumber numper of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param dponp @f$ \frac{dp}{p} @f$
 * @param twissheader Twiss Header Map
 * @param twissdata Twiss Table Map
 * @param r0 Classical particle radius
 * @param aatom Atomic Mass Number
 *
 * @return IBS amplitude growth rates (longitudinal, horizontal, vertical)
 * @note PRSTAB 8, 064403 (2005)
 *
 */
double *Nagaitsevtailcut(double pnumber, double ex, double ey, double sigs,
                         double dponp, map<string, double> &twissheader,
                         map<string, vector<double>> &twissdata, double r0,
                         double aatom);
/*
================================================================================

IBS MODELS USING DEDICATED INTEGRATORS (SEE INTEGRATORS)

================================================================================
*/
/**
 * Zimmerman model element weighted with Coulomb logs
 * calculated element by element using scaled Simpson Decade.
 *
 * @param pnumber numper of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param sige @f$ \frac{dE}{E} @f$
 * @param twissheader Twiss Header Map
 * @param twissdata Twiss Table Map
 * @param r0 Classical particle radius
 * @param printout Boolean flag to switch verbose mode on.
 *
 * @return IBS amplitude growth rates (longitudinal, horizontal, vertical)
 * @note CERN NOTE AB-2006--002
 *
 */
double *ibsmadx(double pnumber, double ex, double ey, double sigs, double sige,
                map<string, double> &twissheader,
                map<string, vector<double>> &twissdata, double r0,
                bool printout);

/**
 * Zimmerman model element weighted with Coulomb logs
 * calculated element by element taking tailcut into account using scaled
 * Simpson Decade.
 *
 * @param pnumber numper of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param sige @f$ \frac{dE}{E} @f$
 * @param twissheader Twiss Header Map
 * @param twissdata Twiss Table Map
 * @param r0 Classical particle radius
 * @param aatom Atomic Mass Number
 *
 * @return IBS amplitude growth rates (longitudinal, horizontal, vertical)
 * @note CERN NOTE AB-2006--002
 *
 */
double *ibsmadxtailcut(double pnumber, double ex, double ey, double sigs,
                       double sige, map<string, double> &twissheader,
                       map<string, vector<double>> &twissdata, double r0,
                       double aatom);
/**
 * Bjorken-Mtingwa model element weighted with Coulomb logs
 * calculated element by element using standard Simpson integration.
 *
 * @param pnumber numper of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param sige @f$ \frac{dE}{E} @f$
 * @param twissheader Twiss Header Map
 * @param twissdata Twiss Table Map
 * @param r0 Classical particle radius
 *
 * @return IBS amplitude growth rates (longitudinal, horizontal, vertical)
 * @note CERN NOTE AB-2006--002
 *
 */
double *BjorkenMtingwa2(double pnumber, double ex, double ey, double sigs,
                        double dponp, map<string, double> &twissheader,
                        map<string, vector<double>> &twissdata, double r0);
/**
 * Bjorken-Mtingwa model element weighted with Coulomb logs
 * calculated element by element using Simpson Decade integration.
 *
 * @param pnumber numper of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param dponp @f$ \frac{dp}{p} @f$
 * @param twissheader Twiss Header Map
 * @param twissdata Twiss Table Map
 * @param r0 Classical particle radius
 *
 * @return IBS amplitude growth rates (longitudinal, horizontal, vertical)
 * @note CERN NOTE AB-2006--002
 *
 */
double *BjorkenMtingwa(double pnumber, double ex, double ey, double sigs,
                       double dponp, map<string, double> &twissheader,
                       map<string, vector<double>> &twissdata, double r0);
/**
 * Bjorken-Mtingwa model element weighted with Coulomb logs
 * calculated element by element taking tailcut into account using Simpson
 * Decade integration.
 *
 * @param pnumber numper of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param dponp @f$ \frac{dp}{p} @f$
 * @param twissheader Twiss Header Map
 * @param twissdata Twiss Table Map
 * @param r0 Classical particle radius
 * @param aatom Atomic Mass Number
 *
 * @return IBS amplitude growth rates (longitudinal, horizontal, vertical)
 * @note CERN NOTE AB-2006--002
 *
 */
double *BjorkenMtingwatailcut(double pnumber, double ex, double ey, double sigs,
                              double dponp, map<string, double> &twissheader,
                              map<string, vector<double>> &twissdata, double r0,
                              double aatom);
/**
 * Conte-Martini model element weighted with Coulomb logs
 * calculated element by element using Simpson
 * Decade integration.
 *
 * @param pnumber numper of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param dponp @f$ \frac{dp}{p} @f$
 * @param twissheader Twiss Header Map
 * @param twissdata Twiss Table Map
 * @param r0 Classical particle radius
 *
 * @return IBS amplitude growth rates (longitudinal, horizontal, vertical)
 * @note CERN NOTE AB-2006--002
 *
 */
double *ConteMartini(double pnumber, double ex, double ey, double sigs,
                     double dponp, map<string, double> &twissheader,
                     map<string, vector<double>> &twissdata, double r0);
/**
 * Conte-Martini model element weighted with Coulomb logs
 * calculated element by element taking tailcut into account using Simpson
 * Decade integration.
 *
 * @param pnumber numper of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param dponp @f$ \frac{dp}{p} @f$
 * @param twissheader Twiss Header Map
 * @param twissdata Twiss Table Map
 * @param r0 Classical particle radius
 * @param aatom Atomic Mass Number
 *
 * @return IBS amplitude growth rates (longitudinal, horizontal, vertical)
 * @note CERN NOTE AB-2006--002
 *
 */
double *ConteMartinitailcut(double pnumber, double ex, double ey, double sigs,
                            double dponp, map<string, double> &twissheader,
                            map<string, vector<double>> &twissdata, double r0,
                            double aatom);
/**
 * Zimmerman model element weighted with Coulomb logs
 * calculated element by element using Simpson
 * Decade integration.
 *
 * @param pnumber numper of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param dponp @f$ \frac{dp}{p} @f$
 * @param twissheader Twiss Header Map
 * @param twissdata Twiss Table Map
 * @param r0 Classical particle radius
 *
 * @return IBS amplitude growth rates (longitudinal, horizontal, vertical)
 * @note CERN NOTE AB-2006--002
 *
 */
double *MadxIBS(double pnumber, double ex, double ey, double sigs, double dponp,
                map<string, double> &twissheader,
                map<string, vector<double>> &twissdata, double r0);