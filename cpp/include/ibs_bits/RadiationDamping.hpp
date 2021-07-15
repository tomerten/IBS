#include "NumericFunctions.hpp"
#include <iostream>
#include <map>
#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

/**
 * Method to print out radiation integrals.
 *
 * @param out Array output of Radation Damping Methods (Approx and Lattice)
 */
void printradint(double out[7]);

/**
 *
 * Method for calculating approximate radiation integrals. The approximation is
 due to using ring averages for the bending radius and beta functions.
 *
 * @param latticeLength Length of the accelerator lattice
 * @param gamma Relativistic gamma
 * @param gammaTransition Transition relativistic gamma
 * @param dipoleBendingRadius Main dipole bending radius
 * @param betax Ring averaged horizontal beta
 * @param betay Ring averaged vertical beta
 *
 * @return Approximate Radiation Damping integrals: @f$ I_1, I_2, I_3, I_{4x},
 * I_{4y}, I_{5x}, I_{5y} @f$
 *
 * @note Based on formulas in the Accelerator Handbook for Physicists and
 Engineers.
 */
double *RadiationDampingApprox(double latticeLength, double gamma,
                               double gammaTransition,
                               double dipoleBendingRadius, double betax,
                               double betay);
/**
 * Method for calculating element weighted radiation integrals. The
 integrals are calculated element by element and then weighted with the
 element's length.
 *
 * @param table Twiss Table data in Map format.
 * @see GetTwissTableAsMap
 *
 * @return Element weighted Radiation Damping integrals: @f$ I_1, I_2, I_3,
 I_{4x},
 * I_{4y}, I_{5x}, I_{5y} @f$
 *
 * @note Based on formulas in the Accelerator Handbook for Physicists and
 */
double *RadiationDampingLattice(map<string, vector<double>> &table);

/**
 * Method that calculates the radiation damping times, the radiation partition
 numbers and radiation and quantum excitation equilibria for horizontal
 emittance, vertical emittance, energy spread squared and bunchlength.
 *
 * @param twissheadermap Twiss Header Map
 * @param radiationIntegrals Output of the Radiation Damping Integrals Methods
 * @param aatom Atomic Mass Number
 * @param qs Synchrotron Tune
 *
 * @see RadiationDampingApprox, RadiationDampingLattice
 *
 * @return @f$ \tau_x, \tau_y, \tau_{dEoE2}, \epsilon_x^\infty,
 \epsilon_y^\infty, \sigma_E^2, \sigma_s, J_x, J_y @f$
 *
 * @note Based on formulas in the Accelerator Handbook for Physicists and
 */
double *RadiationDampingLifeTimesAndEquilibriumEmittancesWithPartitionNumbers(
    map<string, double> &twissheadermap, double radiationIntegrals[7],
    double aatom, double qs);

/**
 * Method to calculated the radiation losses per turn [eV].
 *
 * @param twiss Twiss Header Map
 * @param I2 Second radiation integral
 * @param aatom Atomic Mass Number
 *
 * @return Radiation Losses per turn [eV]
 *
 * @note Based on formulas in the Accelerator Handbook for Physicists and
 */
double RadiationLossesPerTurn(map<string, double> &twiss, double I2,
                              double aatom);

/**
 * CALCULATE THE CRITICAL FREQUENCY AND ANGLE OF EMITTED RADIATION
 *
 * @param rho Average bending radius
 * @param gamma Relativistic gamma
 * @param omega Angular revolution frequency
 *
 * @return Radiation data as an array containing
 *  0 -> omega critical
 *  1 -> theta critical
 *  2 -> critical photon energy for given omega
 *  3 -> average energy per photon
 *  4 -> average number of photons emitted per turn
 *
 * @note Based on Wolski Accelerator Lectures.
 */
double *RadiationCriticalEnergy(double rho, double gamma, double omega);

/**
 *
 * Method to calculate Radiation Damping integrals for single lattice element.
 *
 * @param l Length of the element
 * @param bx Horizontal beta function at the element
 * @param by Vertical beta function at the element
 * @param dx Horizontal Dispersion at the element
 * @param dpx Derivative of the horizontal dispersion at the element
 * @param dy Vertical dispersion at the element
 * @param dpy Derivative of the vertical dispersion at the element
 * @param ax Horizontal alpha
 * @param ay Vertical alpha
 * @param angle Bending angle
 * @param k1l Quadrupole component
 * @param k1sl Skew quadrupole component
 *
 * @return Radiation Damping integrals for single lattice element.
 *
 * @note Based on formulas in the Accelerator Handbook for Physicists and
 */
double *RadiationDampingElement(double l, double bx, double by, double dx,
                                double dpx, double dy, double dpy, double ax,
                                double ay, double angle, double k1l,
                                double k1sl);