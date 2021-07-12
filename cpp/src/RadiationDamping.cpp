#include "../include/ibs_bits/RadiationDamping.hpp"
#include "../include/ibs_bits/NumericFunctions.hpp"
#include <iostream>
#include <map>
#include <math.h>
#include <numeric>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

/*
================================================================================
================================================================================
METHOD TO PRINT RADIATION INTEGRALS FROM ARRAY.

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

================================================================================
  Arguments:
  ----------
    - double[] out
        array with radiation integrals.

  Returns:
  --------
    void

================================================================================
================================================================================
*/
void printradint(double out[7]) {
  printf("\n");
  printf("    Radiation Integrals \n");
  printf("    I1  = %15.6f\n", out[0]);
  printf("    I2  = %15.6f\n", out[1]);
  printf("    I3  = %15.6f\n", out[2]);
  printf("    I4x = %15.6f\n", out[3]);
  printf("    I4y = %15.6f\n", out[4]);
  printf("    I5x = %15.6f\n", out[5]);
  printf("    I5y = %15.6f\n", out[6]);
}
/*
================================================================================
================================================================================
METHOD TO CALCULATE RADIATION INTEGRALS USING SMOOTH LATTICE APPROXIMATION.

================================================================================
  AUTHORS:
    - TOM MERTENS
    - MIKE BLASKIEWISC
    - RODERIK BRUCE
    - MICHAELA SCHAUMANN

  HISTORY:
    - 04/02/2021 : UPDATE BUILD C LIBRARY FOR IBS
    - 08/06/2021 : initial cpp version (Tom)

================================================================================
  Arguments:
  ----------
    - double latticelength
        lattice length
    - double gamma
        relativistic gamma
    - double gammaTransition
        relativistic gamma transition
    - double dipoleBendingRadius
        dipole bending radius
    - double betax
        ring averaged beta x
    - double betay
        ring averaged beta y

  Returns:
  --------
    double [7] radiationIntegrals
      0 -> I1
      1 -> I2
      2 -> I3
      3 -> I4x
      4 -> I4y
      5 -> I5x
      6 -> I5y

================================================================================
================================================================================
*/

double *RadiationDampingApprox(double latticeLength, double gamma,
                               double gammaTransition,
                               double dipoleBendingRadius, double betax,
                               double betay) {
  static double radiationIntegrals[7];

  // Courant-Snyder optical functions
  double alphax = 0.0;
  double alphay = 0.0;
  double gammax = (1.0 + alphax * alphax) / betax;
  double gammay = (1.0 + alphay * alphay) / betay;

  // average dispersion
  double dx = latticeLength / (2.0 * pi * gammaTransition * gammaTransition);
  double dy = 0.0;
  double dpx = 0.1;
  double dpy = 0.0;

  // curly H
  double hx = betax * dpx + 2.0 * alphax * dx * dpx + gammax * dx;
  double hy = betay * dpy + 2.0 * alphay * dy * dpy + gammay * dy;

  // calculate radiation integrals
  radiationIntegrals[0] = dx / dipoleBendingRadius * latticeLength;
  radiationIntegrals[1] = 2.0 * pi / dipoleBendingRadius;
  radiationIntegrals[2] =
      2.0 * pi / (dipoleBendingRadius * dipoleBendingRadius);
  radiationIntegrals[3] = 0.0;
  radiationIntegrals[4] = 0.0;
  radiationIntegrals[5] =
      hx * 2.0 * pi / (dipoleBendingRadius * dipoleBendingRadius);
  radiationIntegrals[6] =
      hy * 2.0 * pi / (dipoleBendingRadius * dipoleBendingRadius);

  return radiationIntegrals;
}

/*
================================================================================
================================================================================
METHOD TO CALCULATE RADIATION INTEGRALS

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

================================================================================
  Arguments:
  ----------
    - map<string, vector<double>> table
        Twiss table after applying updateTwiss (NumericFunctions.hpp)

  Returns:
  --------
    double [7] radiationIntegrals
      0 -> I1
      1 -> I2
      2 -> I3
      3 -> I4x
      4 -> I4y
      5 -> I5x
      6 -> I5y


================================================================================
================================================================================
*/
double *RadiationDampingLattice(map<string, vector<double>> &table) {
  static double radiationIntegrals[7];
  radiationIntegrals[0] =
      accumulate(table["I1"].begin(), table["I1"].end(), 0.0);
  radiationIntegrals[1] =
      accumulate(table["I2"].begin(), table["I2"].end(), 0.0);
  radiationIntegrals[2] =
      accumulate(table["I3"].begin(), table["I3"].end(), 0.0);
  radiationIntegrals[3] =
      accumulate(table["I4x"].begin(), table["I4x"].end(), 0.0);
  radiationIntegrals[4] =
      accumulate(table["I4y"].begin(), table["I4y"].end(), 0.0);
  radiationIntegrals[5] =
      accumulate(table["I5x"].begin(), table["I5x"].end(), 0.0);
  radiationIntegrals[6] =
      accumulate(table["I5y"].begin(), table["I5y"].end(), 0.0);

  return radiationIntegrals;
}

/*
================================================================================
================================================================================
METHOD TO CALCULATE RADIATION LIFETIMES, RADIATION EQUILIBRIA AND PARTITION
NUMBERS.

================================================================================
  AUTHORS:
    - TOM MERTENS
    - MIKE BLASKIEWISC
    - RODERIK BRUCE
    - MICHAELA SCHAUMANN

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)
    - 11/06/2021 : removed I1 dependency and changed it for Synch Tune dep.

================================================================================
  Arguments:
  ----------
    - map<string, double> twissheadermap
        MADX twiss header map
    - double[7] radiationIntegrals
        radiation integrals array
    - double aatom
        atomic number of the particle
    - double qs
        synchrotron tune

  Returns:
  --------
    double [8] output
    0 -> alphax
    1 -> alphay
    2 -> alphas
    3 -> exinf
    4 -> eyinf
    5 -> sigEoE2
    6 -> sigs
    7 -> jx
    8 -> jy


================================================================================
================================================================================
*/

double *RadiationDampingLifeTimesAndEquilibriumEmittancesWithPartitionNumbers(
    map<string, double> &twissheadermap, double radiationIntegrals[7],
    double aatom, double qs) {

  const double c = clight;
  const double hbar = 1.0545718176461565e-34;
  const double electron_volt_joule_relationship = 1.602176634e-19;

  static double output[9];

  double gamma = twissheadermap["GAMMA"];
  double gammatr = twissheadermap["GAMMATR"];
  double p0 = twissheadermap["PC"] * 1.0e9;
  double len = twissheadermap["LENGTH"] * 1.0;
  double restE = twissheadermap["MASS"] * 1.0e9;
  double charge = twissheadermap["CHARGE"] * 1.0;
  double q1 = twissheadermap["Q1"];
  double i1 = radiationIntegrals[0];
  double i2 = radiationIntegrals[1];
  double i3 = radiationIntegrals[2];
  double i4x = radiationIntegrals[3];
  double i4y = radiationIntegrals[4];
  double i5x = radiationIntegrals[5];
  double i5y = radiationIntegrals[6];

  double particle_radius = charge * charge / aatom * 1.54e-18;

  double CalphaEC = particle_radius * c / (3.0 * restE * restE * restE) *
                    (p0 * p0 * p0 / len);

  // transverse partition numbers
  double jx = 1.0 - i4x / i2;
  double jy = 1.0 - i4y / i2;
  double alphax = 2.0 * CalphaEC * i2 * jx;
  double alphay = 2.0 * CalphaEC * i2 * jy;
  double alphas = 2.0 * CalphaEC * i2 * (jx + jy);

  // mc**2 expressed in Joule to match units of cq
  double mass = restE * electron_volt_joule_relationship;
  double cq = 55.0 / (32.0 * sqrt(3.0)) * (hbar * c) / mass;

  double sigE0E2 = cq * gamma * gamma * i3 / (2.0 * i2 + i4x + i4y);
  // ! = deltaE/E_0 see wiedemann p. 302,
  // and Wolski: E/(p0*c) - 1/beta0 = (E - E0)/(p0*c) = \Delta E/E0*beta0 with
  // E0 = p0*c/beta0 therefore:
  double betar = BetaRelativisticFromGamma(gamma);
  double dpop = dee_to_dpp(sqrt(sigE0E2), betar);
  // I1 calculated from twiss madx is to inaccurate
  // TODO: add option to get it from twiss header
  // double alfap = i1 / len;
  // double sigs = clight * alfap * sqrt(sigE0E2) / omegas;
  double sigs = dpop * len * eta(gamma, gammatr) / (2 * pi * qs);
  double exinf = cq * gamma * gamma * i5x / (jx * i2);
  double eyinf = cq * gamma * gamma * i5y / (jy * i2);

  double betaAvg = len / (q1 * 2.0 * pi);

  eyinf = (eyinf == 0.0) ? cq * betaAvg * i3 / (2.0 * jy * i2) : eyinf;

  output[0] = 1.0 / alphax;
  output[1] = 1.0 / alphay;
  output[2] = 1.0 / alphas;
  output[3] = exinf;
  output[4] = eyinf;
  output[5] = sigE0E2;
  output[6] = sigs;
  output[7] = jx;
  output[8] = jy;

  return output;
}

/*
================================================================================
================================================================================
METHOD TO CALCULATE RADIATION LOSSES PER TURN.

================================================================================
  AUTHORS:
    - TOM MERTENS
    - MIKE BLASKIEWISC
    - RODERIK BRUCE
    - MICHAELA SCHAUMANN

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

================================================================================
  Arguments:
  ----------
    - map<string, double> twiss
        madx twiss header map
    - double I2
        second radiation integral
    - double aatom
        atomic mass number of the particle

  Returns:
  --------
    double
      Energy loss per turn in eV.
================================================================================
================================================================================
*/

double RadiationLossesPerTurn(map<string, double> &twiss, double I2,
                              double aatom) {
  double gamma = twiss["GAMMA"];
  double p0 = twiss["PC"];
  double len = twiss["LENGTH"];
  double mass = twiss["MASS"];
  double charge = twiss["CHARGE"];

  double particle_radius = charge * charge / aatom * 1.54e-18;
  double cgamma = (4.0 * pi / 3.0) * (particle_radius / (mass * mass * mass));
  double betar = sqrt(1.0 - 1.0 / (gamma * gamma));
  double vrev = clight * betar;
  double trev = len / vrev;

  return (clight * cgamma) / (2.0 * pi * len) * p0 * p0 * p0 * p0 * I2 * 1.0e9 *
         trev;
}

/*
================================================================================
================================================================================
CALCULATE THE CRITICAL FREQUENCY AND ANGLE OF EMITTED RADIATION IN

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

================================================================================
  Arguments:
  ----------
    - double rho
        dipole bending radius
    - double gamma
        relativistic gamma
    - double omega
        angular frequency

  Retruns:
  --------
  double[5] output
    0 -> omega critical
    1 -> theta critical
    2 -> critical photon energy for given omega
    3 -> average energy per photon
    4 -> average number of photons emitted per turn
================================================================================
================================================================================
*/

double *RadiationCriticalEnergy(double rho, double gamma, double omega) {
  const double c = clight;
  const double alphafine = 7.297352569300000e-03;
  const double h = 6.626070150000000e-34;
  const double twoOthree = 2.0 / 3.0;
  const double gamma3 = gamma * gamma * gamma;

  static double output[5];

  output[0] = (1.0 / twoOthree) * c / rho * gamma3;
  output[1] = 1.0 / gamma * pow(output[0] / omega, 1.0 / 3.0);
  output[2] = h * output[0];
  output[3] = 1.0 / 3.0 * output[2];
  output[4] = 2 * pi * alphafine * gamma;

  return output;
}

/*
================================================================================
================================================================================
METHOD TO CALCULATE RADIATION INTEGRALS FOR SINGLE ELEMENT

================================================================================
  AUTHORS:
    - TOM MERTENS
    - MIKE BLASKIEWISC
    - RODERIK BRUCE
    - MICHAELA SCHAUMANN

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

================================================================================
  Arguments:
  ----------
    - double l
    - double bx
    - double by
    - double dx
    - double dpx
    - double dy
    - double dpy
    - double ax
    - double ay
    - double angle
    - double k1l
    - double k1sl
  Returns:
  --------
    double [6] radiationIntegrals
    0 -> I2
    1 -> I3
    2 -> I4x
    3 -> I4y
    4 -> I5x
    5 -> I5y

================================================================================
================================================================================
*/
double *RadiationDampingElement(double l, double bx, double by, double dx,
                                double dpx, double dy, double dpy, double ax,
                                double ay, double angle, double k1l,
                                double k1sl) {
  static double radiationIntegrals[6];
  double I2 = 0.0;
  double I3 = 0.0;
  double I4x = 0.0;
  double I4y = 0.0;
  double I5x = 0.0;
  double I5y = 0.0;

  double rhoi, ki, gammax, gammay, hx, hy;
  double rhoi2, rhoi3;
  // calculate local bending radius
  rhoi = (angle == 0.0) ? 0.0 : l / angle;
  rhoi2 = rhoi * rhoi;
  rhoi3 = rhoi2 * rhoi;

  // strength per length unit
  ki = (l == 0.0) ? 0.0 : k1l / l;

  // first for integrals
  I2 = (rhoi == 0.0) ? 0.0 : l / rhoi2;
  I3 = (rhoi == 0.0) ? 0.0 : l / rhoi3;
  I4x = (rhoi == 0.0)
            ? 0.0
            : (dx / rhoi3) * l + (2.0 / rhoi) * (ki * dx + k1sl * dy * l);
  I4y = 0.0;

  // Courant-Snyder gamma
  gammax = (1.0 + ax * ax) / bx;
  gammay = (1.0 + ay * ay) / by;

  // curly H
  hx = bx * dpx * dpx + 2.0 * ax * dx * dpx + gammax * dx * dx;
  hy = by * dpy * dpy + 2.0 * ay * dy * dpy + gammay * dy * dy;

  I5x = (rhoi == 0) ? 0.0 : hx * l / rhoi3;
  I5y = (rhoi == 0) ? 0.0 : hy * l / rhoi3;

  radiationIntegrals[0] = I2;
  radiationIntegrals[1] = I3;
  radiationIntegrals[2] = I4x;
  radiationIntegrals[3] = I4y;
  radiationIntegrals[4] = I5x;
  radiationIntegrals[5] = I5y;

  return radiationIntegrals;
}