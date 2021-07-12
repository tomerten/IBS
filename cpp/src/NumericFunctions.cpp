#include "../include/ibs_bits/NumericFunctions.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdio.h>
#include <string>
#include <vector>
using namespace std;
/*
================================================================================
================================================================================
METHOD TO CALCULATE dE/E FROM sigma_s.

REMARK:
-------
  dE/E (do not confuse with dp/p, there is a factor beta**2) De/E = beta**2 dp/p

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 16/02/2021 : initial version (Tom)

================================================================================
  Arguments:
  ----------
    - double : qs
        longitunal tune
    - double : sigs
        bunch length [m]
    - double : omega0
        angular RF frequency
    - double : gamma
        relativistic gamma
    - double : gammatr
        relativistic gamma transition

  Returns:
  --------
    double
      dE/E

================================================================================
================================================================================
*/
double sigefromsigs(double omega0, double sigs, double qs, double gamma,
                    double gammatr) {
  // dE/E = Beta**2 dp/p
  double beta2 = BetaRelativisticFromGamma(gamma);
  beta2 *= beta2;
  return beta2 * qs * omega0 * (sigs / (fabs(eta(gamma, gammatr)) * clight));
}

/*
================================================================================
================================================================================
METHOD TO CALCULATE sigma_s FROM dE/E.

REMARK:
-------
  dE/E (do not confuse with dp/p, there is a factor beta**2) De/E = beta**2 dp/p

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial version (Tom)

================================================================================
  Arguments:
  ----------
    - double : sige
        dE/E
    - double : gamma
        relativistic gamma
    - double : gammatr
        relativistic gamma for transition
    - double : omegas
        synchrotron angular frequency (omega0 * qs)


  Returns:
  --------
    double
      sigma_s

================================================================================
================================================================================
*/

double sigsfromsige(double sige, double gamma, double gammatr, double omegas) {
  // dE/E = Beta**2 dp/p
  double beta2 = BetaRelativisticFromGamma(gamma);
  beta2 *= beta2;
  return clight * fabs(eta(gamma, gammatr)) / omegas * sige / beta2;
}
/*
================================================================================
================================================================================
METHOD TO CALCULATE PHASE SLIP FACTOR ETA.

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial version (Tom)

================================================================================
  Arguments:
  ----------
    - double : gamma
        relativistic gamma
    - double : gammatr
        relativistic gamma for transition

  Returns:
  --------
    double
      eta - phase slip factor (1/ gammatr**2  - 1/ gamma**2)

================================================================================
================================================================================
*/

double eta(double gamma, double gammatr) {
  return 1.0 / (gammatr * gammatr - 1.0 / (gamma * gamma));
}

/*
================================================================================
================================================================================
METHOD TO CALCULATE FMOHL.

================================================================================
  AUTHORS:
    - TOM MERTENS
    - MIKE BLASKIEWISC
    - RODERIK BRUCE
    - MICHAELA SCHAUMANN

  HISTORY:
    - 05/02/2021 COPYRIGHT : CERN / HZB
    - 08/06/2021 : initial cpp version (Tom)

  REFS:
    - HANDBOOK FOR ACCELERATOR PHYSICISTS AND ENGINEERS
    - BANE: A SIMPLIFIED MODEL OF INTRABEAM SCATTERING (EPAC 2002)
    - PIWINSKI  Tech. Rep. HEAC 74, Stanford, 1974.

================================================================================
  Arguments:
  ----------
    - double : a
        parameter 1
    - double : b
        parameter 2
    - double : q
        parameter 3
    - int : n

  Returns:
  --------
    double
      Fmohl

================================================================================
================================================================================
*/
double fmohl(double a, double b, double q, int n) {

  double u, cp, cq;
  double sum = 0.0;
  double const du = 1.0 / n; // careful needs to be 1.0 not 1

#pragma omp parallel for reduction(+ : sum)
  for (int i = 0; i <= n; i++) {
    double dsum = 0.0;
    u = i * du;
    cp = sqrt(a * a + (1.0 - a * a) * u * u);
    cq = sqrt(b * b + (1.0 - b * b) * u * u);
    dsum = 2.0 * log(q * (1.0 / cp + 1.0 / cq) / 2.0) - euler;
    dsum *= (1.0 - 3.0 * u * u) / (cp * cq);
    dsum = (i == 0) ? dsum /= 2.0 : dsum;
    dsum = (i == n) ? dsum /= 2.0 : dsum;
    sum += dsum;
  }

  sum *= 8 * pi * du;

  return sum;
}

/*
================================================================================
================================================================================
METHOD TO CALCULATE CLASSICAL PARTICLE RADIUS.

================================================================================
  AUTHORS:
    - TOM MERTENS
    - MIKE BLASKIEWISC
    - RODERIK BRUCE
    - MICHAELA SCHAUMANN

  HISTORY:
    - 05/02/2021 COPYRIGHT : CERN / HZB
    - 08/06/2021 : initial cpp version (Tom)

  REF:
    HANDBOOK FOR ACCELERATOR PHYSICISTS AND ENGINEERS (p. )

================================================================================
  Arguments:
  ----------
    - double charge
        particle charge
    - double aatom
        atomic number (A)

  Returns:
  --------
    double
      classical particle radius

================================================================================
================================================================================
*/

double ParticleRadius(double charge, double aatom) {
  return charge * charge / aatom * prad;
}

/*
================================================================================
================================================================================
METHOD TO CALCULATE RELATIVISTIC BETA FROM RELATIVISTIC GAMMA

================================================================================
  AUTHORS:
    - TOM MERTENS
    - MIKE BLASKIEWISC
    - RODERIK BRUCE
    - MICHAELA SCHAUMANN

  HISTORY:
    - 05/02/2021 COPYRIGHT : CERN / HZB
    - 08/06/2021 : initial cpp version (Tom)

================================================================================
  Arguments:
  ----------
    - double gamma
        relativistic gamma

  Returns:
  --------
    double
      relativistic beta

================================================================================
================================================================================
*/
double BetaRelativisticFromGamma(double gamma) {
  return sqrt(1 - (1 / (gamma * gamma)));
}

/*
================================================================================
================================================================================
MAX AND MIN METHODS

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF :
    https://stackoverflow.com/questions/3437404/min-and-max-in-c

================================================================================
================================================================================
*/

#define max(a, b)                                                              \
  ({                                                                           \
    __typeof__(a) _a = (a);                                                    \
    __typeof__(b) _b = (b);                                                    \
    _a > _b ? _a : _b;                                                         \
  })

#define min(a, b)                                                              \
  ({                                                                           \
    __typeof__(a) _a = (a);                                                    \
    __typeof__(b) _b = (b);                                                    \
    _a < _b ? _a : _b;                                                         \
  })

/*
================================================================================
================================================================================
RDS FUNCTION NAGAITSEV IBS PAPER.

================================================================================
  AUTHORS:
    - TOM MERTENS
    - MIKE BLASKIEWISC
    - RODERIK BRUCE
    - MICHAELA SCHAUMANN

  HISTORY:
    - 05/02/2021 COPYRIGHT : CERN / HZB
    - 08/06/2021 : initial cpp version (Tom)

  REF:
      PRSTAB 8, 064403 (2005)
================================================================================
  Arguments:
  ----------
    - double x
    - double y
    - double z

  Returns:
  --------
    double
      rds(x,y,z)

================================================================================
================================================================================
*/

double rds(double x, double y, double z) {
  // init
  double errtol = 0.05;
  // double tiny   = 1.0e-25;
  // double big    = 4.5e21;
  double c1 = 3.0 / 14.0;
  double c2 = 1.0 / 6.0;
  double c3 = 9.0 / 22.0;
  double c4 = 3.0 / 26.0;
  double c5 = 0.25 * c3;
  double c6 = 1.50 * c4;

  double xt = x;
  double yt = y;
  double zt = z;
  double sum = 0.0;
  double fac = 1.0;
  int iter = 0;

  double sqrtx, sqrty, sqrtz, alamb;
  double ave, delx, dely, delz;

  do {
    iter = iter + 1;
    sqrtx = sqrt(xt);
    sqrty = sqrt(yt);
    sqrtz = sqrt(zt);

    alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
    sum = sum + fac / (sqrtz * (zt + alamb));

    fac = 0.25 * fac;
    xt = 0.25 * (xt + alamb);
    yt = 0.25 * (yt + alamb);
    zt = 0.25 * (zt + alamb);
    ave = 0.20 * (xt + yt + 3.0 * zt);
    delx = (ave - xt) / ave;
    dely = (ave - yt) / ave;
    delz = (ave - zt) / ave;
  } while (max(max(fabs(delx), fabs(dely)), fabs(delz)) >= errtol);

  double ea = delx * dely;
  double eb = delz * delz;
  double ec = ea - eb;
  double ed = ea - 6.0 * eb;
  double ee = ed + ec + ec;

  return 3.0 * sum +
         fac *
             (1.0 + ed * (-c1 + c5 * ed - c6 * delz * ee) +
              delz * (c2 * ee + delz * (-c3 * ec + delz * c4 * ea))) /
             (ave * sqrt(ave));
}

/*
================================================================================
================================================================================
METHOD RETURNS EFFECTIVE RF VOLTAGE TIMES CHARGE (eV)

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF:
      USPAS ACCELERATOR COURSE
      LEE
================================================================================
  Arguments:
  ----------
    - double phi
        RF phase in deg
    - double charge
        particle charge
    - double nrf
        number of RF systems
    - doube[] harmon
        corresponding harmonic numbers
    - double[] voltages
        corresponding voltages

  Returns:
  --------
    double
      Effective Voltage times particle charge - units (eV)

================================================================================
================================================================================
*/

double EffectiveRFVoltageInElectronVolt(double phi, double charge, int nrf,
                                        double harmon[], double voltages[]) {
  // init
  double vrf = voltages[0] * sin(phi / 180.0 * pi);

  // add the rest taking harmonic numbers into account
  for (int i = 1; i < nrf; i++) {
    vrf += voltages[i] * sin((harmon[i] / harmon[0]) * phi / 180.0 * pi);
  }

  // multiply with charge
  vrf *= charge;

  return vrf;
}

/*
================================================================================
================================================================================
METHOD RETURNS DERIVATIVE WITH RESPECT TO PHASE OF EFFECTIVE RF VOLTAGE
TIMES CHARGE (eV).

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF:
      USPAS ACCELERATOR COURSE
      LEE
================================================================================
  Arguments:
  ----------
    - double phi
        RF phase
    - double charge
        particle charge
    - double nrf
        number of RF systems
    - doube[] harmon
        corresponding harmonic numbers
    - double[] voltages
        corresponding voltages

  Returns:
  --------
    double
      Effective Voltage times particle charge - units (eV)

================================================================================
================================================================================
*/
double EffectiveRFVoltageInElectronVoltPrime(double phi, double charge, int nrf,
                                             double harmon[],
                                             double voltages[]) {
  // init
  double vrf = voltages[0] * cos(phi / 180.0 * pi);

  // add other rfs
  for (int i = 1; i < nrf; i++) {
    vrf += voltages[i] * (harmon[i] / harmon[0]) *
           cos((harmon[i] / harmon[0]) * phi / 180.0 * pi);
  }

  // V -> eV
  vrf *= charge;

  return vrf;
}

/*
================================================================================
================================================================================
METHOD RETURNS EFFECTIVE RF VOLTAGE TIMES CHARGE MINUS ENERGY LOSS PER TURN.
ALL UNITS ARE IN eV.

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF:
      USPAS ACCELERATOR COURSE
      LEE
================================================================================
  Arguments:
  ----------
    - double phi
        RF phase
    - double U0
        Energy loss per turn in eV
    - double charge
        particle charge
    - double nrf
        number of RF systems
    - doube[] harmon
        corresponding harmonic numbers
    - double[] voltages
        corresponding voltages

  Returns:
  --------
    double
      Effective Voltage times particle charge minus energy loss per turn
      units (eV).

================================================================================
================================================================================
*/
double VeffRFeVRadlosses(double phi, double U0, double charge, int nrf,
                         double harmon[], double voltages[]) {
  double vrf =
      EffectiveRFVoltageInElectronVolt(phi, charge, nrf, harmon, voltages) - U0;

  return vrf;
}

/*
================================================================================
================================================================================
METHOD TO GET THE SYNCHRONUOUS PHASE.

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF :
    https://www.quantstart.com/articles/Implied-Volatility-in-C-using-Template-Functions-and-Newton-Raphson/

================================================================================
  Arguments:
  ----------
    - double target
        target value for the test function (zero)
    - double init_phi
        initial phase to start searching from
    - double U0
        Energy loss per turn in eV
    - double charge
        particle charge
    - double nrf
        number of RF systems
    - doube[] harmon
        corresponding harmonic numbers
    - double[] voltages
        corresponding voltages
    - double epsilon
        threshold accuracy for the return value

  Returns:
  --------
    double
      synchronuous phase

================================================================================
================================================================================
*/

double SynchronuousPhase(double target, double init_phi, double U0,
                         double charge, int nrf, double harmon[],
                         double voltages[], double epsilon) {
  // Set the initial option prices and volatility
  double y = VeffRFeVRadlosses(init_phi, U0, charge, nrf, harmon, voltages);
  double x = init_phi;

  // Newton Raphson
  while (fabs(y - target) > epsilon) {
    double d_x =
        EffectiveRFVoltageInElectronVoltPrime(x, charge, nrf, harmon, voltages);
    x += (target - y) / d_x;
    y = VeffRFeVRadlosses(x, U0, charge, nrf, harmon, voltages);
  }

  return x;
}

/*
================================================================================
================================================================================
METHOD TO GET THE EFFECTIVE RF VOLTAGE TIMES CHARGE WITH POTENTIAL WELL
DISTORTION. (units eV).

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF :
    Intrabeam Scattering Studies at the Cornell Electron-positron Storage Ring
Test Accelerator

================================================================================
  Arguments:
  ----------
    - double phi
        phase
    - double U0
        Energy loss per turn in eV
    - double charge
        particle charge
    - double nrf
        number of RF systems
    - doube[] harmon
        corresponding harmonic numbers
    - double[] voltages
        corresponding voltages
    - double L
    - double N
    - double sigs
        bunch length sigma s
    - double pc
        beam momentum
    - double epsilon
        threshold accuracy for the return value

  Returns:
  --------
    double
      Effective RF Voltage times charge in the presence of potential well
      distortion.
================================================================================
================================================================================
*/
double VeffRFeVPotentialWellDistortion(double phi, double U0, double charge,
                                       int nrf, double harmon[],
                                       double voltages[], double L, double N,
                                       double sigs, double pc) {

  return VeffRFeVRadlosses(phi, U0, charge, nrf, harmon, voltages) +
         charge * L * N * ec * clight * clight /
             (sqrt(2 * pi) * sigs * sigs * sigs * pc * 1e9);
}

/*
================================================================================
================================================================================
METHOD TO GET THE DERIVATIVE OF THE  EFFECTIVE RF VOLTAGE TIMES CHARGE WITH
POTENTIAL WELL DISTORTION WITH RESPECT OT THE PHASE. (units eV).

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF :
    Intrabeam Scattering Studies at the Cornell Electron-positron Storage Ring
Test Accelerator

================================================================================
  Arguments:
  ----------
    - double phi
        phase
    - double U0
        Energy loss per turn in eV
    - double charge
        particle charge
    - double nrf
        number of RF systems
    - doube[] harmon
        corresponding harmonic numbers
    - double[] voltages
        corresponding voltages
    - double L
    - double N
    - double sigs
        bunch length sigma s
    - double pc
        beam momentum

  Returns:
  --------
    double
      Effective RF Voltage times charge in the presence of potential well
      distortion derivative with respect to the phase.
================================================================================
================================================================================
*/
double VeffRFeVPotentialWellDistortionPrime(double phi, double U0,
                                            double charge, int nrf,
                                            double harmon[], double voltages[],
                                            double L, double N, double sigs,
                                            double pc) {
  return EffectiveRFVoltageInElectronVoltPrime(phi, charge, nrf, harmon,
                                               voltages) +
         charge * ec * L * N * clight * clight /
             (sqrt(2.0 * pi) * sigs * sigs * sigs * pc * 1e9);
}

/*
================================================================================
================================================================================
METHOD TO GET THE SYNCHRONUOUS PHASE IN THE PRESENCE OF POTENTIAL WELL
DISTORTION.

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF :
    https://www.quantstart.com/articles/Implied-Volatility-in-C-using-Template-Functions-and-Newton-Raphson/

================================================================================
  Arguments:
  ----------
    - double target
        target value for the test function (zero)
    - double init_phi
        initial phase to start searching from
    - double U0
        Energy loss per turn in eV
    - double charge
        particle charge
    - double nrf
        number of RF systems
    - doube[] harmon
        corresponding harmonic numbers
    - double[] voltages
        corresponding voltages
    - double epsilon
        threshold accuracy for the return value
    - double L
    - double N
    - double sigs
        bunch length sigma s
    - double pc
        beam momentum

  Returns:
  --------
    double
      synchronuous phase with potential well distortion.

================================================================================
================================================================================
*/
double SynchronuousPhaseWithPWD(double target, double init_phi, double U0,
                                double charge, int nrf, double harmon[],
                                double voltages[], double L, double N,
                                double sigs, double pc, double epsilon) {
  // Set the initial option prices and volatility
  double y = VeffRFeVPotentialWellDistortion(init_phi, U0, charge, nrf, harmon,
                                             voltages, L, N, sigs, pc);
  double x = init_phi;

  while (fabs(y - target) > epsilon) {
    double d_x = VeffRFeVPotentialWellDistortionPrime(
        x, U0, charge, nrf, harmon, voltages, L, N, sigs, pc);
    x += (target - y) / d_x;
    y = VeffRFeVPotentialWellDistortion(x, U0, charge, nrf, harmon, voltages, L,
                                        N, sigs, pc);
  }
  return x;
}
/*
================================================================================
================================================================================
METHOD TO GET THE SYNCHROTRON TUNE.

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

================================================================================
  Arguments:
  ----------
    - double omega0
        accelerator angular frequency
    - double U0
        Energy loss per turn in eV
    - double charge
        particle charge
    - double nrf
        number of RF systems
    - doube[] harmon
        corresponding harmonic numbers
    - double[] voltages
        corresponding voltages
    - double eta
        phase slip factor
    - double phis
        synchronuous phase
    - double pc
        beam momentum in GeV

  Returns:
  --------
    double
      Synchrotron Tune

================================================================================
================================================================================
*/
double SynchrotronTune(double omega0, double U0, double charge, int nrf,
                       double harmon[], double voltages[], double phis,
                       double eta, double pc) {
  return sqrt(harmon[0] * eta *
              fabs(EffectiveRFVoltageInElectronVoltPrime(phis, charge, nrf,
                                                         harmon, voltages)) /
              (2 * pi * pc * 1e9));
}

/*
================================================================================
================================================================================
METHOD TO GET THE SYNCHROTRON TUNE IN THE PRESENCE OF POTENTIAL WELL DISTORTION.

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

================================================================================
  Arguments:
  ----------
    - double omega0
        accelerator angular frequency
    - double U0
        Energy loss per turn in eV
    - double charge
        particle charge
    - double nrf
        number of RF systems
    - doube[] harmon
        corresponding harmonic numbers
    - double[] voltages
        corresponding voltages
    - double L
    - double N
    - double sigs
    - double eta
        phase slip factor
    - double phis
        synchronuous phase
    - double pc
        beam momentum in GeV

  Returns:
  --------
    double
      Synchrotron Tune

================================================================================
================================================================================
*/
double SynchrotronTunePWD(double omega0, double U0, double charge, int nrf,
                          double harmon[], double voltages[], double L,
                          double N, double sigs, double phis, double eta,
                          double pc) {
  return sqrt(harmon[0] * eta *
              fabs(VeffRFeVPotentialWellDistortionPrime(
                  phis, U0, charge, nrf, harmon, voltages, L, N, sigs, pc)) /
              (2 * pi * pc * 1e9));
}

/*
================================================================================
================================================================================
METHOD TO GET THE SYNCHROTRON TUNE IN THE PRESENCE OF POTENTIAL WELL DISTORTION.

================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

================================================================================
  Arguments:
  ----------
    - double omega0
        accelerator angular frequency
    - double U0
        Energy loss per turn in eV
    - double charge
        particle charge
    - double nrf
        number of RF systems
    - doube[] harmon
        corresponding harmonic numbers
    - double[] voltages
        corresponding voltages
    - double L
    - double N
    - double sigs
    - double eta
        phase slip factor
    - double phis
        synchronuous phase
    - double pc
        beam momentum in GeV

  Returns:
  --------
    double
      Synchrotron Tune

================================================================================
================================================================================
*/
double SigeFromRFAndSigs(double sigs, double U0, double charge, int nrf,
                         double harmon[], double voltages[], double gamma,
                         double gammatr, double pc, double circ, double phis,
                         bool printout) {

  double betar = BetaRelativisticFromGamma(gamma);

  double trev = circ / (betar * clight);
  double frev = 1.0 / trev;
  double neta = eta(gamma, gammatr);
  double omega0 = 2.0 * pi * frev;

  double nus = SynchrotronTune(omega0, U0, charge, nrf, harmon, voltages, phis,
                               neta, pc);
  double sige; // = nus * omega0 * sigs / (clight * eta);

  sige = 2.0 * pi * nus * sigs / (neta * circ);
  sige = omega0 * nus * sigs / (neta * clight);

  if (printout) {
    printf("Synchrotron Tune : %12.6e\n ", nus);
    printf("Synchrotron freq : %12.6e\n ", nus * omega0);
    printf("Sige             : %12.6e\n ", sige);
  }
  return sige;
}

/*
================================================================================
================================================================================
METHOD TO UPDATE THE MADX TWISS TABLE WITH COLUMNS (KEYS IN THE MAP) WITH
RECURRING VALUES IN THE IBS AND RADIATION DAMPING ROUTINES.

REMARK:
-------
  The code below is not optimized for speed but for readability.
  In principle this function is run once per simulation so the
  time gain of optimization would be minimal.
================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 06/08/2021 : initial version (Tom)

================================================================================
  Arguments:
  ----------
    - map<string, <double>>& : table
      Twiss table

  Returns:
  --------
    void
      updates the Twiss table (as map) with new keys

================================================================================
================================================================================
 */
void updateTwiss(map<string, vector<double>> &table) {
  // get length of table to reserve the vector sizes
  int size = table["L"].size();

  vector<double> rho(size), k(size), I1(size), I2(size), I3(size), I4x(size),
      I4y(size), I5x(size), I5y(size), gammax(size), gammay(size), hx(size),
      hy(size);

  complex<double> kc, k2c, klc;

  // calculate the new columns
  for (int i = 0; i < size; i++) {
    double angle = table["ANGLE"][i];
    double l = table["L"][i];
    double bx = table["BETX"][i];
    double by = table["BETY"][i];
    double ax = table["ALFX"][i];
    double ay = table["ALFY"][i];
    double dx = table["DX"][i];
    double dpx = table["DPX"][i];
    double dy = table["DY"][i];
    double dpy = table["DPY"][i];
    double k1l = table["K1L"][i];
    double k1sl = table["K1SL"][i];
    double e1 = angle / 2.0;
    double e2 = angle / 2.0;
    double rhoi2, rhoi3;

    // local bending radius
    rho[i] = (angle == 0.0) ? 0.0 : l / angle;
    rhoi2 = rho[i] * rho[i];
    rhoi3 = rhoi2 * rho[i];

    if (rho[i] != 0.0) {
      double rhoinv = 1.0 / rho[i];

      // effect of poleface rotation
      double alfx = ax - bx * tan(e1) * rhoinv;
      double dpxx = dpx + dx * tan(e1) * rhoinv;
      double gamx = (1.0 + alfx * alfx) / bx;

      // global gradient combining weak focusing and dipole gradient
      // k2 can be positive or negative and k can be real or imaginary
      k2c = rhoinv * rhoinv + 2.0 * k1l / l;
      kc = sqrt(k2c);
      klc = kc * l;

      // propagation of dispersion at exit
      double dx2 = real(dx * cos(klc) + dpxx * sin(klc) / kc +
                        rhoinv * (1.0 - cos(klc)) / (kc * kc));

      double dispaverage =
          real(dx * sin(klc) / klc + dpxx * (1.0 - cos(klc)) / (kc * klc) +
               rhoinv * (klc - sin(klc)) / (k2c * klc));

      /*
      double curlyhaverage = real(
          gamx * dx * dx + 2 * alfx * dx * dpxx + bx * dpxx * dpxx +
          2 * rhoinv * l *
              (-(gamx * dx + alfx * dpxx) * (klc - sin(klc)) /
                   (klc * klc * kc) +
               (alfx * dx + bx * dpxx) * (1.0 - cos(klc)) / (klc * klc)) +
          l * l * rhoinv * rhoinv *
              (gamx * (3.0 * klc - 4.0 * sin(klc) + sin(klc) * cos(klc)) /
                   (2.0 * k2c * klc * klc * klc) -
               alfx * (1.0 - cos(klc)) * (1.0 - cos(klc)) /
                   (kc * klc * klc * klc) +
               bx * (klc - cos(klc) * sin(klc)) / (2.0 * klc * klc * klc)));
               */

      k[i] = k1l / l;

      // 11/06/2021 - updated to match madx integrals
      I1[i] = (dx * rhoinv) * l;
      I2[i] = l / rhoi2;
      I3[i] = l / fabs(rhoi3);
      I4x[i] = dispaverage * rhoinv * (rhoinv * rhoinv + 2 * k[i]) * l -
               rhoinv * rhoinv * (dx * tan(e1) + dx2 * tan(e2)) +
               2.0 * rhoinv * k1sl * dy;

      hx[i] = bx * dpxx * dpxx + 2.0 * alfx * dx * dpxx + gamx * dx * dx;
      // hx[i] = curlyhaverage;
      I5x[i] = hx[i] * l / rhoi3;
      // I5x[i] = curlyhaverage * l / rhoi3;
    } else {
      I1[i] = 0.0;
      I2[i] = 0.0;
      I3[i] = 0.0;
      I4x[i] = 0.0;
    }

    I4y[i] = 0.0;

    // Courant-Snyder gamma
    gammay[i] = (1.0 + ay * ay) / by;

    // curly H
    hy[i] = by * dpy * dpy + 2.0 * ay * dy * dpy + gammay[i] * dy * dy;

    I5y[i] = (rho[i] == 0) ? 0.0 : hy[i] * l / rhoi3;
  }

  // extend the map table
  table["rho"] = rho;
  table["k"] = k;

  table["gammax"] = gammax;
  table["gammay"] = gammay;
  table["hx"] = hx;
  table["hy"] = hy;

  table["I1"] = I1;
  table["I2"] = I2;
  table["I3"] = I3;
  table["I4x"] = I4x;
  table["I4y"] = I4y;
  table["I5x"] = I5x;
  table["I5y"] = I5y;
}

/*
================================================================================
================================================================================
METHOD TO PRINT A COLUMN FROM THE UPDATED TWISS TABLE (DEBUGGING)
================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 06/08/2021 : initial version (Tom)

================================================================================
  Arguments:
  ----------
    - string : key
      key in the map - column name  of the Twiss table
    - map<string, <double>> : table
      Twiss or extended Twiss table

  Returns:
  --------
    void

================================================================================
================================================================================
*/
void printTwissMap(string key, map<string, vector<double>> &table) {
  vector<double> toprint = table[key];

  for (auto elem : toprint) {
    cout << elem << endl;
  }
}

/*
================================================================================
================================================================================
METHOD CONVERT DE/E TO DP/P IN ABSENCE OF ELECTRIC POTENTIAL
================================================================================
  AUTHORS:
    - TOM MERTENS
    - MALTE TITZE

  HISTORY:
    - 08/07/2021 : initial version (Tom and Malte)

================================================================================
  Arguments:
  ----------
    - double : de/e
      dE/E
    - double : beta0
      relativistic beta

  Returns:
  --------
    double:
      DP/P

================================================================================
================================================================================
*/

double dee_to_dpp(double dee, double beta0) {
  return sqrt(((dee + 1.0) * (dee + 1.0) - 1.0) / (beta0 * beta0) + 1.0) - 1.0;
}

/*
================================================================================
================================================================================
METHOD CONVERT DP/P TO DE/E IN ABSENCE OF ELECTRIC POTENTIAL
================================================================================
  AUTHORS:
    - TOM MERTENS
    - MALTE TITZE

  HISTORY:
    - 08/07/2021 : initial version (Tom and Malte)

================================================================================
  Arguments:
  ----------
    - double : dP/P
      dP/P
    - double : beta0
      relativistic beta

  Returns:
  --------
    double:
      DE/E

================================================================================
================================================================================
*/

double dpp_to_dee(double dpp, double beta0) {
  return sqrt((dpp + 1.0) * (dpp + 1.0) - (beta0 * beta0)) - 1.0;
}