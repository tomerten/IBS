#include "../include/ibs_bits/CoulombLogFunctions.hpp"
#include "../include/ibs_bits/Integrators.hpp"
#include "../include/ibs_bits/Models.hpp"
#include "../include/ibs_bits/NumericFunctions.hpp"
#include "../include/ibs_bits/RadiationDamping.hpp"
#include "../include/ibs_bits/twiss.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

// text color defs for printf
void red() { printf("\033[1;31m"); }
void yellow() { printf("\033[1;33m"); }
void green() { printf("\033[1;32m"); }
void blue() { printf("\033[1;34m"); }
void cyan() { printf("\033[1;36m"); }
void reset() { printf("\033[0m"); }

// method to use printf to print a formatted line
void printline(string key, double value, string units) {
  printf("%-20s : %20.6e (%s)\n", key.c_str(), value, units.c_str());
}

// write data to file
void WriteToFile(string filename, vector<double> &t, vector<double> &ex,
                 vector<double> &ey, vector<double> &sigs) {

  ofstream csvfile(filename);

  if (csvfile.is_open()) {

    // min number of rows
    int num_of_rows = min({t.size(), ex.size(), ey.size(), sigs.size()});

    // column headers
    csvfile << "t,ex,ey,sigs" << endl;

    // write data line by line
    for (int i = 0; i < num_of_rows; i++) {
      csvfile << t[i] << "," << ex[i] << "," << ey[i] << "," << sigs[i] << endl;
    }
  } else {
    cout << "File could not be opened";
  }
  csvfile.close();
}

/*
================================================================================
================================================================================
METHOD TO PERORM ODE SIMULATIONS TO GET ESTIMATES FOR EQUILIBRIUM
EX, EY AND SIGS TAKING RADIATION DAMPING, QUANTUM EXCITATION AND IBS
INTO ACCOUNT.

AS EQUILIBRIA DO NOT NECESSARILY EXIST A SAFETY MAX STEPS OF 10000 IS APPLIED.
================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : INITIAL VERSION - ISSUES WITH COUPLING AND HOR EMIT
    - 12/07/2021 : UPDATE TO FIX PREVIOUS ISSUES
                    * ADDED TWO METHODS TO PERFORM SIMULATION
                      + USING RELAXATION (EQ. 47 IN REF)
                      + USING DERIVATIVES (BMAD REF)

  REFS:
    - BMAD SOURCE CODE - based on ibs_mod.f90 and ibs_ring.f90
    - PHYS REV ST ACCEL BEAMS 081001 (2005) (EQ 47)

  NOTES:
    - coupling uses eq. 47 from above where the vertical ey0 is calculated
      as : max(ex0 * coupling, ey0_rad) with (0.0 < coupling <=1.0)
    - coupling is very rudimentary proceed with caution when using these values.
    - EQ 47 is only valid in betatron dominated coupling, if vertical dispersion
      is large the approximation is not valid.

================================================================================
  Arguments:
  ----------
    - map<string, double> &twiss
        twiss header
    - map<string, vector<double>> &twissdata
        twiss data as map of double vectors
    - int nrf
        number of rf systems
    - double[] harmon
        array of the harmonic numbers of the rf systems
    - double[] voltages
        array of the voltages of the rf systems
    - vector<double> &t
        vector of timestamps - as input : single initial value in the vector
    - vector<double> &ex
        vector of horizontal emittance - as input : single initial value in the
        vector
    - vector<double> &ey
        vector of vertical emittance - as input : single initial value in the
        vector
    - vector<double> &sigs
        vector of bunch lengths sigma s - as input : single initial value in the
        vector
    - vector<double> &sige
        vector of energy spreads sigma E - as input : single initial value in
        the vector
    - int model
        integer to select the IBS models
    - double pnumber
        number of particles in the bunch
    - int couplingpercentage
        horizontal betatron coupling
    - double threshold
        while loop cutoff for relative changes in the values
    - string method
        method to use : rlx or der

  Returns:
  --------
    - vector<double> &t
        vector of timestamps
    - vector<double> &ex
        vector of horizontal emittance
        vector
    - vector<double> &ey
        vector of vertical emittance
        vector
    - vector<double> &sigs
        vector of bunch lengths sigma s
    - vector<double> &sige
        vector of energy spreads sigma E
================================================================================
================================================================================
*/

void ODE(map<string, double> &twiss, map<string, vector<double>> &twissdata,
         int nrf, double harmon[], double voltages[], vector<double> &t,
         vector<double> &ex, vector<double> &ey, vector<double> &sigs,
         vector<double> sige, int model, double pnumber, int couplingpercentage,
         double threshold, string method) {

  // safetey max steps
  int MaxSteps = 10000;

  // sanitize limit settings
  if (!(method == "rlx" || method == "der")) {
    method = "der";

    red();
    printf("%-20s : (%s)\n", "Warning method set to ", method.c_str());
    reset();
  }

  if (threshold > 1.0 || threshold < 1.0e-6) {
    threshold = 1e-4;
  }

  if (couplingpercentage > 100 || couplingpercentage < 0) {
    couplingpercentage = 0;
  }

  // rescale coupling percentage
  double coupling = (double)couplingpercentage / 100.0;

  // Radiation integrals
  double twiss_rad[6];
  double *radint;

  double gamma = twiss["GAMMA"];
  double pc = twiss["PC"];
  double gammatr = twiss["GAMMATR"];
  double mass = twiss["MASS"];
  double charge = twiss["CHARGE"];
  double q1 = twiss["Q1"];
  double len = twiss["LENGTH"];

  double aatom = emass / pmass;
  double betar = BetaRelativisticFromGamma(gamma);
  double r0 = ParticleRadius(1, aatom);
  double trev = len / (betar * clight);
  double frev = 1.0 / trev;
  double omega = 2.0 * pi * frev;
  double neta = eta(gamma, gammatr);
  double epsilon = 1.0e-6;

  // get radiation integrals
  radint = RadiationDampingLattice(twissdata);

  // Longitudinal Parameters
  double U0 = RadiationLossesPerTurn(twiss, radint[1], aatom);
  double phis =
      SynchronuousPhase(0.0, 173, U0, charge, nrf, harmon, voltages, epsilon);
  double qs =
      SynchrotronTune(omega, U0, charge, nrf, harmon, voltages, phis, neta, pc);
  double omegas = qs * omega;

  // equilibria
  double *equi =
      RadiationDampingLifeTimesAndEquilibriumEmittancesWithPartitionNumbers(
          twiss, radint, aatom, qs);

  double tauradx, taurady, taurads, sigeoe2;

  tauradx = equi[0];
  taurady = equi[1];
  taurads = equi[2];
  sigeoe2 = equi[5];

  double ey0_coupled = max(coupling * equi[3], equi[4]);
  double sige0 = sigefromsigs(omega, equi[6], qs, gamma, gammatr);

  cyan();
  printf("Radiation Damping Times\n");
  printf("=======================\n");
  printline("Tau_rad_x", tauradx, "s");
  printline("Tau_rad_y", taurady, "s");
  printline("Tau_rad_s", tauradx, "s");

  blue();
  printf("\nLongitudinal Parameters\n");
  printf("=======================\n");
  printline("Synchrotron Tune", qs, "");
  printline("Synchrotron Freq", omegas, "Hz");
  printline("SigEOE2", sigeoe2, "");
  printline("SigEOE ", sqrt(sigeoe2), "");
  printline("eta", eta(gamma, gammatr), "");
  printline("Sigs", sigs[0], "");
  printline("Sigs_inf ", equi[6], "");
  printline("SigE0 ", sige0, "");

  sige0 = SigeFromRFAndSigs(equi[6], U0, charge, nrf, harmon, voltages, gamma,
                            gammatr, pc, len, phis, false);

  // check value
  cyan();
  printf("%-20s : %20.6e (%s)\n", "Sige0 - check", sige0, "");
  reset();

  // write first sige and sige2
  vector<double> sige2;
  sige.push_back(sige0);
  sige2.push_back(sige[0] * sige[0]);

  // loop variable
  int i = 0;

  // ibs growth rates
  double *ibs;
  double aes, aex, aey;

  // temp emit vectors to allow for correcting with quantum excitation
  vector<double> extemp, eytemp, sige2temp;
  extemp.push_back(ex[0]);
  eytemp.push_back(ey[0]);
  sige2temp.push_back(sige2[0]);

  // initial ibs growth rates
  switch (model) {
  case 1:
    ibs = PiwinskiSmooth(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss, r0);
    break;
  case 2:
    ibs = PiwinskiLattice(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                          twissdata, r0);
    break;
  case 3:
    ibs = PiwinskiLatticeModified(pnumber, ex[0], ey[0], sigs[0], sige[0],
                                  twiss, twissdata, r0);
    break;
  case 4:
    ibs = Nagaitsev(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss, twissdata,
                    r0);
    break;
  case 5:
    ibs = Nagaitsevtailcut(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                           twissdata, r0, aatom);
    break;
  case 6:
    ibs = ibsmadx(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss, twissdata, r0,
                  false);
    break;
  case 7:
    ibs = ibsmadxtailcut(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                         twissdata, r0, aatom);
    break;
  case 8:
    ibs = BjorkenMtingwa2(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                          twissdata, r0);
    break;
  case 9:
    ibs = BjorkenMtingwa(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                         twissdata, r0);
    break;
  case 10:
    ibs = BjorkenMtingwatailcut(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                                twissdata, r0, aatom);
    break;
  case 11:
    ibs = ConteMartini(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                       twissdata, r0);
    break;
  case 12:
    ibs = ConteMartinitailcut(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                              twissdata, r0, aatom);
    break;
  case 13:
    ibs =
        MadxIBS(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss, twissdata, r0);
    break;
  }

  // get max tau limited to max 1.0 sec
  double taum = max(tauradx, taurady);
  taum = max(taum, taurads);
  taum = max(taum, 1.0 / ibs[0]);
  taum = max(taum, 1.0 / ibs[1]);
  taum = max(taum, 1.0 / ibs[2]);
  taum = min(taum, 1.0);

  // get min to auto derive stepsize
  double ddt = min(tauradx, taurady);
  ddt = min(ddt, taurads);
  ddt = min(ddt, 1.0 / ibs[0]);
  ddt = min(ddt, 1.0 / ibs[1]);
  ddt = min(ddt, 1.0 / ibs[2]);

  // reduce the step still a bit
  // ddt /= 2.0;

  // define max numer of steps
  int ms = (int)(10 * taum / ddt);
  ms = min(ms, MaxSteps);

  // print initial IBS summary
  printouts(ibs);

  // print run parameters
  red();
  printf("\nMax tau : %12.6e\n", taum);
  printf("dt      : %12.6e\n", ddt);
  printf("Max step: %i\n\n", ms);
  printf("Coupling: %12.6f\n\n", coupling);
  reset();

  /*
  ================================================================================
  MAIN LOOP
  ================================================================================
  */
  // progressbar
  int barWidth = 70;
  do {
    std::cout << "[";
    int progress = (double)i / ms * barWidth;
    for (int j = 0; j < barWidth; ++j) {
      if (j < progress)
        std::cout << "=";
      else if (j == progress)
        std::cout << ">";
      else
        std::cout << " ";
    }
    std::cout << "]" << int((double)i / ms * 100) << " %\r";
    std::cout.flush();
    // update timestep
    ddt = min(tauradx, taurady);
    ddt = min(ddt, taurads);
    ddt = min(ddt, 1.0 / ibs[0]);
    ddt = min(ddt, 1.0 / ibs[1]);
    ddt = min(ddt, 1.0 / ibs[2]);
    // ddt /= 2.0;

    // ibs growth rates update
    switch (model) {
    case 1:
      ibs = PiwinskiSmooth(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 2:
      ibs = PiwinskiLattice(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                            twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 3:
      ibs = PiwinskiLatticeModified(pnumber, ex[i], ey[i], sigs[i], sige[i],
                                    twiss, twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 4:
      ibs = Nagaitsev(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss, twissdata,
                      r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 5:
      ibs = Nagaitsevtailcut(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                             twissdata, r0, aatom);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 6:
      ibs = ibsmadx(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss, twissdata,
                    r0, false);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 7:
      ibs = ibsmadxtailcut(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                           twissdata, r0, aatom);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 8:
      ibs = BjorkenMtingwa2(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                            twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 9:
      ibs = BjorkenMtingwa(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                           twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 10:
      ibs = BjorkenMtingwatailcut(pnumber, ex[i], ey[i], sigs[i], sige[i],
                                  twiss, twissdata, r0, aatom);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 11:
      ibs = ConteMartini(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                         twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 12:
      ibs = ConteMartinitailcut(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                                twissdata, r0, aatom);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 13:
      ibs = MadxIBS(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss, twissdata,
                    r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    }

    i++;

    if (method == "rlx") {
      ddt *= 2.;
      double ratio_x = tauradx * aex;
      double ratio_y = taurady * aey;
      double ratio_s = taurads * aes;

      double xfactor = 1.0 / (1.0 - ratio_x);
      double yfactor = 1.0 / (1.0 - ratio_y);
      double sfactor = 1.0 / (1.0 - ratio_s);

      t.push_back(t[i - 1] + ddt);
      ex.push_back(ex[i - 1] + ddt * (xfactor * equi[3] - ex[i - 1]));
      ey.push_back(ey[i - 1] +
                   ddt * (((1.0 - coupling) * yfactor + coupling * xfactor) *
                              ey0_coupled -
                          ey[i - 1]));
      sige.push_back(sige[i - 1] +
                     ddt * (sfactor * sqrt(equi[5]) - sige[i - 1]));
      sige2.push_back(sige[i] * sige[i]);
      sigs.push_back(sigsfromsige(sige[i], gamma, gammatr, omegas));
    } else {
      double dxdt =
          -(ex[i - 1] - equi[3]) * 2. / tauradx + ex[i - 1] * 2.0 * aex;
      double dydt =
          -(ey[i - 1] - ey0_coupled) * 2. / taurady + ey[i - 1] * 2.0 * aey;
      double dedt =
          -(sige[i - 1] - sqrt(equi[5])) / taurads + sige[i - 1] * aes;

      t.push_back(t[i - 1] + ddt);
      ex.push_back(ex[i - 1] + ddt * dxdt);
      ey.push_back(ey[i - 1] + ddt * dydt);
      sige.push_back(sige[i - 1] + ddt * dedt);
      sige2.push_back(sige[i] * sige[i]);
      sigs.push_back(sigsfromsige(sige[i], gamma, gammatr, omegas));
    }

    /*
      double tmpx = aex;
      aex =
          (2.0 * aex * aey) / ((2.0 - coupling) * aey +
                               coupling * aex); // full coupling -> (tx + ty) /
  2 aey = (2.0 * aex * aey) / ((2.0 - coupling) * tmpx + coupling * aey); //
  full coupling -> (tx + ty) / 2

      // update loop variable
      i++;

      // add ddt to previous time stamp
      t.push_back(t[i - 1] + ddt);

      // horizontal emittance
      ex.push_back(extemp[i - 1] * exp(2 * ddt * (-1 / tauradx + aex)) +
                   equi[3] * (1 - exp(-2 * ddt * i / tauradx)));
      extemp.push_back(extemp[i - 1] * exp(2 * ddt * (-1 / tauradx + aex)));

      // EQ. 47 in PHYS REV ST ACCEL BEAMS 081001 (2005)
      if (couplingpercentage > 0) {
        ey.push_back(((1 - coupling) * (1.0 / aey) / ((1.0 / aey) - taurady) +
                      coupling * (1.0 / aex) / ((1.0 / aex) - tauradx)) *
                     equi[4]);
      } else {
        // vertical emittance
        ey.push_back(eytemp[i - 1] * exp(2 * ddt * (-1 / taurady + aey)) +
                     equi[4] * (1 - exp(-2 * ddt * i / taurady)));
        eytemp.push_back(eytemp[i - 1] * exp(2 * ddt * (-1 / taurady + aey)));
      }

      // bunch length
      sige2.push_back(sige2temp[i - 1] * exp(2 * ddt * (-1 / taurads + aes)) +
                      equi[5] * (1 - exp(-2 * i * ddt / taurads)));
      sige2temp.push_back(sige2temp[i - 1] * exp(2 * ddt * (-1 / taurads +
  aes)));

      sige.push_back(sqrt(sige2[i]));
      sigs.push_back(sigsfromsige(sige[i], gamma, gammatr, omegas));

  */
    // while condition
  } while (i < ms && (fabs((ex[i] - ex[i - 1]) / ex[i - 1]) > threshold ||
                      fabs((ey[i] - ey[i - 1]) / ey[i - 1]) > threshold ||
                      fabs((sigs[i] - sigs[i - 1]) / sigs[i - 1]) > threshold));
  // end progressbar
  std::cout << std::endl;

  // print final values
  blue();
  printf("%-20s : %12.6e\n", "Final ex", ex[ex.size() - 1]);
  printf("%-20s : %12.6e\n", "Final ey", ey[ey.size() - 1]);
  printf("%-20s : %12.6e\n", "Final sigs", sigs[sigs.size() - 1]);

  printf("%-20s : %12.6e\n", "Final tau_ibs_x", 1.0 / ibs[1]);
  printf("%-20s : %12.6e\n", "Final tau_ibs_y", 1.0 / ibs[2]);
  printf("%-20s : %12.6e\n", "Final tau_ibs_s", 1.0 / ibs[0]);
  reset();
}

void ODE(map<string, double> &twiss, map<string, vector<double>> &twissdata,
         int nrf, double harmon[], double voltages[], vector<double> &t,
         vector<double> &ex, vector<double> &ey, vector<double> &sigs,
         vector<double> sige, int model, double pnumber, int nsteps,
         double stepsize, int couplingpercentage) {
  // sanitize limit settings
  if (couplingpercentage > 100 || couplingpercentage < 0) {
    couplingpercentage = 0;
  }

  double coupling = (double)couplingpercentage / 100.0;

  // Radiation integrals
  double gamma = twiss["GAMMA"];
  double pc = twiss["PC"];
  double gammatr = twiss["GAMMATR"];
  double mass = twiss["MASS"];
  double charge = twiss["CHARGE"];
  double q1 = twiss["Q1"];
  double len = twiss["LENGTH"];
  double twiss_rad[6];

  double aatom = emass / pmass;
  double betar = BetaRelativisticFromGamma(gamma);
  double r0 = ParticleRadius(1, aatom);
  double trev = len / (betar * clight);
  double frev = 1.0 / trev;
  double omega = 2.0 * pi * frev;
  double neta = eta(gamma, gammatr);
  double epsilon = 1.0e-6;

  double *radint;
  radint = RadiationDampingLattice(twissdata);

  // Longitudinal Parameters
  double U0 = RadiationLossesPerTurn(twiss, radint[1], aatom);
  double phis =
      SynchronuousPhase(0.0, 173, U0, charge, nrf, harmon, voltages, epsilon);
  double qs =
      SynchrotronTune(omega, U0, charge, nrf, harmon, voltages, phis, neta, pc);
  double omegas = qs * omega;

  // equilibria
  double *equi =
      RadiationDampingLifeTimesAndEquilibriumEmittancesWithPartitionNumbers(
          twiss, radint, aatom, qs);

  double tauradx, taurady, taurads;
  tauradx = equi[0];
  taurady = equi[1];
  taurads = equi[2];
  double sigeoe2 = equi[5];

  cyan();
  printf("Radiation Damping Times\n");
  printf("=======================\n");
  printf("%-30s %20.6e (%s)\n", "Tx :", tauradx, "");
  printf("%-30s %20.6e (%s)\n", "Ty :", taurady, "");
  printf("%-30s %20.6e (%s)\n", "Ts :", taurads, "");

  blue();
  printf("\nLongitudinal Parameters\n");
  printf("=======================\n");
  printf("%-20s : %20.6e (%s)\n", "qs", qs, "");
  printf("%-20s : %20.6e (%s)\n", "synch freq", omegas, "");
  printf("%-20s : %20.6e (%s)\n", "SigEOE2", sigeoe2, "");
  printf("%-20s : %20.6e (%s)\n", "SigEOE", sqrt(sigeoe2), "");
  printf("%-20s : %20.6e (%s)\n", "eta", eta(gamma, gammatr), "");
  printf("%-20s : %20.6e (%s)\n", "Sigs", sigs[0], "");
  printf("%-20s : %20.6e (%s)\n", "Sigsinf", equi[6], "");

  double sige0 = sigefromsigs(omega, equi[6], qs, gamma, gammatr);
  printf("%-20s : %20.6e (%s)\n", "Sige0", sige0, "");
  sige0 = SigeFromRFAndSigs(equi[6], U0, charge, nrf, harmon, voltages, gamma,
                            gammatr, pc, len, phis, false);
  printf("%-20s : %20.6e (%s)\n", "Sige0 - check", sige0, "");
  reset();

  sige0 = SigeFromRFAndSigs(sigs[0], U0, charge, nrf, harmon, voltages, gamma,
                            gammatr, pc, len, phis, false);

  // write first sige and sige2
  vector<double> sige2;
  sige.push_back(sige0);
  sige2.push_back(sige[0] * sige[0]);

  // loop variable
  int i = 0;
  // ibs growth rates
  double *ibs;
  double aes, aex, aey;

  // temp emit vectors to allow for correcting with quantum excitation
  vector<double> extemp, eytemp, sige2temp;
  extemp.push_back(ex[0]);
  eytemp.push_back(ey[0]);
  sige2temp.push_back(sige2[0]);

  // initial ibs growth rates
  switch (model) {
  case 1:
    ibs = PiwinskiSmooth(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss, r0);
    break;
  case 2:
    ibs = PiwinskiLattice(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                          twissdata, r0);
    break;
  case 3:
    ibs = PiwinskiLatticeModified(pnumber, ex[0], ey[0], sigs[0], sige[0],
                                  twiss, twissdata, r0);
    break;
  case 4:
    ibs = Nagaitsev(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss, twissdata,
                    r0);
    break;
  case 5:
    ibs = Nagaitsevtailcut(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                           twissdata, r0, aatom);
    break;
  case 6:
    ibs = ibsmadx(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss, twissdata, r0,
                  false);
    break;
  case 7:
    ibs = ibsmadxtailcut(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                         twissdata, r0, aatom);
    break;
  case 8:
    ibs = BjorkenMtingwa2(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                          twissdata, r0);
    break;
  case 9:
    ibs = BjorkenMtingwa(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                         twissdata, r0);
    break;
  case 10:
    ibs = BjorkenMtingwatailcut(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                                twissdata, r0, aatom);
    break;
  case 11:
    ibs = ConteMartini(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                       twissdata, r0);
    break;
  case 12:
    ibs = ConteMartinitailcut(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                              twissdata, r0, aatom);
    break;
  case 13:
    ibs =
        MadxIBS(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss, twissdata, r0);
    break;
  }

  printouts(ibs);

  /*
  ================================================================================
  MAIN LOOP
  ================================================================================
  */
  do {
    // ibs growth rates update
    switch (model) {
    case 1:
      ibs = PiwinskiSmooth(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 2:
      ibs = PiwinskiLattice(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                            twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 3:
      ibs = PiwinskiLatticeModified(pnumber, ex[i], ey[i], sigs[i], sige[i],
                                    twiss, twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 4:
      ibs = Nagaitsev(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss, twissdata,
                      r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 5:
      ibs = Nagaitsevtailcut(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                             twissdata, r0, aatom);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 6:
      ibs = ibsmadx(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss, twissdata,
                    r0, false);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 7:
      ibs = ibsmadxtailcut(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                           twissdata, r0, aatom);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 8:
      ibs = BjorkenMtingwa2(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                            twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 9:
      ibs = BjorkenMtingwa(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                           twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 10:
      ibs = BjorkenMtingwatailcut(pnumber, ex[i], ey[i], sigs[i], sige[i],
                                  twiss, twissdata, r0, aatom);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 11:
      ibs = ConteMartini(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                         twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 12:
      ibs = ConteMartinitailcut(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                                twissdata, r0, aatom);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 13:
      ibs = MadxIBS(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss, twissdata,
                    r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    }

    // update transverse growth rates with coupling
    double tmpx = aex;
    aex =
        (2.0 * aex * aey) / ((2.0 - coupling) * aey +
                             coupling * aex); // full coupling -> (tx + ty) / 2
    aey =
        (2.0 * aex * aey) / ((2.0 - coupling) * tmpx +
                             coupling * aey); // full coupling -> (tx + ty) / 2

    i++;

    t.push_back(t[i - 1] + stepsize);
    ex.push_back(extemp[i - 1] * exp(2 * stepsize * (-1 / tauradx + aex)) +
                 equi[3] * (1 - exp(-2 * stepsize * i / tauradx)));
    extemp.push_back(extemp[i - 1] * exp(2 * stepsize * (-1 / tauradx + aex)));

    ey.push_back(eytemp[i - 1] * exp(2 * stepsize * (-1 / taurady + aey)) +
                 equi[4] * (1 - exp(-2 * stepsize * i / taurady)));
    eytemp.push_back(eytemp[i - 1] * exp(2 * stepsize * (-1 / taurady + aey)));

    sige2.push_back(sige2temp[i - 1] *
                        exp(2 * stepsize * (-1 / taurads + aes)) +
                    equi[5] * (1 - exp(-2 * i * stepsize / taurads)));
    sige2temp.push_back(sige2temp[i - 1] *
                        exp(2 * stepsize * (-1 / taurads + aes)));

    sige.push_back(sqrt(sige2[i]));
    sigs.push_back(sigsfromsige(sige[i], gamma, gammatr, omegas));

    /*
    // DEBUG
    printf("--------------------\n");
    printf("step : %10i\n", i);
    cyan();
    printf("ex   : %12.6e\n", ex[i]);
    printf("ey   : %12.6e\n", ey[i]);
    printf("sigs : %12.6e\n", sigs[i]);
    green();
    printf("exdiff abs  : %12.6e\n", (ex[i] - ex[i - 1]));
    printf("eydiff abs  : %12.6e\n", (ey[i] - ey[i - 1]));
    printf("sigsdiff abs: %12.6e\n", (sigs[i] - sigs[i - 1]));
    yellow();
    printf("exdiff   : %12.6e\n", fabs((ex[i] - ex[i - 1]) / ex[i - 1]));
    printf("eydiff   : %12.6e\n", fabs((ey[i] - ey[i - 1]) / ey[i - 1]));
    printf("sigsdiff : %12.6e\n", fabs((sigs[i] - sigs[i - 1]) / sigs[i -
    1])); reset();
    */
  } while (i < nsteps);
  blue();
  printf("%-20s : %12.6e\n", "Final ex", ex[ex.size() - 1]);
  printf("%-20s : %12.6e\n", "Final ey", ey[ey.size() - 1]);
  printf("%-20s : %12.6e\n", "Final sigs", sigs[sigs.size() - 1]);

  printf("%-20s : %12.6e\n", "Final tau_ibs_x", 1.0 / ibs[1]);
  printf("%-20s : %12.6e\n", "Final tau_ibs_y", 1.0 / ibs[2]);
  printf("%-20s : %12.6e\n", "Final tau_ibs_s", 1.0 / ibs[0]);
  reset();
}

void ODE_BMAD(map<string, double> &twiss,
              map<string, vector<double>> &twissdata, int nrf, double harmon[],
              double voltages[], vector<double> &t, vector<double> &ex,
              vector<double> &ey, vector<double> &sigs, vector<double> sige,
              int model, double pnumber, int nsteps, double stepsize,
              int couplingpercentage) {
  // sanitize limit settings
  if (couplingpercentage > 100 || couplingpercentage < 0) {
    couplingpercentage = 0;
  }

  double coupling = (double)couplingpercentage / 100.0;

  // Radiation integrals
  double gamma = twiss["GAMMA"];
  double pc = twiss["PC"];
  double gammatr = twiss["GAMMATR"];
  double mass = twiss["MASS"];
  double charge = twiss["CHARGE"];
  double q1 = twiss["Q1"];
  double len = twiss["LENGTH"];
  double twiss_rad[6];

  double aatom = emass / pmass;
  double betar = BetaRelativisticFromGamma(gamma);
  double r0 = ParticleRadius(1, aatom);
  double trev = len / (betar * clight);
  double frev = 1.0 / trev;
  double omega = 2.0 * pi * frev;
  double neta = eta(gamma, gammatr);
  double epsilon = 1.0e-6;

  double *radint;
  radint = RadiationDampingLattice(twissdata);

  // Longitudinal Parameters
  double U0 = RadiationLossesPerTurn(twiss, radint[1], aatom);
  double phis =
      SynchronuousPhase(0.0, 173, U0, charge, nrf, harmon, voltages, epsilon);
  double qs =
      SynchrotronTune(omega, U0, charge, nrf, harmon, voltages, phis, neta, pc);
  double omegas = qs * omega;

  // equilibria
  double *equi =
      RadiationDampingLifeTimesAndEquilibriumEmittancesWithPartitionNumbers(
          twiss, radint, aatom, qs);

  double tauradx, taurady, taurads;
  tauradx = equi[0];
  taurady = equi[1];
  taurads = equi[2];
  double sigeoe2 = equi[5];

  cyan();
  printf("Radiation Damping Times\n");
  printf("=======================\n");
  printf("%-30s %20.6e (%s)\n", "Tx :", tauradx, "");
  printf("%-30s %20.6e (%s)\n", "Ty :", taurady, "");
  printf("%-30s %20.6e (%s)\n", "Ts :", taurads, "");

  blue();
  printf("\nLongitudinal Parameters\n");
  printf("=======================\n");
  printf("%-20s : %20.6e (%s)\n", "qs", qs, "");
  printf("%-20s : %20.6e (%s)\n", "synch freq", omegas, "");
  printf("%-20s : %20.6e (%s)\n", "SigEOE2", sigeoe2, "");
  printf("%-20s : %20.6e (%s)\n", "SigEOE", sqrt(sigeoe2), "");
  printf("%-20s : %20.6e (%s)\n", "eta", eta(gamma, gammatr), "");
  printf("%-20s : %20.6e (%s)\n", "Sigs", sigs[0], "");
  printf("%-20s : %20.6e (%s)\n", "Sigsinf", equi[6], "");

  double sige0 = sigefromsigs(omega, equi[6], qs, gamma, gammatr);
  printf("%-20s : %20.6e (%s)\n", "Sige0", sige0, "");
  sige0 = SigeFromRFAndSigs(equi[6], U0, charge, nrf, harmon, voltages, gamma,
                            gammatr, pc, len, phis, false);
  printf("%-20s : %20.6e (%s)\n", "Sige0 - check", sige0, "");
  reset();

  sige0 = SigeFromRFAndSigs(sigs[0], U0, charge, nrf, harmon, voltages, gamma,
                            gammatr, pc, len, phis, false);

  // write first sige and sige2
  vector<double> sige2;
  sige.push_back(sige0);
  sige2.push_back(sige[0] * sige[0]);

  // loop variable
  int i = 0;
  // ibs growth rates
  double *ibs;
  double aes, aex, aey;

  // temp emit vectors to allow for correcting with quantum excitation
  vector<double> extemp, eytemp, sige2temp;
  extemp.push_back(ex[0]);
  eytemp.push_back(ey[0]);
  sige2temp.push_back(sige2[0]);

  // initial ibs growth rates
  switch (model) {
  case 1:
    ibs = PiwinskiSmooth(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss, r0);
    break;
  case 2:
    ibs = PiwinskiLattice(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                          twissdata, r0);
    break;
  case 3:
    ibs = PiwinskiLatticeModified(pnumber, ex[0], ey[0], sigs[0], sige[0],
                                  twiss, twissdata, r0);
    break;
  case 4:
    ibs = Nagaitsev(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss, twissdata,
                    r0);
    break;
  case 5:
    ibs = Nagaitsevtailcut(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                           twissdata, r0, aatom);
    break;
  case 6:
    ibs = ibsmadx(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss, twissdata, r0,
                  false);
    break;
  case 7:
    ibs = ibsmadxtailcut(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                         twissdata, r0, aatom);
    break;
  case 8:
    ibs = BjorkenMtingwa2(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                          twissdata, r0);
    break;
  case 9:
    ibs = BjorkenMtingwa(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                         twissdata, r0);
    break;
  case 10:
    ibs = BjorkenMtingwatailcut(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                                twissdata, r0, aatom);
    break;
  case 11:
    ibs = ConteMartini(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                       twissdata, r0);
    break;
  case 12:
    ibs = ConteMartinitailcut(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss,
                              twissdata, r0, aatom);
    break;
  case 13:
    ibs =
        MadxIBS(pnumber, ex[0], ey[0], sigs[0], sige[0], twiss, twissdata, r0);
    break;
  }

  printouts(ibs);

  /*
  ================================================================================
  MAIN LOOP
  ================================================================================
  */
  do {
    // ibs growth rates update
    switch (model) {
    case 1:
      ibs = PiwinskiSmooth(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 2:
      ibs = PiwinskiLattice(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                            twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 3:
      ibs = PiwinskiLatticeModified(pnumber, ex[i], ey[i], sigs[i], sige[i],
                                    twiss, twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 4:
      ibs = Nagaitsev(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss, twissdata,
                      r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 5:
      ibs = Nagaitsevtailcut(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                             twissdata, r0, aatom);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 6:
      ibs = ibsmadx(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss, twissdata,
                    r0, false);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 7:
      ibs = ibsmadxtailcut(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                           twissdata, r0, aatom);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 8:
      ibs = BjorkenMtingwa2(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                            twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 9:
      ibs = BjorkenMtingwa(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                           twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 10:
      ibs = BjorkenMtingwatailcut(pnumber, ex[i], ey[i], sigs[i], sige[i],
                                  twiss, twissdata, r0, aatom);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 11:
      ibs = ConteMartini(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                         twissdata, r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 12:
      ibs = ConteMartinitailcut(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss,
                                twissdata, r0, aatom);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    case 13:
      ibs = MadxIBS(pnumber, ex[i], ey[i], sigs[i], sige[i], twiss, twissdata,
                    r0);
      aes = ibs[0];
      aex = ibs[1];
      aey = ibs[2];
      break;
    }

    double ratio_x = tauradx * aex;
    double ratio_y = taurady * aey;
    double ratio_s = taurads * aes;

    printf("%-20s : %12.6e\n", "Ratio x", ratio_x);
    printf("%-20s : %12.6e\n", "Ratio y", ratio_y);
    printf("%-20s : %12.6e\n", "Ratio s", ratio_s);

    double xfactor = 1.0 / (1.0 - ratio_x);
    double yfactor = 1.0 / (1.0 - ratio_y);
    double sfactor = 1.0 / (1.0 - ratio_s);

    if (ratio_x > 1.0 || ratio_y > 1.0 || ratio_s > 1.0) {
      stepsize /= 2.0;
    }
    i++;
    t.push_back(t[i - 1] + stepsize);
    ex.push_back(ex[i - 1] + stepsize * (xfactor * equi[3] - ex[i - 1]));
    ey.push_back(ey[i - 1] +
                 stepsize * (((1.0 - coupling) * yfactor + coupling * xfactor) *
                                 equi[4] -
                             ey[i - 1]));
    sige.push_back(sige[i - 1] +
                   stepsize * (sfactor * sqrt(equi[5]) - sige[i - 1]));

    sigs.push_back(sigsfromsige(sige[i], gamma, gammatr, omegas));

  } while (i < nsteps);
  blue();
  printf("%-20s : %12.6e\n", "Final ex", ex[ex.size() - 1]);
  printf("%-20s : %12.6e\n", "Final ey", ey[ey.size() - 1]);
  printf("%-20s : %12.6e\n", "Final sigs", sigs[sigs.size() - 1]);

  printf("%-20s : %12.6e\n", "Final tau_ibs_x", 1.0 / ibs[1]);
  printf("%-20s : %12.6e\n", "Final tau_ibs_y", 1.0 / ibs[2]);
  printf("%-20s : %12.6e\n", "Final tau_ibs_s", 1.0 / ibs[0]);
  reset();
}
