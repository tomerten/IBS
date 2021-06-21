#include <ibs>
#include <map>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

void red() { printf("\033[1;31m"); }
void yellow() { printf("\033[1;33m"); }
void green() { printf("\033[1;32m"); }
void blue() { printf("\033[1;34m"); }
void cyan() { printf("\033[1;36m"); }
void reset() { printf("\033[0m"); }

int main() {
  string twissfilename = "../src/b2_design_lattice_1996.twiss";
  map<string, double> twissheadermap;
  map<string, vector<double>> twisstablemap;

  twissheadermap = GetTwissHeader(twissfilename);
  twisstablemap = GetTwissTableAsMap(twissfilename);

  double aatom = emass / pmass;
  // VeffRFeV
  int nrf = 1;
  double harmon[1];
  double voltages[1];
  harmon[0] = 400.;
  voltages[0] = -4. * 375e3;

  /*
  ================================================================================
  BASIC NUMERIC FUNCTIONS
  ================================================================================
  */
  // parameters needed in this section
  double gammar = twissheadermap["GAMMA"];
  double gammatr = twissheadermap["GAMMATR"];
  double pc = twissheadermap["PC"];
  double len = twissheadermap["LENGTH"];
  // fmohl
  double a = 5.709563671168914e-04;
  double b = 2.329156389696222e-01;
  double q = 2.272866910079534e00;
  int npp = 1000;
  // synchronuous phase eps
  double epsilon = 1.0e-6;
  double rho = 4.35;

  // updateTwiss
  updateTwiss(twisstablemap);
  // printTwissMap("I2", twisstablemap);
  double betar = BetaRelativisticFromGamma(gammar);
  double r0 = ParticleRadius(1, aatom);
  double trev = len / (betar * clight);
  double frev = 1.0 / trev;
  double omega = 2.0 * pi * frev;
  double Lpwd = 1.0e6;
  double neta = eta(gammar, gammatr);
  double sige0 = sigefromsigs(2.0 * pi * 1.2e6, 0.001, 5e-3, gammar, gammatr);
  double VrfEffeV =
      EffectiveRFVoltageInElectronVolt(173, -1, 1, harmon, voltages);
  double VrfEffeVp =
      EffectiveRFVoltageInElectronVoltPrime(173, -1, 1, harmon, voltages);
  double *radint;
  radint = RadiationDampingLattice(twisstablemap);
  double U0 = RadiationLossesPerTurn(twissheadermap, radint[1], emass / pmass);
  double VrfEffeVU0 = VeffRFeVRadlosses(173, U0, -1, 1, harmon, voltages);
  double phis =
      SynchronuousPhase(0.0, 173, U0, -1, 1, harmon, voltages, epsilon);
  double qs =
      SynchrotronTune(omega, U0, -1, 1, harmon, voltages, phis, neta, pc);
  double omegas = qs * omega;
  double VrfEffeVPWD = VeffRFeVPotentialWellDistortion(
      173, U0, -1, 1, harmon, voltages, Lpwd, 1, 0.005, pc);
  double VrfEffeVPWDp = VeffRFeVPotentialWellDistortionPrime(
      173, U0, -1, 1, harmon, voltages, Lpwd, 1, 0.005, pc);
  double phisPWD = SynchronuousPhaseWithPWD(
      0.0, 173, U0, -1, 1, harmon, voltages, Lpwd, 1, 0.005, pc, epsilon);
  double qsPWD = SynchrotronTunePWD(omega, U0, -1, 1, harmon, voltages, Lpwd, 1,
                                    0.005, phisPWD, neta, pc);
  double bxavg = len / (2.0 * pi * twissheadermap["Q1"]);
  double byavg = len / (2.0 * pi * twissheadermap["Q2"]);

  blue();
  printf("\nRadiation Damping\n");
  printf("=================\n");

  green();
  printf("\nRadiation Smooth Ring Approximation\n");
  reset();

  radint = RadiationDampingApprox(len, gammar, gammatr, 4.35, bxavg, byavg);
  printradint(radint);

  green();
  printf("\nRadiation Damping element by element\n");
  reset();

  radint = RadiationDampingLattice(twisstablemap);
  printradint(radint);

  phis = 173.49249640825525;
  qs = SynchrotronTune(omega, U0, -1, 1, harmon, voltages, phis, neta, pc);
  omegas = qs * omega;
  double *equi =
      RadiationDampingLifeTimesAndEquilibriumEmittancesWithPartitionNumbers(
          twissheadermap, radint, aatom, qs);
  green();
  printf("\nRadiation Damping Equib \n");
  reset();

  printf("%-30s %10.6e (%s)\n", "Taux :", equi[0], "s");
  printf("%-30s %10.6e (%s)\n", "Tauy :", equi[1], "s");
  printf("%-30s %10.6e (%s)\n", "Taus :", equi[2], "s");
  printf("%-30s %10.6e (%s)\n", "exinf :", equi[3], "");
  printf("%-30s %10.6e (%s)\n", "eyinf :", equi[4], "");
  printf("%-30s %10.6e (%s)\n", "sigeoe2 :", equi[5], "");
  printf("%-30s %10.6e (%s)\n", "sigt :", equi[6], "");
  printf("%-30s %10.6e (%s)\n", "jx :", equi[7], "");
  printf("%-30s %10.6e (%s)\n", "jy :", equi[8], "");

  green();
  printf("\nCritical Energy Calculations. \n");
  reset();
  double *critical;
  critical = RadiationCriticalEnergy(rho, gammar, omega);

  printf("%-30s %10.6e (%s)\n", "omega_crit :", critical[0], "");
  printf("%-30s %10.6e (%s)\n", "Theta_crit :", critical[1], "");
  printf("%-30s %10.6e (%s)\n", "E_crit :", critical[2], "eV");
  printf("%-30s %10.6e (%s)\n", "E_per_photon_avg :", critical[3], "eV/photon");
  printf("%-30s %10.6e (%s)\n", "N_photons_avg_per_turn :", critical[4],
         "1/turn");

  return 0;
}