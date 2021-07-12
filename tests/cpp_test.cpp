#include <ibs>
#include <map>
#include <stdio.h>
#include <string>
#include <vector>

void red() { printf("\033[1;31m"); }
void yellow() { printf("\033[1;33m"); }
void green() { printf("\033[1;32m"); }
void blue() { printf("\033[1;34m"); }
void cyan() { printf("\033[1;36m"); }
void reset() { printf("\033[0m"); }

int main() {
  string twissfilename = "b2_design_lattice_1996.twiss";
  map<string, double> twiss;
  map<string, vector<double>> twissdata;

  twiss = GetTwissHeader(twissfilename);
  twissdata = GetTwissTableAsMap(twissfilename);
  updateTwiss(twissdata);
  double aatom = emass / pmass;

  int nrf = 1;
  double harmon[1];
  double voltages[1];

  harmon[0] = 400.;
  voltages[0] = -4. * 375e3;

  double twiss_rad[6];
  double *radint;

  double gamma = twiss["GAMMA"];
  double pc = twiss["PC"];
  double gammatr = twiss["GAMMATR"];
  double mass = twiss["MASS"];
  double charge = twiss["CHARGE"];
  double q1 = twiss["Q1"];
  double len = twiss["LENGTH"];
  double epsilon = 1.0e-6;

  double betar = BetaRelativisticFromGamma(1.0);
  double r0 = ParticleRadius(1, aatom);
  double trev = len / (betar * clight);
  double frev = 1.0 / trev;
  double omega = 2.0 * pi * frev;
  double neta = eta(gamma, gammatr);

  double sigt = 0.005;

  radint = RadiationDampingLattice(twissdata);
  double U0 = RadiationLossesPerTurn(twiss, radint[1], aatom);
  double phis =
      SynchronuousPhase(0.0, 173, U0, charge, nrf, harmon, voltages, epsilon);
  double qs =
      SynchrotronTune(omega, U0, charge, nrf, harmon, voltages, phis, neta, pc);
  double omegas = qs * omega;
  double *equi =
      RadiationDampingLifeTimesAndEquilibriumEmittancesWithPartitionNumbers(
          twiss, radint, aatom, qs);

  printradint(radint);

  blue();
  printf("%-20s : %12.6f\n", "U0", U0);
  printf("%-20s : %12.6f\n", "Phis", phis);
  reset();

  return 0;
}
