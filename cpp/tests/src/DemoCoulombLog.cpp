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
  double sigs = 0.005;
  double neta = eta(gammar, gammatr);
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
  double sige0 = sigefromsigs(omega, sigs, qs, gammar, gammatr);
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

  // va
  double pnumber = 1e10;
  double ex = 5e-9;
  double ey = 1e-10;
  double dxavg = twissheadermap["DXRMS"];
  double dyavg = twissheadermap["DYRMS"];
  double charge = twissheadermap["CHARGE"];
  double en0 = twissheadermap["ENERGY"];
  double mass = twissheadermap["MASS"];
  double sigt = sigs / clight;
  bool printout = true;
  radint = RadiationDampingLattice(twisstablemap);
  double *equi =
      RadiationDampingLifeTimesAndEquilibriumEmittancesWithPartitionNumbers(
          twissheadermap, radint, aatom, omegas);

  blue();
  printf("\nCoulombLog Functions\n");
  printf("====================\n");
  reset();

  yellow();
  printf("\nSettings used\n");
  reset();

  printf("%-30s %20.6e (%s)\n", "# particles :", pnumber, "");
  printf("%-30s %20.6e (%s)\n", "ex :", ex, "rad m pi");
  printf("%-30s %20.6e (%s)\n", "ey :", ey, "rad m pi");
  printf("%-30s %20.6e (%s)\n", "en0 :", en0, "GeV");
  printf("%-30s %20.6e (%s)\n", "mass :", mass, "GeV");

  double clog[2];
  green();
  printf("\nCoulombLog using ring averages...\n");
  cyan();
  printf("\nWithout Tailcut...\n");
  reset();

  CoulombLog(pnumber, equi[3], equi[4], twissheadermap, sige0, sigs, r0,
             printout, clog);

  cyan();
  printf("\nWith Tailcut...\n");
  reset();
  TailCutCoulombLog(pnumber, equi[3], equi[4], twissheadermap, sige0, sigs,
                    equi[0], equi[1], equi[2], r0, printout, clog);

  green();
  printf("\nAveraging element by element. \n");
  reset();

  double avgclog = 0.0;
  for (int i = 0; i < twisstablemap["L"].size(); ++i) {
    twclog(pnumber, twisstablemap["BETX"][i], twisstablemap["BETY"][i],
           twisstablemap["DX"][i], twisstablemap["DY"][i], equi[3], equi[4], r0,
           gammar, charge, en0, mass, sige0, sigs, clog);
    avgclog += twisstablemap["L"][i] * clog[0];
  }
  avgclog /= twissheadermap["LENGTH"];
  cyan();
  printf("\nWithout Tailcut...\n");
  reset();
  printf("CoulombLog avg : %12.6e\n", avgclog);

  avgclog = 0.0;
  for (int i = 0; i < twisstablemap["L"].size(); ++i) {
    twclogtail(pnumber, twisstablemap["L"][i], twisstablemap["BETX"][i],
               twisstablemap["BETY"][i], twisstablemap["DX"][i],
               twisstablemap["DPX"][i], twisstablemap["DY"][i],
               twisstablemap["DPY"][i], twisstablemap["ALFX"][i],
               twisstablemap["ALFY"][i], twisstablemap["ANGLE"][i],
               twisstablemap["K1L"][i], twisstablemap["K1SL"][i], equi[3],
               equi[4], r0, aatom, gammar, en0, len, mass, charge, sige0, sigs,
               clog);
    avgclog += twisstablemap["L"][i] * clog[0];
  }
  avgclog /= twissheadermap["LENGTH"];

  cyan();
  printf("\nWith Tailcut...\n");
  reset();
  printf("CoulombLog avg : %12.6e\n", avgclog);
  return 0;
}