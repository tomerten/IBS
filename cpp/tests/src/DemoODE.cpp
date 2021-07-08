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
  /*
  ================================================================================
  INPUT
  ================================================================================
  */
  // Twiss
  string twissfilename = "../src/b2_design_lattice_1996.twiss";

  // atomic mass
  double aatom = emass / pmass;

  // RF
  int nrf = 1;
  double harmon[1];
  double voltages[1];

  harmon[0] = 400.;
  voltages[0] = -4. * 375e3;

  /*
  ================================================================================
  READ FROM INPUT FILES
  ================================================================================
  */
  // read twiss header as map and the table as map of vectors
  map<string, double> twissheadermap;
  map<string, vector<double>> twisstablemap;

  twissheadermap = GetTwissHeader(twissfilename);
  twisstablemap = GetTwissTableAsMap(twissfilename);

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

  // synchronuous phase eps
  double epsilon = 1.0e-6;

  // updateTwiss
  updateTwiss(twisstablemap);

  // printTwissMap("I2", twisstablemap);
  double betar = BetaRelativisticFromGamma(gammar);
  double r0 = ParticleRadius(1, aatom);
  double trev = len / (betar * clight);
  double frev = 1.0 / trev;
  double omega = 2.0 * pi * frev;
  double neta = eta(gammar, gammatr);

  // bunch length
  double sigt = 0.005;

  // impedance parameter
  // TODO: add explanation
  double Lpwd = 1.0e6;

  // Radiation integrals
  double *radint;
  radint = RadiationDampingLattice(twisstablemap);

  // Longitudinal parameters
  double U0 = RadiationLossesPerTurn(twissheadermap, radint[1], emass / pmass);
  double phis =
      SynchronuousPhase(0.0, 173, U0, -1, 1, harmon, voltages, epsilon);
  double qs =
      SynchrotronTune(omega, U0, -1, 1, harmon, voltages, phis, neta, pc);
  double omegas = qs * omega;
  double sige0 = sigefromsigs(omega, sigt, qs, gammar, gammatr);

  // Radiation - Excitation equilbria
  double *equi =
      RadiationDampingLifeTimesAndEquilibriumEmittancesWithPartitionNumbers(
          twissheadermap, radint, aatom, omegas);

  double pnumber = 1e10;
  double ex = 5e-9;
  double ey = 1e-10;
  double bxavg = len / (2.0 * pi * twissheadermap["Q1"]);
  double byavg = len / (2.0 * pi * twissheadermap["Q2"]);
  double dxavg = twissheadermap["DXRMS"];
  double dyavg = twissheadermap["DYRMS"];
  double charge = twissheadermap["CHARGE"];
  double sige = 1e-4;
  double en0 = twissheadermap["ENERGY"];
  double mass = twissheadermap["MASS"];
  bool printout = true;

  /*
  ================================================================================
  ODE
  ================================================================================
  */
  printf("\n\n");
  blue();
  printf("ODE Model\n");
  printf("==========\n");
  reset();

  vector<double> t;
  vector<double> exa;
  vector<double> eya;
  vector<double> sigsa;
  vector<double> sigea;

  int maxsteps = 2;
  double timestep = 1.0e-3;

  t.push_back(0.0);
  exa.push_back(7.5e-9);
  eya.push_back(1e-9);
  sigsa.push_back(5e-3);

  ODE(twissheadermap, twisstablemap, 1, harmon, voltages, t, exa, eya, sigsa,
      sigea, 1, pnumber, 0, 1e-4);
  WriteToFile("ODE_test_output_piwinski_smooth.csv", t, exa, eya, sigsa);

  t.clear();
  exa.clear();
  eya.clear();
  sigsa.clear();
  t.push_back(0.0);
  exa.push_back(7.5e-9);
  eya.push_back(1e-9);
  sigsa.push_back(5e-3);

  ODE(twissheadermap, twisstablemap, 1, harmon, voltages, t, exa, eya, sigsa,
      sigea, 2, pnumber, 0, 1e-4);
  WriteToFile("ODE_test_output_piwinski_lattice.csv", t, exa, eya, sigsa);

  t.clear();
  exa.clear();
  eya.clear();
  sigsa.clear();
  t.push_back(0.0);
  exa.push_back(7.5e-9);
  eya.push_back(1e-9);
  sigsa.push_back(5e-3);

  ODE(twissheadermap, twisstablemap, 1, harmon, voltages, t, exa, eya, sigsa,
      sigea, 3, pnumber, 0, 1e-4);
  WriteToFile("ODE_test_output_piwinski_latticemodified.csv", t, exa, eya,
              sigsa);

  t.clear();
  exa.clear();
  eya.clear();
  sigsa.clear();
  t.push_back(0.0);
  exa.push_back(7.5e-9);
  eya.push_back(1e-9);
  sigsa.push_back(5e-3);

  ODE(twissheadermap, twisstablemap, 1, harmon, voltages, t, exa, eya, sigsa,
      sigea, 4, pnumber, 0, 1e-4);
  WriteToFile("ODE_test_output_nagaitsev.csv", t, exa, eya, sigsa);

  t.clear();
  exa.clear();
  eya.clear();
  sigsa.clear();
  t.push_back(0.0);
  exa.push_back(7.5e-9);
  eya.push_back(1e-9);
  sigsa.push_back(5e-3);

  ODE(twissheadermap, twisstablemap, 1, harmon, voltages, t, exa, eya, sigsa,
      sigea, 5, pnumber, 0, 1e-4);
  WriteToFile("ODE_test_output_nagaitsevtailcut.csv", t, exa, eya, sigsa);

  t.clear();
  exa.clear();
  eya.clear();
  sigsa.clear();
  t.push_back(0.0);
  exa.push_back(7.5e-9);
  eya.push_back(1e-9);
  sigsa.push_back(5e-3);

  ODE(twissheadermap, twisstablemap, 1, harmon, voltages, t, exa, eya, sigsa,
      sigea, 6, pnumber, 0, 1e-4);
  WriteToFile("ODE_test_output_madx.csv", t, exa, eya, sigsa);

  t.clear();
  exa.clear();
  eya.clear();
  sigsa.clear();
  t.push_back(0.0);
  exa.push_back(7.5e-9);
  eya.push_back(1e-9);
  sigsa.push_back(5e-3);

  ODE(twissheadermap, twisstablemap, 1, harmon, voltages, t, exa, eya, sigsa,
      sigea, 7, pnumber, 0, 1e-4);
  WriteToFile("ODE_test_output_madxtailcut.csv", t, exa, eya, sigsa);

  t.clear();
  exa.clear();
  eya.clear();
  sigsa.clear();
  t.push_back(0.0);
  exa.push_back(7.5e-9);
  eya.push_back(1e-9);
  sigsa.push_back(5e-3);

  ODE(twissheadermap, twisstablemap, 1, harmon, voltages, t, exa, eya, sigsa,
      sigea, 8, pnumber, 0, 1e-4);
  WriteToFile("ODE_test_output_bjorken_mtingwa2.csv", t, exa, eya, sigsa);

  t.clear();
  exa.clear();
  eya.clear();
  sigsa.clear();
  t.push_back(0.0);
  exa.push_back(7.5e-9);
  eya.push_back(1e-9);
  sigsa.push_back(5e-3);

  ODE(twissheadermap, twisstablemap, 1, harmon, voltages, t, exa, eya, sigsa,
      sigea, 9, pnumber, 0, 1e-4);
  WriteToFile("ODE_test_output_bjorken_mtingwa.csv", t, exa, eya, sigsa);

  t.clear();
  exa.clear();
  eya.clear();
  sigsa.clear();
  t.push_back(0.0);
  exa.push_back(7.5e-9);
  eya.push_back(1e-9);
  sigsa.push_back(5e-3);

  ODE(twissheadermap, twisstablemap, 1, harmon, voltages, t, exa, eya, sigsa,
      sigea, 10, pnumber, 0, 1e-4);
  WriteToFile("ODE_test_output_bjorken_mtingwatailcut.csv", t, exa, eya, sigsa);

  t.clear();
  exa.clear();
  eya.clear();
  sigsa.clear();
  t.push_back(0.0);
  exa.push_back(7.5e-9);
  eya.push_back(1e-9);
  sigsa.push_back(5e-3);

  ODE(twissheadermap, twisstablemap, 1, harmon, voltages, t, exa, eya, sigsa,
      sigea, 11, pnumber, 0, 1e-4);
  WriteToFile("ODE_test_output_conte_martini.csv", t, exa, eya, sigsa);

  t.clear();
  exa.clear();
  eya.clear();
  sigsa.clear();
  t.push_back(0.0);
  exa.push_back(7.5e-9);
  eya.push_back(1e-9);
  sigsa.push_back(5e-3);

  ODE(twissheadermap, twisstablemap, 1, harmon, voltages, t, exa, eya, sigsa,
      sigea, 12, pnumber, 0, 1e-4);
  WriteToFile("ODE_test_output_conte_martini_tailcut.csv", t, exa, eya, sigsa);

  t.clear();
  exa.clear();
  eya.clear();
  sigsa.clear();
  t.push_back(0.0);
  exa.push_back(7.5e-9);
  eya.push_back(1e-9);
  sigsa.push_back(5e-3);

  ODE(twissheadermap, twisstablemap, 1, harmon, voltages, t, exa, eya, sigsa,
      sigea, 13, pnumber, 0, 1e-4);
  WriteToFile("ODE_test_output_madxibs.csv", t, exa, eya, sigsa);

  return 0;
}