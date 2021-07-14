#ifndef NUMERIC_FUNCTIONS_HPP
#define NUMERIC_FUNCTIONS_HPP
#include <map>
#include <stdio.h>
#include <string>
#include <vector>
using namespace std;
/*
================================================================================
CONSTANTS
================================================================================
*/
/// Speed of light
const double clight = 299792458.0;
/// Reduced Planck constant  [GeV]
const double hbar = 6.582119569e-25;
/// Electron mass [GeV]
const double emass = 0.51099895000e-3;
/// Proton mass [GeV]
const double pmass = 0.93827208816;
/// Neutron mass [Gev]
const double nmass = 0.93956542052; // GeV CODATA 2018
/// Muon mass [GeV]
const double mumass = 0.1056583755; // GeV CODATA 2018
/// Atomic mass unit  [GeV]
const double atomicmassunit = 0.93149410242; // GeV scipy constants
/// Pi
const double pi = 3.141592653589793;
/// Electric Charge unit [Coulomb]
const double ec = 1.602176634e-19;
/// Euler constant
const double euler = 0.577215664901533;
/// Classical radius [m]
const double erad = 2.8179403262e-15;
/// Classical proton radius [m]
const double prad = erad * emass / pmass;

/*
================================================================================
BASIC FUNCTIONS
================================================================================
*/

/**
 * Add ratiation integrals (I1-I5) per element to the Twiss Table Map. It also
 adds the relevant quantities to calculate them (e.g. curly H, rho, ...).
 *
 * @param[in, out] table Twiss Table Map
 *
 */
void updateTwiss(map<string, vector<double>> &table);

/**
 * Prints the values of the given Twiss column name. Mostly for debugging.
 *
 * @param key Name of the column to print out.
 * @param table Twiss Table Map to print from.
 *
 */
void printTwissMap(string key, map<string, vector<double>> &table);

/**
 * Calculate relativistic beta given relativistic gamma.
 *
 * @param gamma Relativistic gamma factor.
 *
 * @return Relativistic Beta
 */
double BetaRelativisticFromGamma(double gamma);

/**
 * Calculate slip factor eta, from gamma and gamma transition.
 *
 * @param gamma Relativistic Gamma of the particles.
 * @param gammatr Transition Relativistic Gamma of the accelerator lattice.
 *
 * @return Slip factor Eta.
 */
double eta(double gamma, double gammatr);

/**
 * Calculate the classical particle radius from its charge and atomic mass
 * number A.
 *
 * @param charge Charge of the particle.
 * @param aatom Atomic mass number of the particle (A)
 *
 * @return Classical radius of the particle.
 */
double ParticleRadius(double charge, double aatom);

/**
 * Convert dE/E to dp/p.
 *
 * @param dee dE/E
 * @param beta0 Relativistic Beta
 *
 * @return dp/p.
 */
double dee_to_dpp(double dee, double beta0);

/**
 * Convert dp/p to dE/E.
 *
 * @param dpp dp/p
 * @param beta0 Relativistic Beta
 *
 * @return dE/E.
 */
double dpp_to_dee(double dpp, double beta0);

/*
================================================================================
LONGITUDINAL / RF FUNCTIONS
================================================================================
*/
/**
 * Convert @f$\sigma_s @f$ to @f$\sigma_e @f$
 *
 * @param omega0 Particle revolution angular frequency (@f$ 2 \pi f  @f$ )
 * @param sigs @f$\sigma_s @f$
 * @param qs Synchrotron Tune
 * @param gamma Relativistic Gamma
 * @param gammatr Transition Relativistic Gamma
 *
 * @return @f$\sigma_e @f$
 */
double sigefromsigs(double omega0, double sigs, double qs, double gamma,
                    double gammatr);

/**
 * Convert @f$\sigma_e @f$ to @f$\sigma_s @f$
 *
 * @param omega0 Particle revolution angular frequency (@f$ 2 \pi f  @f$ )
 * @param sige @f$\sigma_e @f$
 * @param gamma Relativistic Gamma
 * @param gammatr Transition Relativistic Gamma
 * @param omegas Synchrotron Frequency (@f$ \omega q_s @f$)
 *
 * @return @f$\sigma_s @f$
 */
double sigsfromsige(double sige, double gamma, double gammatr, double omegas);

/**
 * Get the energy gain in eV produced by the RF systems for given harmonic
 * numbers and voltages.
 *
 * @param phi RF phase for the particle
 * @param charge Charge of the particle
 * @param nrf Number of RF systems
 * @param harmon List of harmonic numbers
 * @param voltages List of voltages
 *
 * @return Energy gain [eV]
 */
double EffectiveRFVoltageInElectronVolt(double phi, double charge, int nrf,
                                        double harmon[], double voltages[]);

/**
 * Get the derivative of the energy gain in eV produced by the RF systems for
 * given harmonic numbers and voltages.
 *
 * @param phi RF phase for the particle
 * @param charge Charge of the particle
 * @param nrf Number of RF systems
 * @param harmon List of harmonic numbers
 * @param voltages List of voltages
 *
 * @return Derivative of the Energy gain [eV]
 */
double EffectiveRFVoltageInElectronVoltPrime(double phi, double charge, int nrf,
                                             double harmon[],
                                             double voltages[]);

/**
 * Get the energy gain in eV produced by the RF systems for
 * given harmonic numbers and voltages taking radiation losses into account.
 *
 * @param phi RF phase for the particle
 * @param U0 Radiation losses per turn [eV]
 * @param charge Charge of the particle
 * @param nrf Number of RF systems
 * @param harmon List of harmonic numbers
 * @param voltages List of voltages
 *
 * @return Energy gain [eV] taking radiation losses into account.
 */
double VeffRFeVRadlosses(double phi, double U0, double charge, int nrf,
                         double harmon[], double voltages[]);

/**
 * Get the synchronuous phase for a particle for given energy gain, RF settings
 * and radiation losses.
 *
 * @param target Target energy gain per turn [eV]
 * @param init_phi RF phase to start search from.
 * @param U0 Radiation losses per turn [eV]
 * @param charge Charge of the particle
 * @param nrf Number of RF systems
 * @param harmon List of harmonic numbers
 * @param voltages List of voltages
 * @param epsilon Accuracy threshold for the target.
 *
 * @return Synchronuous phase.
 */
double SynchronuousPhase(double target, double init_phi, double U0,
                         double charge, int nrf, double harmon[],
                         double voltages[], double epsilon);

/**
 * Get the energy gain in eV produced by the RF systems for given harmonic
 * numbers and voltages taking potetial well distortion and radiation losses
 * into account.
 *
 * @param phi RF phase for the particle
 * @param U0 Radiation losses per turn [eV]
 * @param charge Charge of the particle
 * @param nrf Number of RF systems
 * @param harmon List of harmonic numbers
 * @param voltages List of voltages
 * @param L Impedance L
 * @param N Impedance N
 * @param sigs @f$\sigma_s @f$ bunch length
 * @param pc Particle Momentum in GeV
 *
 * @return Energy gain [eV] with Potential Well Distortion and Radiation Losses.
 */
double VeffRFeVPotentialWellDistortion(double phi, double U0, double charge,
                                       int nrf, double harmon[],
                                       double voltages[], double L, double N,
                                       double sigs, double pc);

/**
 * Get the derivative of the energy gain in eV produced by the RF systems for
 * given harmonic numbers and voltages taking potetial well distortion and
 * radiation losses into account.
 *
 * @param phi RF phase for the particle
 * @param U0 Radiation losses per turn [eV]
 * @param charge Charge of the particle
 * @param nrf Number of RF systems
 * @param harmon List of harmonic numbers
 * @param voltages List of voltages
 * @param L Impedance L
 * @param N Impedance N
 * @param sigs @f$\sigma_s @f$ bunch length
 * @param pc Particle Momentum in GeV
 *
 * @return Derivative of the Energy gain [eV] with Potential Well Distortion and
 * Radiation Losses.
 */
double VeffRFeVPotentialWellDistortionPrime(double phi, double U0,
                                            double charge, int nrf,
                                            double harmon[], double voltages[],
                                            double L, double N, double sigs,
                                            double pc);

/**
 * Get the synchronuous phase with Potential Well Distortion for a particle for
 * given energy gain, RF settings and radiation losses.
 *
 * @param target Target energy gain per turn [eV]
 * @param init_phi RF phase to start search from
 * @param U0 Radiation losses per turn [eV]
 * @param charge Charge of the particle
 * @param nrf Number of RF systems
 * @param harmon List of harmonic numbers
 * @param voltages List of voltages
 * @param L Impedance L
 * @param N Impedance N
 * @param sigs @f$\sigma_s @f$ bunch length
 * @param pc Particle Momentum in GeV
 * @param epsilon Accuracy threshold for the target
 *
 * @return Synchronuous phase with  Potential Well Distortion.
 */
double SynchronuousPhaseWithPWD(double target, double init_phi, double U0,
                                double charge, int nrf, double harmon[],
                                double voltages[], double L, double N,
                                double sigs, double pc, double epsilon);

/**
 * Synchrotron Tune
 *
 * @param omega0 Angular revolution frequency of the particle.
 * @param U0 Radiation losses per turn [eV]
 * @param charge Charge of the particle
 * @param nrf Number of RF systems
 * @param harmon List of harmonic numbers
 * @param voltages List of voltages
 * @param phis Synchronuous phase of the particle
 * @param eta Phase slip factor
 * @param pc Particle Momentum [GeV]
 *
 * @return Synchrotron Tune
 */
double SynchrotronTune(double omega0, double U0, double charge, int nrf,
                       double harmon[], double voltages[], double phis,
                       double eta, double pc);

/**
 * Synchrotron Tune with Potential Well Distortion
 *
 * @param omega0 Angular revolution frequency of the particle.
 * @param U0 Radiation losses per turn [eV]
 * @param charge Charge of the particle
 * @param nrf Number of RF systems
 * @param harmon List of harmonic numbers
 * @param voltages List of voltages
 * @param L Impedance L
 * @param N Impedance N
 * @param phis Synchronuous phase of the particle
 * @param eta Phase slip factor
 * @param pc Particle Momentum [GeV]
 *
 * @return Synchrotron Tune with Potential Well Distortion
 */
double SynchrotronTunePWD(double omega0, double U0, double charge, int nrf,
                          double harmon[], double voltages[], double L,
                          double N, double sigs, double phis, double eta,
                          double pc);

/**
 * Convert @f$\sigma_s @f$ to @f$\sigma_e @f$ using RF parameters
 *
 * @param sigs @f$\sigma_s @f$
 * @param U0 Radiation losses per turn [eV]
 * @param charge Charge of the particle
 * @param nrf Number of RF systems
 * @param harmon List of harmonic numbers
 * @param voltages List of voltages
 * @param gamma Relativistic Gamma
 * @param gammatr Transition Relativistic Gamma
 * @param pc Particle Momentum [GeV]
 * @param circ Accelerator length
 * @param phis Synchronuous phase of the particle
 * @param printout Flag to print out the summary.
 *
 * @return @f$\sigma_e @f$
 */
double SigeFromRFAndSigs(double sigs, double U0, double charge, int nrf,
                         double harmon[], double voltages[], double gamma,
                         double gammatr, double pc, double circ, double phis,
                         bool printout);

/**
 * Piwinski IBS approximation fmohl function. Reference in IBS models docstring.
 *
 * @param a
 * @param b
 * @param q
 * @param n Number of iterations for accuracy (standard 1000 is ok)
 *
 * @return Fmohl
 */
double fmohl(double a, double b, double q, int n);

/**
 * Nagaitsev's approximation rds function. Reference in IBS models docstring.
 *
 * @param x
 * @param y
 * @param z
 *
 * @return rds
 */
double rds(double x, double y, double z);

#endif