#include "NumericFunctions.hpp"
#include "RadiationDamping.hpp"
#include <map>
#include <math.h>
#include <stdio.h>
#include <string>

/**
 * Method to calculate Coulomb Log and constant for single lattice element.
 *
 * @note Adapted from original code in MADX (copyright CERN). Based on the
 * formulae in AIP physics vade mecum p.264 (1981).
 *
 * @param pnumber Number of particles in the bunch
 * @param bx Horizontal beta function
 * @param by Vertical beta function
 * @param dx Horizontal dispersion
 * @param dy Vertical dispersion
 * @param ex Horizontal emittance
 * @param ey Vertical emittance
 * @param r0 Classical radius of the particles in the bunch
 * @param gamma Relativistic Gamma
 * @param charge Particle electric charge
 * @param en0 Energy of the particle
 * @param amass Rest Energy of the partice [GeV]
 * @param sige Energy spread
 * @param sigt Bunch length
 * @param[out] clog Ouput array
 *  0 -> Coulomb Log
 *  1 -> IBS Coulomb pre-factor
 *
 */
void twclog(double pnumber, double bx, double by, double dx, double dy,
            double ex, double ey, double r0, double gamma, double charge,
            double en0, double amass, double sige, double sigt, double *clog);

/**
 * Method to calculate Coulomb Log and constant for single lattice element with
 * Tailcut (taking radiation damping into account).
 *
 * @note Based on the
 * formulae in AIP physics vade mecum p.264 (1981).
 *
 * @param pnumber Number of particles in the bunch
 * @param l Element length
 * @param bx Horizontal beta function
 * @param by Vertical beta function
 * @param dx Horizontal dispersion
 * @param dpx Horizontal dispersion derivative
 * @param dy Vertical dispersion
 * @param dpy Vertical dispersion derivative
 * @param ax Horizontal alpha
 * @param ay Vertical alpha
 * @param angle Bending angle
 * @param k1l Quadrupole component
 * @param k1sl Skew Quadrupole component
 * @param ex Horizontal emittance
 * @param ey Vertical emittance
 * @param r0 Classical radius of the particles in the bunch
 * @param aatom Atomic Mass Number
 * @param gamma Relativistic Gamma
 * @param en0 Energy of the particle
 * @param len Accelerator length
 * @param amass Rest Energy of the partice [GeV]
 * @param charge Particle electric charge
 * @param sige Energy spread
 * @param sigt Bunch length
 * @param[out] clog Ouput array
 *  0 -> Coulomb Log
 *  1 -> IBS Coulomb pre-factor
 *
 */
void twclogtail(double pnumber, double l, double bx, double by, double dx,
                double dpx, double dy, double dpy, double ax, double ay,
                double angle, double k1l, double k1sl, double ex, double ey,
                double r0, double aatom, double gamma, double en0, double len,
                double amass, double charge, double sige, double sigt,
                double *clog);

/**
 * Method to calculate Coulomb Log and constant using accelerator ring averages.
 *
 * @note Based on the
 * formulae in AIP physics vade mecum p.264 (1981).
 *
 * @param pnumber Number of particles in the bunch
 * @param ex Horizontal emittance
 * @param ey Vertical emittance
 * @param twissheader Twiss Header Map
 * @param sige Energy spread
 * @param sigt Bunch length
 * @param r0 Classical radius of the particles in the bunch
 * @param printout Boolean flag to switch verbose mode on or off
 * @param[out] clog Ouput array
 *  0 -> Coulomb Log
 *  1 -> IBS Coulomb pre-factor
 *
 */
void CoulombLog(double pnumber, double ex, double ey,
                map<string, double> &twissheader, double sige, double sigt,
                double r0, bool printout, double *clog);

/**
 * Method to calculate Coulomb Log and constant using accelerator ring averages
 * with tailcut.
 *
 * @note Based on the
 * formulae in AIP physics vade mecum p.264 (1981).
 *
 * @param pnumber Number of particles in the bunch
 * @param ex Horizontal emittance
 * @param ey Vertical emittance
 * @param twissheader Twiss Header Map
 * @param sige Energy spread
 * @param sigt Bunch length
 * @param tauradx Horizontal Radiation Damping Time
 * @param taurady Vertical Radiation Damping Time
 * @param taurads Longitudinal Radiation Damping Time
 * @param r0 Classical radius of the particles in the bunch
 * @param printout Boolean flag to switch verbose mode on or off
 * @param[out] clog Ouput array
 *  0 -> Coulomb Log
 *  1 -> IBS Coulomb pre-factor
 *
 */

void TailCutCoulombLog(double pnumber, double ex, double ey,
                       map<string, double> &twissheader, double sige,
                       double sigt, double tauradx, double taurady,
                       double taurads, double r0, bool printout, double *clog);