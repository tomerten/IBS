#include <algorithm>
#include <map>
#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;
/**
 * Method to write simulation output to file.
 *
 * @param filename name of the file to write to
 * @param t timesteps
 * @param ex horizontal emit steps
 * @param ey vertical emit steps
 * @param sigs bunch length steps
 *
 */
void WriteToFile(string filename, vector<double> &t, vector<double> &ex,
                 vector<double> &ey, vector<double> &sigs);

/**
 * Run ODE simulation using auto time step.
 *
 * @param twiss Twiss Header Map
 * @param twissdata Twiss Table Map
 * @param nrf number of rf systems
 * @param harmon list of harmonic numbers for the rf systems
 * @param voltages list of voltages for the rf systems
 * @param[in, out] t timesteps
 * @param[in, out] ex horizontal emittance
 * @param[in, out] ey vertical emittance
 * @param[in, out] sigs bunch length
 * @param[in, out] sige energy spread
 * @param model IBS model (1-13)
 * @param pnumber number of particles per bunch
 * @param couplingpercentage hor/ver coupling in percentage
 * @param threshold evolution stop threshold
 * @param method simulation method (rlx or der)
 * @note: Relaxation method based on implementation in BMAD
 *
 */
void ODE(map<string, double> &twiss, map<string, vector<double>> &twissdata,
         int nrf, double harmon[], double voltages[], vector<double> &t,
         vector<double> &ex, vector<double> &ey, vector<double> &sigs,
         vector<double> sige, int model, double pnumber, int couplingpercentage,
         double threshold, string method);
/**
 *
 * Run ODE simulation using auto time step.
 *
 * @param twiss Twiss Header Map
 * @param twissdata Twiss Table Map
 * @param nrf number of rf systems
 * @param harmon list of harmonic numbers for the rf systems
 * @param voltages list of voltages for the rf systems
 * @param[in, out] t timesteps
 * @param[in, out] ex horizontal emittance
 * @param[in, out] ey vertical emittance
 * @param[in, out] sigs bunch length
 * @param[in, out] sige energy spread
 * @param model IBS model (1-13)
 * @param pnumber number of particles per bunch
 * @param nsteps number of simulation steps
 * @param stepsize time step size
 * @param couplingpercentage hor/ver coupling in percentage
 * @param method simulation method (rlx or der)
 *
 */
void ODE(map<string, double> &twiss, map<string, vector<double>> &twissdata,
         int nrf, double harmon[], double voltages[], vector<double> &t,
         vector<double> &ex, vector<double> &ey, vector<double> &sigs,
         vector<double> sige, int model, double pnumber, int nsteps,
         double stepsize, int couplingpercentage, string method);
