#include <functional>

/**
 * Simpson integrator with decade splitting.
 *
 * @param a  lambda**2 coefficient integral denominator
 * @param b  lambda coefficient integral denominator
 * @param c  constant integral denominator
 * @param cl longitudinal growth time factor
 * @param cx horizontal growth time factor
 * @param cy vertical growht time factor
 * @param cprime scaling factor
 * @param cyy scaling factor adapted to sqrt denominator
 * @param tl1 longitudinal lambda coefficient integral numerator
 * @param tl2 longitudinal constant term integral numerator
 * @param tx1 horizontal lambda coefficient integral numerator
 * @param tx2 horizontal constant term integral numerator
 * @param ty1 vertical lambda coefficient integral numerator
 * @param ty2 vertical constant term integral numerator
 * @param[in,out] tau outputArray
 *  0 -> IBS amplitude growth rate longitudinal
 *  1 -> IBS amplitude growth rate horizontal
 *  2 -> IBS amplitude growth rate vertical
 */
void SimpsonDecade(double a, double b, double c, double cl, double cx,
                   double cy, double cprime, double cyy, double tl1, double tl2,
                   double tx1, double tx2, double ty1, double ty2, double *tau);

/**
 * The IBS integral integrand function.
 *
 * @param lambda integration variable
 * @param ax lambda coefficient integral numerator
 * @param bx constant term integral numerator
 * @param a lambda**2 coefficient integral denominator
 * @param b lambda coefficient integral denominator
 * @param c constant integral denominator
 *
 * @return integrand value
 *
 * @note CERN NOTE CERN-AB-2006-002 EQ 8
 */
double IBSIntegralIntegrand(double lambda, double ax, double bx, double a,
                            double b, double c);

/**
 * Standard Simpson integration.
 *
 * @param ibsintegrand IBS integrand as in CERN NOTE CERN-AB-2006-002 EQ 8
 * @param ax lambda coefficient integral numerator
 * @param bx constant term integral numerator
 * @param a lambda**2 coefficient integral denominator
 * @param b lambda coefficient integral denominator
 * @param c constant integral denominator
 * @param al lower integration boundary
 * @param bl higher integration boundary
 * @param n number of interval splits for the Simpson method
 *
 * @return Simpson integral value
 *
 */
double simpson(const std::function<double(double, double, double, double,
                                          double, double)> &ibsintegrand,
               double ax, double bx, double a, double b, double c, double al,
               double bl, int n);

/**
 * Standard Simpson integration.
 *
 * @param BjorkenMtingwaIntegrand IBS Bjorken-Mtingwa integrand as in CERN NOTE
 * CERN-AB-2006-002
 * @param ax horizontal lambda coefficient integral numerator
 * @param bx horizontal constant term integral numerator
 * @param ay vertical lambda coefficient integral numerator
 * @param by vertical constant term integral numerator
 * @param as longitudinal lambda coefficient integral numerator
 * @param bs longitudinal constant term integral numerator
 * @param a lambda**2 coefficient integral denominator
 * @param b lambda coefficient integral denominator
 * @param c constant integral denominator
 * @param[in, out] integral outputArray
 *  0 -> IBS longitudinal amplitude growth rate
 *  1 -> IBS horizontal amplitude growth rate
 *  2 -> IBS vertical amplitude growth rate
 *
 */
void intSimpson(
    const std::function<double(double, double, double, double, double, double)>
        &BjorkenMtingwaIntegrand,
    double ax, double bx, double ay, double by, double as, double bs, double a,
    double b, double ci, double *integral);

/**
 * Zimmerman IBS model growth rates adapted from MADX implementation.
 *
 * @param pnumber number of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param sige energy spread dE/E
 * @param gammas relativistic gamma
 * @param betax horizontal beta function
 * @param betay vertical beta function
 * @param alx horizontal alpha
 * @param aly vertical alpha
 * @param dx horizontal dispersion
 * @param dpx derivative horizontal dispersion
 * @param dy vertical dispersion
 * @param dpy derivative vertical dispersion
 * @param[in, out] tau output Array
 *  0 -> IBS longitudinal amplitude growth rate
 *  1 -> IBS horizontal amplitude growth rate
 *  2 -> IBS vertical amplitude growth rate
 *
 * @note Based on implementation of twsint in MADX (copyright CERN) CERN NOTE
 * CERN-AB-2006-002 EQ 8.
 *
 */
void twsint(double pnumber, double ex, double ey, double sigs, double sige,
            double gammas, double betax, double betay, double alx, double aly,
            double dx, double dpx, double dy, double dpy, double *tau);

/*
==============================================================================

MODEL INTEGRATORS

==============================================================================
*/
/**
 * BJorken-Mtingwa IBS model growth rates using Simpson Decade Integratino.
 *
 * @param pnumber number of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param sige energy spread dE/E
 * @param gammas relativistic gamma
 * @param betax horizontal beta function
 * @param betay vertical beta function
 * @param alx horizontal alpha
 * @param aly vertical alpha
 * @param dx horizontal dispersion
 * @param dpx derivative horizontal dispersion
 * @param dy vertical dispersion
 * @param dpy derivative vertical dispersion
 * @param[in, out] tau output Array
 *  0 -> IBS longitudinal amplitude growth rate
 *  1 -> IBS horizontal amplitude growth rate
 *  2 -> IBS vertical amplitude growth rate
 *
 * @note CERN NOTE CERN-AB-2006-002 EQ 8.
 *
 */
void BjorkenMtingwaInt(double pnumber, double ex, double ey, double sigs,
                       double sige, double gammas, double betx, double bety,
                       double alx, double aly, double dx, double dpx, double dy,
                       double dpy, double *tau);

/**
 * Conte-Martini IBS model growth rates using Simpson Decade Integratino.
 *
 * @param pnumber number of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param sige energy spread dE/E
 * @param gammas relativistic gamma
 * @param betax horizontal beta function
 * @param betay vertical beta function
 * @param alx horizontal alpha
 * @param aly vertical alpha
 * @param dx horizontal dispersion
 * @param dpx derivative horizontal dispersion
 * @param dy vertical dispersion
 * @param dpy derivative vertical dispersion
 * @param[in, out] tau output Array
 *  0 -> IBS longitudinal amplitude growth rate
 *  1 -> IBS horizontal amplitude growth rate
 *  2 -> IBS vertical amplitude growth rate
 *
 * @note CERN NOTE CERN-AB-2006-002 EQ 8.
 *
 */
void ConteMartiniInt(double pnumber, double ex, double ey, double sigs,
                     double sige, double gammas, double betx, double bety,
                     double alx, double aly, double dx, double dpx, double dy,
                     double dpy, double *tau);

/**
 * Zimmerman IBS model growth rates using Simpson Decade Integration.
 *
 * @param pnumber number of real particles in the bunch
 * @param ex horizontal emittance
 * @param ey vertical emittance
 * @param sigs bunch length
 * @param sige energy spread dE/E
 * @param gammas relativistic gamma
 * @param betax horizontal beta function
 * @param betay vertical beta function
 * @param alx horizontal alpha
 * @param aly vertical alpha
 * @param dx horizontal dispersion
 * @param dpx derivative horizontal dispersion
 * @param dy vertical dispersion
 * @param dpy derivative vertical dispersion
 * @param[in, out] tau output Array
 *  0 -> IBS longitudinal amplitude growth rate
 *  1 -> IBS horizontal amplitude growth rate
 *  2 -> IBS vertical amplitude growth rate
 *
 * @note  CERN NOTE CERN-AB-2006-002 EQ 8.
 *
 */
void MadxInt(double pnumber, double ex, double ey, double sigs, double sige,
             double gammas, double betx, double bety, double alx, double aly,
             double dx, double dpx, double dy, double dpy, double *tau);
