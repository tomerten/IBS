#include <functional>
void SimpsonDecade(double a, double b, double c, double cl, double cx,
                   double cy, double cprime, double cyy, double tl1, double tl2,
                   double tx1, double tx2, double ty1, double ty2, double *tau);

void twsint(double pnumber, double ex, double ey, double sigs, double sige,
            double gammas, double betax, double betay, double alx, double aly,
            double dx, double dpx, double dy, double dpy, double *tau);

double IBSIntegralIntegrand(double lambda, double ax, double bx, double a,
                            double b, double c);

double simpson(const std::function<double(double, double, double, double,
                                          double, double)> &ibsintegrand,
               double ax, double bx, double a, double b, double c, double al,
               double bl, int n);

void intSimpson(const std::function<double(double, double, double, double,
                                           double, double)> &ibsintegrand,
                double ax, double bx, double ay, double by, double as,
                double bs, double a, double b, double ci, double *integral);

/*
==============================================================================

MODEL INTEGRATORS

==============================================================================
*/
void BjorkenMtingwaInt(double pnumber, double ex, double ey, double sigs,
                       double sige, double gammas, double betx, double bety,
                       double alx, double aly, double dx, double dpx, double dy,
                       double dpy, double *tau);

void ConteMartiniInt(double pnumber, double ex, double ey, double sigs,
                     double sige, double gammas, double betx, double bety,
                     double alx, double aly, double dx, double dpx, double dy,
                     double dpy, double *tau);

void MadxInt(double pnumber, double ex, double ey, double sigs, double sige,
             double gammas, double betx, double bety, double alx, double aly,
             double dx, double dpx, double dy, double dpy, double *tau);

void twsint(double pnumber, double ex, double ey, double sigs, double sige,
            double gammas, double betax, double betay, double alx, double aly,
            double dx, double dpx, double dy, double dpy, double *tau);
