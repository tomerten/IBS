#include "../include/ibs_bits/CoulombLogFunctions.hpp"
#include "../include/ibs_bits/NumericFunctions.hpp"
#include <functional>
#include <math.h>
#include <stdio.h>

/*
================================================================================
================================================================================
SIMPSON INTEGRATOR WITH DECADE SPLITTING TO BE ABLE TO ADD
CONTE-MARTINI AND BJORKEN-MTINGWA.
================================================================================
  AUTHORS:
    - MADX AUTHORS COPYRIGHT CERN
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF:
    - BASED ON MADX ORIGINAL SOURCE CODE IN TWSINT FUNCTION
    - CERN NOTE CERN-AB-2006-002 EQ 8

================================================================================
  Arguments:
  ----------
    - double a
        lambda**2 coefficient integral denominator
    - double b
        lambda coefficient integral denominator
    - double c
        constant integral denominator
    - double cl
        longitudinal growth time factor
    - double cx
        horizontal growth time factor
    - double cy
        vertical growht time factor
    - double cprime (double)
        scaling factor
    - double cyy
        scaling factor adapted to sqrt denominator
    - double tl1
        lambda coefficient integral numerator
    - double tl2
        constant term integral numerator
    - double tx1
        lambda coefficient integral numerator
    - double tx2
        constant term integral numerator
    - double ty1
        lambda coefficient integral numerator
    - double ty2
        constant term integral numerator

  Returns:
  --------
    double [3] tau
      IBS growth rates for AMPLITUDES (SIGMA NOT EMIT -> MULTIPLY WITH TWO FOR
      EMIT)
      0 -> al
      1 -> ax
      2 -> ay

================================================================================
================================================================================
*/
void SimpsonDecade(double a, double b, double c, double cl, double cx,
                   double cy, double cprime, double cyy, double tl1, double tl2,
                   double tx1, double tx2, double ty1, double ty2,
                   double *tau) {
  const int maxdec = 30, ns = 50;

  const double ten = 10.0;
  const double three = 3.0;
  const double test = 1.0e-7;
  const double coeff[2] = {2.0, 4.0};

  double ar[31], br[30];

  double h, aloop, alam, term, cof, f;
  double suml, sumx, sumy, tmpl, tmpx, tmpy;
  double polyl, polyx, polyy, func;
  bool flag = 0;

  double zintl = 0.0;
  double zintx = 0.0;
  double zinty = 0.0;

  ar[0] = 0.0;

  for (int iloop = 0; iloop < maxdec; iloop++) {
    br[iloop] = pow(ten, iloop);
    ar[iloop + 1] = br[iloop];
    h = (br[iloop] - ar[iloop]) / ns;
    aloop = ar[iloop];

    term = sqrt(cyy * cyy *
                (aloop * aloop * aloop + a * aloop * aloop + b * aloop + c));

    func = sqrt(aloop) / (term * term * term);
    polyl = tl1 * aloop + tl2;
    polyx = tx1 * aloop + tx2;
    polyy = ty1 * aloop + ty2;
    suml = func * polyl;
    sumx = func * polyx;
    sumy = func * polyy;

    for (int iiz = 0; iiz < ns; iiz++) {
      alam = aloop + iiz * h;
      cof = coeff[iiz % 2];
      term = sqrt(cyy * cyy *
                  (alam * alam * alam + a * alam * alam + b * alam + c));

      f = sqrt(alam) / (term * term * term);
      polyl = tl1 * alam + tl2;
      polyx = tx1 * alam + tx2;
      polyy = ty1 * alam + ty2;

      suml = suml + cof * f * polyl;
      sumx = sumx + cof * f * polyx;
      sumy = sumy + cof * f * polyy;
    }

    suml = suml - f * polyl;
    sumx = sumx - f * polyx;
    sumy = sumy - f * polyy;
    tmpl = (suml / three) * h;
    tmpx = (sumx / three) * h;
    tmpy = (sumy / three) * h;
    zintl = zintl + tmpl;
    zintx = zintx + tmpx;
    zinty = zinty + tmpy;

    if ((fabs(tmpl / zintl) < test) & (fabs(tmpx / zintx) < test) &
        (fabs(tmpy / zinty) < test)) {
      // ---- Divide answers by cprime to account for scaling.
      tau[1] = cx * (zintx / cprime);
      tau[0] = cl * (zintl / cprime);
      tau[2] = cy * (zinty / cprime);

      flag = 1;
      break;
    }
  }

  if (flag == 0) {
    tau[0] = 0.0;
    tau[1] = 0.0;
    tau[2] = 0.0;
  }
}
/*
================================================================================
================================================================================
IBS INTEGRALS INTEGRAND AS FUNCTION TO BE ABLE TO ADD
CONTE-MARTINI AND BJORKEN-MTINGWA.
================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

================================================================================
  Arguments:
  ----------
    - double a
        lambda**2 coefficient integral denominator
    - double b
        lambda coefficient integral denominator
    - double c
        constant integral denominator
    - double ax
        lambda coefficient integral numerator
    - double bx
        constant term integral numerator
    - double lambda
        integration variable

  Returns:
  --------
    double integral
        integrand for IBS integral
================================================================================
================================================================================
*/
double IBSIntegralIntegrand(double lambda, double ax, double bx, double a,
                            double b, double c) {
  double num = sqrt(lambda) * (ax * lambda + bx);
  double term = lambda * lambda * lambda + a * lambda * lambda + b * lambda + c;
  double denom = sqrt(term * term * term);
  return num / denom;
}

/*
================================================================================
================================================================================
STANDARD SIMPSON INTEGRATION METHOD

NOTE:
-----
NORMAL SIMPSON GIVES TOO
LARGE DEVIATIONS COMPARED TO SIMPSON DECADE.
================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

================================================================================
  Arguments:
  ----------
    - double ibsintegrand(double, double, double, double, double,
                                   double)
    - double a
        lambda**2 coefficient integral denominator
    - double b
        lambda coefficient integral denominator
    - double c
        constant integral denominator
    - double ax
        lambda coefficient integral numerator
    - double bx
        constant term integral numerator
    - double ax
        lambda coefficient integral numerator
    - double bx
        constant term integral numerator
    - double al
        lower integration boundary
    - double bl
        higher integration boundary
    - int n
        number of interval splits for the Simpson method

  Returns:
  --------
    double integral
        value for the integral
================================================================================
================================================================================
*/
double simpson(const std::function<double(double, double, double, double,
                                          double, double)> &ibsintegrand,
               double ax, double bx, double a, double b, double c, double al,
               double bl, int n) {
  double h, integral, x, sum = 0.0;
  int i;

  h = fabs(bl - al) / n;
  for (i = 0; i < n; i++) {
    x = al + i * h;
    if (i % 2 == 0) {
      sum += 2.0 * ibsintegrand(x, ax, bx, a, b, c);
    } else {
      sum += 4.0 * ibsintegrand(x, ax, bx, a, b, c);
    }
  }

  double fa = ibsintegrand(al, ax, bx, a, b, c);
  double fb = ibsintegrand(bl, ax, bx, a, b, c);

  integral = (h / 3.0) * (fa + fb + sum);

  return integral;
}

/*
================================================================================
================================================================================
STANDARD SIMPSON INTEGRAL FOR BJORKEN-MTINGWA INTEGRAND.
================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF:
    - CERN NOTE CERN-AB-2006-002 EQ 8

================================================================================
  Arguments:
  ----------
    - double BjorkenMtingwaIntegrand(double, double, double, double, double,
                                      double)
    - double ax
        lambda coefficient integral numerator x
    - double bx
        constant coefficient integral numerator x
    - double ay
        lambda coefficient integral numerator y
    - double by
        constant coefficient integral numerator y
    - double as
        lambda coefficient integral numerator s
    - double bs
        constant coefficient integral numerator s
    - double a
        lambda**2 coefficient integral denominator
    - double b
        lambda coefficient integral denominator
    - double c
        constant integral denominator
    - double *integral
        return variable

  Returns:
  --------
    double [3] integral
      IBS growth rates for AMPLITUDES (SIGMA NOT EMIT -> MULTIPLY WITH TWO FOR
      EMIT)
      0 -> al
      1 -> ax
      2 -> ay

================================================================================
================================================================================
*/
void intSimpson(
    const std::function<double(double, double, double, double, double, double)>
        &BjorkenMtingwaIntegrand,
    double ax, double bx, double ay, double by, double as, double bs, double a,
    double b, double c, double *integral) {
  double al[31], bl[30], aloop, bloop;

  int maxdec = 30, ns = 50;
  double test = 1.0e-7;
  bool flag = 0;
  // double integral[3] = {0.0,0.0,0.0};

  // init integral
  integral[0] = 0.0;
  integral[1] = 0.0;
  integral[2] = 0.0;

  al[0] = 0.0;
  for (int iloop = 0; iloop < maxdec; iloop++) {
    bl[iloop] = pow(10.0, iloop);
    al[iloop + 1] = bl[iloop];
    aloop = al[iloop];
    bloop = bl[iloop];

    double nintx =
        simpson(BjorkenMtingwaIntegrand, ax, bx, a, b, c, aloop, bloop, ns);
    double ninty =
        simpson(BjorkenMtingwaIntegrand, ay, by, a, b, c, aloop, bloop, ns);
    double nints =
        simpson(BjorkenMtingwaIntegrand, as, bs, a, b, c, aloop, bloop, ns);

    integral[0] += nints;
    integral[1] += nintx;
    integral[2] += ninty;

    if ((fabs(nints / integral[0]) < test) &
        (fabs(nintx / integral[1]) < test) & (fabs(ninty / integral[2]) < test)

    ) {
      flag = 1;
      break;
    }
  }
  if (flag == 0) {
    integral[0] = 0.0;
    integral[1] = 0.0;
    integral[2] = 0.0;
  }
  // return integral;
}

/*
================================================================================
================================================================================
================================================================================
MODELS BELOW ARE MORE COMPLETE
EQUATIONS ARE MORE COMPLICATED AND
THE CODE MIGHT BE NOT COMPLETELY BUG FREE AT THIS POINT

FOR EQUATIONS CHECK CERN NOTE
CERN-AB-2006-002
================================================================================
================================================================================
================================================================================
*/

/*
================================================================================
================================================================================
SIMPSON DECADE INTEGRAL FOR BJORKEN-MTINGWA INTEGRAND.
================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF:
    - CERN NOTE CERN-AB-2006-002 EQ 8

================================================================================
  Arguments:
  ----------
    - double pnumber
        number of real particles in the bunch
    - double ex
        horizontal emittance
    - double ey
        vertical emittance
    - double sigs
        bunch length
    - double sige
        energy spread dE/E
    - double gammas
        relativistic gamma
    - double betx
        beta x
    - double bety
        beta y
    - double alx
        alfa x
    - double aly
        alfa y
    - double dx
        dispersion x
    - double dpx
        dispersion derivative x
    - double dy
        dispersion y
    - double dpy
        dispersion derivative y

  Returns:
  --------
    double [3] tau
      IBS growth rates for AMPLITUDES (SIGMA NOT EMIT -> MULTIPLY WITH TWO FOR
      EMIT) FOR BJORKEN MTINGWA MODEL
      0 -> al
      1 -> ax
      2 -> ay

================================================================================
================================================================================
*/
void BjorkenMtingwaInt(double pnumber, double ex, double ey, double sigs,
                       double sige, double gammas, double betx, double bety,
                       double alx, double aly, double dx, double dpx, double dy,
                       double dpy, double *tau) {
  const double one = 1.0;
  const double two = 2.0;
  const double three = 3.0;
  const double ten = 10.0;
  const double tstlog = 74.0, power = -two / three;

  double a, b, c, ax, bx, ay, by, al, bl, c1, c3, c2, cx, cy, cyy;
  double cscale, chklog, cprime;
  double phix, tl1, tl2, tx1, tx2, ty1, ty2;
  double cl;

  // temp variable for reducing typing
  phix = dpx + (alx * dx / betx);

  c1 = (gammas * dx) * (gammas * dx) / (ex * betx);
  c3 = betx / ex;
  c2 = c3 * (gammas * phix) * (gammas * phix);
  cx = c1 + c2;
  cy = bety / ey;
  cl = (gammas / sige) * (gammas / sige);

  // a, b and c coefficients of BM integral
  a = cx + cl;
  b = (c3 + cy) * (c1 + cl) + cy * c2;
  c = (c3 * cy) * (c1 + cl);

  // Numerator coefficients
  // horizontal
  ax = two * a;
  bx = b;

  // vertical
  ay = -a;
  by = bx - three * c3 * (c1 + cl);

  // longitudinal
  al = ax;
  bl = b;

  // rescaling
  cscale = one;
  chklog = log10(c3) + log10(cy) + log10(c1 + cl);
  if (chklog > tstlog)
    cscale = pow(ten, (tstlog - chklog));
  cprime = c * cscale;

  //---- Constants for integration loop.
  //     To keep the numbers reasonable, the numerator is
  //     scaled by 1/CPRIME and the denominator by 1/CPRIME**2.
  //     The extra factor of CPRIME is accounted for after integrating
  cyy = pow(cprime, power);

  tl1 = al / cprime;
  tl2 = bl / cprime;

  ty1 = ay / cprime;
  ty2 = by / cprime;

  tx1 = ax / cprime;
  tx2 = bx / cprime;

  SimpsonDecade(a, b, c, cl, cx, cy, cprime, cyy, tl1, tl2, tx1, tx2, ty1, ty2,
                tau);
}

/*
================================================================================
================================================================================
SIMPSON DECADE INTEGRAL FOR CONTE-MARTINI INTEGRAND.
================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF:
    - CERN NOTE CERN-AB-2006-002 EQ 8

================================================================================
  Arguments:
  ----------
    - double pnumber
        number of real particles in the bunch
    - double ex
        horizontal emittance
    - double ey
        vertical emittance
    - double sigs
        bunch length
    - double sige
        energy spread dE/E
    - double gammas
        relativistic gamma
    - double betx
        beta x
    - double bety
        beta y
    - double alx
        alfa x
    - double aly
        alfa y
    - double dx
        dispersion x
    - double dpx
        dispersion derivative x
    - double dy
        dispersion y
    - double dpy
        dispersion derivative y

  Returns:
  --------
    double [3] tau
      IBS growth rates for AMPLITUDES (SIGMA NOT EMIT -> MULTIPLY WITH TWO FOR
      EMIT) FOR CONTE-MARTINI MODEL
      0 -> al
      1 -> ax
      2 -> ay

================================================================================
================================================================================
*/
void ConteMartiniInt(double pnumber, double ex, double ey, double sigs,
                     double sige, double gammas, double betx, double bety,
                     double alx, double aly, double dx, double dpx, double dy,
                     double dpy, double *tau) {

  // const double zero = 0.0;
  const double one = 1.0;
  const double two = 2.0;
  const double three = 3.0;
  // const double four = 4.0;
  // const double six = 6.0;
  // const double eight = 8.0;
  const double ten = 10.0;

  // const double onetominus20 = 1.0e-20;
  const double tstlog = 74.0, power = -two / three;
  // const double coeff[2] = {2.0, 4.0};

  // bool flag = 0;

  double a, b, c, ax, bx, ay, by, al, bl, c1, c3, c2, cx, cy, cyy;
  double cscale, chklog, cprime;
  double phix, tl1, tl2, tx1, tx2, ty1, ty2;
  // double zintl, zintx, zinty, term;
  // double aarr[31], barr[30], h, aloop;
  double cl, hx, hxg2, hxg2obx, hxg2oex;

  double gammas2;

  // temp variable for reducing typing
  phix = dpx + (alx * dx / betx);
  // phiy = dpy + (aly * dy / bety);

  gammas2 = gammas * gammas;
  // phix2 = phix * phix;

  hx = (dx * dx + betx * betx * phix * phix) / betx;
  // hy = (dy * dy + bety * bety * phiy * phiy) / bety;

  c1 = (gammas * dx) * (gammas * dx) / (ex * betx);
  c3 = betx / ex;
  c2 = c3 * (gammas * phix) * (gammas * phix);
  cx = c1 + c2;
  cy = bety / ey;
  cl = (gammas / sige) * (gammas / sige);

  hxg2 = hx * gammas2;

  hxg2obx = hxg2 / betx;
  hxg2oex = hxg2 / ex;

  a = hxg2oex + cl + c3 + cy;
  b = (c3 + cy) * (c1 + cl) + cy * c2 + c3 * cy;
  c = (c3 * cy) * (c1 + cl);

  ax = 2.0 * a - 3.0 * (c3 + cy) - c3 +
       (1.0 / hxg2obx) * (6.0 * c2 - cl + 2.0 * c3 - cy);

  bx = (c3 + cy) * (hxg2oex + cl) - c3 * c2 + c3 * (c3 - 4.0 * cy) +
       (1.0 / hxg2obx) *
           (6.0 * cy * c2 + c3 * cy - gammas * gammas * c3 * c3 * phix * phix +
            cl * (c3 - 2.0 * cy));
  ay = -(gammas2 * hx / ex + cl + c3 - 2.0 * cy);
  by = (c3 + cy) * (c1 + cl) + (cy * c2) - (3.0 * c3 * (c1 + cl)) + (c3 * cy);

  al = 2.0 * (hxg2oex + cl) - c3 - cy;
  bl = (c3 + cy) * (c1 + cl) + (cy * c2) - (2.0 * c3 * cy);

  // rescaling
  cscale = one;
  chklog = log10(c3) + log10(cy) + log10(c1 + cl);
  if (chklog > tstlog)
    cscale = pow(ten, (tstlog - chklog));
  cprime = c * cscale;

  //---- Constants for integration loop.
  //     To keep the numbers reasonable, the numerator is
  //     scaled by 1/CPRIME and the denominator by 1/CPRIME**2.
  //     The extra factor of CPRIME is accounted for after integrating
  cyy = pow(cprime, power);

  tl1 = al / cprime;
  tl2 = bl / cprime;

  ty1 = ay / cprime;
  ty2 = by / cprime;

  tx1 = ax / cprime;
  tx2 = bx / cprime;

  SimpsonDecade(a, b, c, cl, cx, cy, cprime, cyy, tl1, tl2, tx1, tx2, ty1, ty2,
                tau);
}

/*
================================================================================
================================================================================
SIMPSON DECADE INTEGRAL FOR MADX ZIMMERMAN INTEGRAND.
================================================================================
  AUTHORS:
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF:
    - CERN NOTE CERN-AB-2006-002 EQ 8

================================================================================
  Arguments:
  ----------
    - double pnumber
        number of real particles in the bunch
    - double ex
        horizontal emittance
    - double ey
        vertical emittance
    - double sigs
        bunch length
    - double sige
        energy spread dE/E
    - double gammas
        relativistic gamma
    - double betx
        beta x
    - double bety
        beta y
    - double alx
        alfa x
    - double aly
        alfa y
    - double dx
        dispersion x
    - double dpx
        dispersion derivative x
    - double dy
        dispersion y
    - double dpy
        dispersion derivative y

  Returns:
  --------
    double [3] tau
      IBS growth rates for AMPLITUDES (SIGMA NOT EMIT -> MULTIPLY WITH TWO FOR
      EMIT) FOR MADX ZIMMERMAN MODEL
      0 -> al
      1 -> ax
      2 -> ay

================================================================================
================================================================================
*/
void MadxInt(double pnumber, double ex, double ey, double sigs, double sige,
             double gammas, double betx, double bety, double alx, double aly,
             double dx, double dpx, double dy, double dpy, double *tau) {
  // const int maxdec = 30, ns = 50;

  // const double zero = 0.0;
  const double one = 1.0;
  const double two = 2.0;
  const double three = 3.0;
  // const double four = 4.0;
  // const double six = 6.0;
  // const double eight = 8.0;
  const double ten = 10.0;

  // const double onetominus20 = 1.0e-20;
  const double tstlog = 74.0, power = -two / three;
  // const double coeff[2] = {2.0, 4.0};

  // bool flag = 0;

  double a, b, c, ax, bx, ay, by, al, bl, c1, c3, c2, cx, cy, cyy;
  double cscale, chklog, cprime;
  double phix, phiy, tl1, tl2, tx1, tx2, ty1, ty2;
  double c1y, cl, hx, hy, hxg2, hyg2, hxg2obx, hxg2oex, hyg2oby, hyg2oey, c2y;

  // temp variable for reducing typing
  phix = dpx + (alx * dx / betx);
  phiy = dpy + (aly * dy / bety);

  hx = (dx * dx + betx * betx * phix * phix) / betx;
  hy = (dy * dy + bety * bety * phiy * phiy) / bety;

  c1 = (gammas * dx) * (gammas * dx) / (ex * betx);
  c1y = (gammas * dy) * (gammas * dy) / (ey * bety);
  c3 = betx / ex;
  c2 = c3 * (gammas * phix) * (gammas * phix);
  cx = c1 + c2;
  cy = bety / ey;
  cl = (gammas / sige) * (gammas / sige);
  c2y = cy * (gammas * phiy) * (gammas * phiy);

  hxg2 = hx * gammas * gammas;
  hyg2 = hy * gammas * gammas;

  hxg2obx = hxg2 / betx;
  hxg2oex = hxg2 / ex;

  hyg2oby = hyg2 / bety;
  hyg2oey = hyg2 / ey;

  a = hxg2oex + hyg2oey + cl + c3 + cy;
  b = (c3 + cy) * (c1 + cl + c1y) + cy * c2 + c3 * c2y + c3 * cy;
  c = c3 * cy * (c1 + cl + c1y);

  al = (2.0 * a - 3.0 * cy - 3.0 * c3);
  bl = (c3 + cy) * (hxg2oex + hyg2oey + cl) - 2.0 * c3 * cy -
       (c3 * c2 + cy * c2y);

  ax = 2.0 * (hxg2oex + hyg2oey + cl) - (hyg2oey / hxg2obx) +
       (1.0 / hxg2obx) * (2.0 * c3 - cy - cl);
  bx = (c3 + cy) * (hxg2oex + hyg2oey + cl) - (c3 * c2 + cy * c2y) +
       c3 * (c3 - 4.0 * cy) +
       (1.0 / hxg2obx) *
           (c3 * cy - gammas * gammas * c3 * c3 * phix * phix +
            2.0 * gammas * gammas * cy * cy * phiy * phiy +
            cl * (c3 - 2.0 * cy) + (hyg2oey / hxg2obx) * (c3 - 2.0 * cy));

  ay = -(hxg2oex + 2.0 * hyg2oey + cl + hxg2oex * betx / bety) +
       2.0 * hyg2oby * (hxg2oex + hyg2oey) + 2.0 * hyg2oby * cl -
       (c3 - 2.0 * cy);
  by = (cy - 2.0 * c3) * (hxg2oex + cl) + (cy - 4.0 * c3) * hyg2oey + c3 * cy +
       2.0 * c3 * c2 - cy * c2y + hyg2oby * (c3 + cy) * (hyg2oey + cl) +
       (c3 + cy) * (hxg2oex * hyg2oby) - hyg2oby * (c3 * c2 + cy * c2y);

  // rescaling
  cscale = one;
  chklog = log10(c3) + log10(cy) + log10(c1 + cl);
  if (chklog > tstlog)
    cscale = pow(ten, (tstlog - chklog));
  cprime = c3 * cy * cscale * (c1 + cl + c1y);

  //---- Constants for integration loop.
  //     To keep the numbers reasonable, the numerator is
  //     scaled by 1/CPRIME and the denominator by 1/CPRIME**2.
  //     The extra factor of CPRIME is accounted for after integrating
  cyy = pow(cprime, power);

  tl1 = al / cprime;
  tl2 = bl / cprime;

  ty1 = ay / cprime;
  ty2 = by / cprime;

  tx1 = ax / cprime;
  tx2 = bx / cprime;

  SimpsonDecade(a, b, c, cl, cx, cy, cprime, cyy, tl1, tl2, tx1, tx2, ty1, ty2,
                tau);
}

/*
================================================================================
================================================================================
SIMPSON INTEGRATOR SPLITTING THE INTEGRAL IN DECADES FOR FASTER CONVERGENCE
WITH SCALING.

DETAILS:
        MAX DECADES = 30 (10**30)
        INTERVAL SPLITS (NS) = 50

REF:
        MADX ORIGINAL SOURCE CODE IN TWSINT FUNCTION
        CERN NOTE CERN-AB-2006-002 EQ 8
================================================================================
  AUTHORS:
    - MADX ORIGINAL AUTHORS
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF:
    - CERN NOTE CERN-AB-2006-002 EQ 8

================================================================================
  Arguments:
  ----------
    - double pnumber
        number of real particles in the bunch
    - double ex
        horizontal emittance
    - double ey
        vertical emittance
    - double sigs
        bunch length
    - double sige
        energy spread dE/E
    - double gammas
        relativistic gamma
    - double betx
        beta x
    - double bety
        beta y
    - double alx
        alfa x
    - double aly
        alfa y
    - double dx
        dispersion x
    - double dpx
        dispersion derivative x
    - double dy
        dispersion y
    - double dpy
        dispersion derivative y

  Returns:
  --------
    double [3] tau
      IBS growth rates for AMPLITUDES (SIGMA NOT EMIT -> MULTIPLY WITH TWO FOR
      EMIT) FOR MADX ZIMMERMAN MODEL
      0 -> al
      1 -> ax
      2 -> ay

================================================================================
================================================================================
*/

void twsint(double pnumber, double ex, double ey, double sigs, double sige,
            double gammas, double betax, double betay, double alx, double aly,
            double dx, double dpx, double dy, double dpy, double *tau) {

  // int iiz, iloop;
  int maxdec = 30, ns = 50;

  double a, b, am, c1, c2, c3, c1y, c2y, chy, cx, cy, cl, r1;
  double cscale, chklog, cprime, ccy;
  double zintl, zintx, zinty;
  double td1, td2, tl1, tl2, tx1, tx2, ty1, ty2;
  double al[31], bl[30], h, aloop;
  double term, func, polyl, polyx, polyy, suml, sumx, sumy;
  double cof, f, alam, phi, phiy, tmpl, tmpx, tmpy;

  double zero = 0.0;
  double one = 1.0;
  double two = 2.0;
  double three = 3.0;
  double four = 4.0;
  double six = 6.0;
  double ten = 10.0;
  double onetominus20 = 1e-20;

  // parameters
  double tstlog = 74.0, power = -two / three, test = 1e-7;
  double coeff[2] = {2.0, 4.0};

  // squares to redude typing
  double gammas2 = gammas * gammas;
  double sige2 = sige * sige;

  // switch for jumping out of the
  // integral loop
  bool breakflag = 0;

  phi = dpx + (alx * dx / betax);
  phiy = dpy + (aly * dy / betay);
  am = one;
  c1 = (gammas * dx) * (gammas * dx) / (ex * betax);
  c1y = (gammas * dy) * (gammas * dy) / (ey * betay);
  c3 = betax / ex;
  c2 = c3 * (gammas * phi) * (gammas * phi);
  cx = c1 + c2;
  cl = am * (gammas / sige) * (gammas / sige);
  cy = betay / ey;
  c2y = cy * (gammas * phiy) * (gammas * phiy);
  chy = c1y + c2y;
  r1 = three / cy;
  a = cx + cl + chy + c3 + cy;
  b = (c3 + cy) * (c1 + cl + c1y) + cy * c2 + c3 * c2y + c3 * cy;

  // rescaling
  cscale = one;
  chklog = log10(c3) + log10(cy) + log10(c1 + cl);
  if (chklog > tstlog)
    cscale = pow(ten, (tstlog - chklog));
  cprime = c3 * cy * cscale * (c1 + cl + c1y);

  // printf("cscale : %f", cscale);
  // printf("chklog : %f", chklog);
  // printf("cprime : %f", cprime);

  //---- Split integral into decades, with NS steps per decade.
  //     variables to save integral segments
  zintl = zero;
  zintx = zero;
  zinty = zero;

  //---- Constants for integration loop.
  //     To keep the numbers reasonable, the numerator is
  //     scaled by 1/CPRIME and the denominator by 1/CPRIME**2.
  //     The extra factor of CPRIME is accounted for after integrating
  ccy = pow(cprime, power);
  td1 = (a - cy) * ccy;
  td2 = one / (sqrt(ccy) * cscale * cy);
  tl1 = (two * a - three * cy - three * c3) / cprime;
  tl2 = (b - three * c3 * cy) / cprime;

  // --- corrected ty1 23.02.2011
  ty1 = (-a + three * cy - chy - chy / cy * (c3 - two * gammas2 / sige2) +
         two * chy * (cx + chy) / cy + six * c2y) /
        cprime;

  // --- corrected ty2 23.02.2011
  ty2 =
      (b - c1y * (c3 + cy) + chy * (cy + chy) +
       chy * ey * (one / ey + betax / (betay * ex)) * gammas2 / sige2 -
       chy * betax / ex * four +
       (one + (betax * ey) / (betay * ex)) * cx * chy +
       (chy * chy) * (betax * ey) / (betay * ex) - chy * ey * c2 * c3 / betay -
       c2y * (cy + c3 + chy) + three * c3 * (two * c2y + c1y)) /
          cprime -
      r1 / cscale;

  //--- corrected tx1 23.02.2011
  tx1 = (two * (a - c3 - cy) * (cx - c3) - cy * cx +
         c3 * (c1y + six * c2 + c2y + two * c3 + cl - cy)) /
        cprime;

  // --- corrected tx2 23.02.2011
  tx2 = (c3 + cx) * ((b - c1y * (c3 + cy)) / cprime) - six / cscale +
        three * c3 * cy * (cl / cprime) +
        (six * c3 * cy * c1y + (betay / ey + betax / ex) * chy * cx +
         chy * (c3 * c3 - two * cy * c3) - c2y * cx * (cy + c3) +
         (two * cy * c3 - c3 * c3) * c2y) /
            cprime;

  // init
  al[0] = zero;

  for (int iloop = 0; iloop < maxdec; iloop++) {
    bl[iloop] = pow(ten, iloop);
    al[iloop + 1] = bl[iloop];
    h = (bl[iloop] - al[iloop]) / ns;
    aloop = al[iloop];

    //---- Evaluate Simpson's rule summation for one interval.
    //     The integrand is calculated in the loop itself
    if (fabs(cy + aloop) > onetominus20) {
      term = sqrt((cy + aloop) * ccy) *
             sqrt((aloop * ccy * aloop + td1 * aloop + td2) +
                  aloop * c2y * (c3 - cy) * ccy / (cy + aloop));
    } else {
      term = sqrt((cy + aloop) * ccy) *
             sqrt((aloop * ccy * aloop + td1 * aloop + td2));
    }

    func = sqrt(aloop) / (term * term * term);
    polyl = tl1 * aloop + tl2;
    polyx = tx1 * aloop + tx2;
    polyy = ty1 * aloop + ty2;
    suml = func * polyl;
    sumx = func * polyx;
    sumy = func * polyy;

    for (int iiz = 0; iiz < ns; iiz++) {
      alam = aloop + iiz * h;
      cof = coeff[iiz % 2];
      if (fabs(cy + alam) > onetominus20) {
        term = sqrt((cy + alam) * ccy) *
               sqrt((alam * ccy * alam + td1 * alam + td2) +
                    alam * c2y * (c3 - cy) * ccy / (cy + alam));
      } else {
        term = sqrt((cy + alam) * ccy) *
               sqrt((alam * ccy * alam + td1 * alam + td2));
      }
      f = sqrt(alam) / (term * term * term);
      polyl = tl1 * alam + tl2;
      polyx = tx1 * alam + tx2;
      polyy = ty1 * alam + ty2;

      suml = suml + cof * f * polyl;
      sumx = sumx + cof * f * polyx;
      sumy = sumy + cof * f * polyy;
    }

    suml = suml - f * polyl;
    sumx = sumx - f * polyx;
    sumy = sumy - f * polyy;
    tmpl = (suml / three) * h;
    tmpx = (sumx / three) * h;
    tmpy = (sumy / three) * h;
    zintl = zintl + tmpl;
    zintx = zintx + tmpx;
    zinty = zinty + tmpy;

    if ((fabs(tmpl / zintl) < test) & (fabs(tmpx / zintx) < test) &
        (fabs(tmpy / zinty) < test)) {
      // ---- Divide answers by cprime to account for scaling.
      tau[1] = (zintx / cprime);
      tau[0] = cl * (zintl / cprime);
      tau[2] = cy * (zinty / cprime);

      breakflag = 1;
      break;
    }
    // printf("%.15e %.15e %.15e\n", tmpl/zintl, tmpx/zintx, tmpy/zinty);
  }
  if (breakflag == 0) {
    tau[0] = 0.0;
    tau[1] = 0.0;
    tau[2] = 0.0;
  }
}