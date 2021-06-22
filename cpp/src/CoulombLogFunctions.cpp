#include "../include/ibs_bits/NumericFunctions.hpp"
#include "../include/ibs_bits/RadiationDamping.hpp"
#include <math.h>
#include <stdio.h>

/*
================================================================================
================================================================================
METHOD TO CALCULATE COULOMBLOG PER ELEMENT OR USING ELEMENT BY ELEMENT
TWISS DATA.

================================================================================
  AUTHORS:
    - MADX COPYRIGHT CERN
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF:
        Calculation of Coulomb logarithm (and print)
        based on the formulae in AIP physics vade mecum p.264 (1981)

================================================================================
  Arguments:
  ----------
    - double pnumber
        number of particles
    - double bx
        beta x
    - double by
        beta y
    - double dx
        dispersion x
    - double dy
        dispersion y
    - double ex
        emittance x
    - double ey
        emittance y
    - double r0
        classical particle radius
    - double gamma
        relativistic gamma
    - double charge
        particle charge
    - double en0
        particle energy
    - double amass
        particle mass
    - double sige
        energy spread
    - double sigt
        bunch length (s)
    - double* clog
        output variable

  Returns:
  --------
  double[2] clog
    1 -> coulomblog
    2 -> coulomblog constant

================================================================================
================================================================================
*/
void twclog(double pnumber, double bx, double by, double dx, double dy,
            double ex, double ey, double r0, double gamma, double charge,
            double en0, double amass, double sige, double sigt, double *clog) {
  // constants
  const double ot2 = 1.0e2;
  const double ft8 = 5.0e8;
  const double ot5 = 1.0e5;
  const double ttm3 = 2.0e-3;
  const double fac1 = 743.4;
  const double fac2 = 1.44e-7;
  const double c = clight;

  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // IMPORTANT: HBAR USED HERE IS REDUCED PLANCK CONSTANT IN GEV!!!!
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // const double hbar = 6.582119569e-25; // 1.0545718176461565e-34;

  const double two = 2.0;
  const double eight = 8.0f;

  double betar = sqrt(1.0 - 1.0 / (gamma * gamma));
  double qion = fabs(charge);
  double etrans = ft8 * (gamma * en0 - amass) * (ex / bx);
  double tempev = two * etrans;

  //---- Calculate beam volume to get density (in cm**-3).
  double sigxcm = ot2 * sqrt(ex * bx + (dx * dx * sige * sige));
  double sigycm = ot2 * sqrt(ey * by + (dy * dy * sige * sige));
  double sigtcm = ot2 * sigt;

  double vol = eight * sqrt(pi * pi * pi) * sigxcm * sigycm * sigtcm;
  double densty = pnumber / vol;

  //---- Calculate RMAX as smaller of SIGXCM and DEBYE length.
  double debyel = fac1 * sqrt(tempev / densty) / qion;
  double rmax = min(sigxcm, debyel);

  //---- Calculate RMIN as larger of classical distance of closest approach
  //     or quantum mechanical diffraction limit from nuclear radius.
  double rmincl = fac2 * qion * qion / tempev;
  double rminqm = hbar * c * ot5 / (two * sqrt(ttm3 * etrans * amass));
  double rmin = max(rmincl, rminqm);
  double coulog = log(rmax / rmin);
  // calculate coulomb log constant pre-factor
  double constt = pnumber * coulog * r0 * r0 * c /
                  (eight * pi * betar * betar * betar * gamma * gamma * gamma *
                   gamma * ex * ey * sige * sigt);

  clog[0] = coulog;
  clog[1] = constt;
}
/*
================================================================================
================================================================================
METHOD TO CALCULATE COULOMBLOG PER ELEMENT OR USING ELEMENT BY ELEMENT TWISS
DATA WITH TAILCUT.

================================================================================
  AUTHORS:
    - MADX COPYRIGHT CERN
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF:
        Calculation of Coulomb logarithm (and print)
        based on the formulae in AIP physics vade mecum p.264 (1981)

================================================================================
  Arguments:
  ----------
    - double pnumber
        number of particles
    - double bx
        beta x
    - double by
        beta y
    - double dx
        dispersion x
    - double dpx
        dispersion derivative x
    - double dy
        dispersion y
    - double dpy
        dispersion derivative y
    - double ax
        alfa x
    - double ay
        alfa y
    - double angle
        bending angle
    - double k1l
        k1l from twiss madx
    - double k1sl
        k1sl from twiss madx
    - double ex
        emittance x
    - double ey
        emittance y
    - double r0
        classical particle radius
    - double aatom
        atomic number A
    - double gamma
        relativistic gamma
    - double charge
        particle charge
    - double en0
        particle energy
    - double len
        acceleratro length
    - double amass
        particle mass
    - double sige
        energy spread
    - double sigt
        bunch length (s)
    - double* clog
        output variable

  Returns:
  --------
  double[2] clog
    1 -> coulomblog
    2 -> coulomblog constant

================================================================================
================================================================================
*/
void twclogtail(double pnumber, double l, double bx, double by, double dx,
                double dpx, double dy, double dpy, double ax, double ay,
                double angle, double k1l, double k1sl, double ex, double ey,
                double r0, double aatom, double gamma, double en0, double len,
                double amass, double charge, double sige, double sigt,
                double *clog) {
  // constants
  const double ot2 = 1.0e2;
  const double ft8 = 5.0e8;
  const double ot5 = 1.0e5;
  const double ttm3 = 2.0e-3;
  const double fac1 = 743.4;
  const double fac2 = 1.44e-7;
  const double c = clight;
  // const double hbar = 1.0545718176461565e-34;
  const double electron_volt_joule_relationship = 1.602176634e-19;

  // For tailcut we need local radiation integrals
  double *radint = RadiationDampingElement(l, bx, by, dx, dpx, dy, dpy, ax, ay,
                                           angle, k1l, k1sl);

  double p0 = en0 * 1.0e9;
  double restE = amass * 1.0e9;
  double particle_radius = charge * charge / aatom * 1.54e-18;
  double CalphaEC = particle_radius * c / (3.0 * restE * restE * restE) *
                    (p0 * p0 * p0 / len);

  // transverse partition numbers
  double jx = 1.0 - radint[3] / radint[0];
  double jy = 1.0 - radint[4] / radint[0];
  double alphax = 2.0 * CalphaEC * radint[0] * jx;
  double alphay = 2.0 * CalphaEC * radint[0] * jy;
  double alphas = 2.0 * CalphaEC * radint[0] * (jx + jy);
  double tauradxbunch = (1.0 / alphax) / gamma;
  double tauradybunch = (1.0 / alphay) / gamma;
  double tauradsbunch = (1.0 / alphas) / gamma;
  double tauradmax = max(max(tauradxbunch, tauradybunch), tauradsbunch);
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // IMPORTANT: HBAR USED HERE IS REDUCED PLANCK CONSTANT IN GEV!!!!
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // const double hbar = 6.582119569e-25; // 1.0545718176461565e-34;

  const double two = 2.0;
  const double eight = 8.0;

  double betar = sqrt(1.0 - 1.0 / (gamma * gamma));
  double qion = fabs(charge);
  double etrans = ft8 * (gamma * en0 - amass) * (ex / bx);
  double tempev = two * etrans;

  //---- Calculate beam volume to get density (in cm**-3).
  double sigxcm = ot2 * sqrt(ex * bx + (dx * dx * sige * sige));
  double sigycm = ot2 * sqrt(ey * by + (dy * dy * sige * sige));
  double sigtcm = ot2 * sigt;

  double vol = eight * sqrt(pi * pi * pi) * sigxcm * sigycm * sigtcm;
  double densty = pnumber / vol;

  //---- Calculate RMAX as smaller of SIGXCM and DEBYE length.
  double debyel = fac1 * sqrt(tempev / densty) / qion;
  double rmax = min(sigxcm, debyel);

  //---- Calculate RMIN as larger of classical distance of closest approach
  //     or quantum mechanical diffraction limit from nuclear radius.
  double rmincl = fac2 * qion * qion / tempev;
  double rminqm = hbar * c * ot5 / (two * sqrt(ttm3 * etrans * amass));
  double rmin;
  if (tauradmax > 0.0) {
    double rmintailcut =
        1.0f / sqrt(pnumber * pi * tauradmax * c * gamma * sqrt(ex / bx));
    rmin = max(max(rmincl, rminqm), rmintailcut);
    // printf("rmax %12.6e   rmin %12.6e rmincl %12.6e rminqm %12.6e rmintailcut
    // "
    //       "%12.6e\n",
    //       rmax, rmin, rmincl, rminqm, rmintailcut);
  } else {
    rmin = max(rmincl, rminqm);
    // printf("rmax %12.6e   rmin %12.6e rmincl %12.6e rminqm %12.6e\n", rmax,
    //       rmin, rmincl, rminqm);
  }
  double coulog = log(rmax / rmin);

  // calculate coulomb log constant pre-factor
  double constt = pnumber * coulog * r0 * r0 * c /
                  (eight * pi * betar * betar * betar * gamma * gamma * gamma *
                   gamma * ex * ey * sige * sigt);

  clog[0] = coulog;
  clog[1] = constt;
}
/*
================================================================================
================================================================================
METHOD TO CALCULATE COULOMBLOG PER ELEMENT OR USING RING AVERAGED VALUES

================================================================================
  AUTHORS:
    - MADX COPYRIGHT CERN
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF:
        Calculation of Coulomb logarithm (and print)
        based on the formulae in AIP physics vade mecum p.264 (1981)

================================================================================
  Arguments:
  ----------
    - double pnumber
        number of particles
    - double ex
        emittance x
    - double ey
        emittance y
    - map<string, double>& twissheader
        madx twiss header
    - double sige
        energy spread
    - double sigt
        bunch length (s)
    - double r0
        classical particle radius
    - bool printout
        printout summary
    - double* clog
        output variable

  Returns:
  --------
  double[2] clog
    1 -> coulomblog
    2 -> coulomblog constant

================================================================================
================================================================================
*/
void CoulombLog(double pnumber, double ex, double ey,
                map<string, double> &twissheader, double sige, double sigt,
                double r0, bool printout, double *clog) {
  // const double one   = 1.0f;
  const double two = 2.0;
  // const double four  = 4.0f;
  const double eight = 8.0f;
  const double ot2 = 1.0e2;
  const double ft8 = 5.0e8;
  const double ot5 = 1.0e5;
  const double ttm3 = 2.0e-3;
  const double fac1 = 743.4;
  const double fac2 = 1.44e-7;

  const double c = clight;
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // IMPORTANT: HBAR USED HERE IS REDUCED PLANCK CONSTANT IN GEV!!!!
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // static double output[2];
  double gamma = twissheader["GAMMA"];
  double charge = twissheader["CHARGE"];
  double len = twissheader["LENGTH"];
  double en0 = twissheader["ENERGY"];
  double amass = twissheader["MASS"];
  double q1 = twissheader["Q1"];
  double q2 = twissheader["Q2"];
  double dxbar = twissheader["DXRMS"];
  double dybar = twissheader["DYRMS"];

  // necessary parameters
  // double r0 = charge * charge / aatom * 1.54e-18;
  double betar = sqrt(1 - 1 / (gamma * gamma));
  double bxbar = len / (2.0 * pi * q1); // avg betax
  double bybar = len / (2.0 * pi * q2); // avgbety

  //---- Calculate transverse temperature as 2*P*X',
  //     i.e., assume the transverse energy is temperature/2.
  double qion = fabs(charge);
  double etrans = ft8 * (gamma * en0 - amass) * (ex / bxbar);
  double tempev = two * etrans;

  //---- Calculate beam volume to get density (in cm**-3).
  double sigxcm = ot2 * sqrt(ex * bxbar + (dxbar * dxbar * sige * sige));
  double sigycm = ot2 * sqrt(ey * bybar + (dybar * dybar * sige * sige));
  double sigtcm = ot2 * sigt;

  double vol = eight * sqrt(pi * pi * pi) * sigxcm * sigycm * sigtcm;
  double densty = pnumber / vol;

  //---- Calculate RMAX as smaller of SIGXCM and DEBYE length.
  double debyel = fac1 * sqrt(tempev / densty) / qion;
  double rmax = min(sigxcm, debyel);

  //---- Calculate RMIN as larger of classical distance of closest approach
  //     or quantum mechanical diffraction limit from nuclear radius.
  double rmincl = fac2 * qion * qion / tempev;
  double rminqm = hbar * c * ot5 / (two * sqrt(ttm3 * etrans * amass));
  double rmin = max(rmincl, rminqm);
  double coulog = log(rmax / rmin);

  double constt = pnumber * coulog * r0 * r0 * c /
                  (eight * pi * betar * betar * betar * gamma * gamma * gamma *
                   gamma * ex * ey * sige * sigt);
  double cbunch = qion * pnumber * ec * betar * c / len;

  if (printout) {
    printf("Radius max          : %.8e\n", rmax);
    printf("Density             : %.8e\n", pow(densty, -1.0 / 3.0));
    printf("CONST               ; %.8e\n", constt);
    printf("ENERGY              ; %.8e GeV\n", en0);
    printf("BETA                ; %.8e\n", betar);
    printf("GAMMA               ; %.8e\n", gamma);
    printf("COULOMB LOG         ; %.8e\n", coulog);
    printf("X-emittance         ; %.8e m*rad\n", ex);
    printf("Y-emittance         ; %.8e m*rad\n", ey);
    printf("Momentum spread     ; %.8e\n", sige);
    printf("Bunch length        ; %.8e m \n", sigt);
    printf("Particles per bunch ; %.8e\n", pnumber);
    printf("Bunch current       ; %.8e A\n", cbunch);
  }

  clog[0] = coulog;
  clog[1] = constt;
}
/*
================================================================================
================================================================================
METHOD TO CALCULATE COULOMBLOG PER ELEMENT OR USING RING AVERAGED VALUES
WITH TAILCUT.

================================================================================
  AUTHORS:
    - MADX COPYRIGHT CERN
    - TOM MERTENS

  HISTORY:
    - 08/06/2021 : initial cpp version (Tom)

  REF:
        Calculation of Coulomb logarithm (and print)
        based on the formulae in AIP physics vade mecum p.264 (1981)

================================================================================
  Arguments:
  ----------
    - double pnumber
        number of particles
    - double ex
        emittance x
    - double ey
        emittance y
    - map<string, double>& twissheader
        madx twiss header
    - double sige
        energy spread
    - double sigt
        bunch length (s)
    - double tauradx
        radiation lifetime x
    - double taurady
        radiation lifetime y
    - double taurads
        radiation lifetime s
    - double r0
        classical particle radius
    - bool printout
        printout summary
    - double* clog
        output variable

  Returns:
  --------
  double[2] clog
    1 -> coulomblog
    2 -> coulomblog constant

================================================================================
================================================================================
*/
void TailCutCoulombLog(double pnumber, double ex, double ey,
                       map<string, double> &twissheader, double sige,
                       double sigt, double tauradx, double taurady,
                       double taurads, double r0, bool printout, double *clog) {
  // const double one   = 1.0f;
  const double two = 2.0f;
  // const double four  = 4.0f;
  const double eight = 8.0f;
  const double ot2 = 1.0e2;
  const double ft8 = 5.0e8;
  const double ot5 = 1.0e5;
  const double ttm3 = 2.0e-3;
  const double fac1 = 743.4;
  const double fac2 = 1.44e-7;

  const double c = clight;
  // const double pi = 3.141592653589793f;
  // const double ec = 1.602176634e-19;
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // IMPORTANT: HBAR USED HERE IS REDUCED PLANCK CONSTANT IN GEV!!!!
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // const double hbar = 6.582119569e-25; // 1.0545718176461565e-34;

  // static double output[2];

  double gamma = twissheader["GAMMA"];
  double charge = twissheader["CHARGE"];
  double len = twissheader["LENGTH"];
  double en0 = twissheader["ENERGY"];
  double amass = twissheader["MASS"];
  double q1 = twissheader["Q1"];
  double q2 = twissheader["Q2"];
  double dxbar = twissheader["DXRMS"];
  double dybar = twissheader["DYRMS"];
  // double gammatr = twissheader[9];

  // necessary parameters
  // double r0 = charge * charge / aatom * 1.54e-18;
  double betar = sqrt(1 - 1 / (gamma * gamma));
  double bxbar = len / (2.0 * pi * q1); // avg betax
  double bybar = len / (2.0 * pi * q2); // avgbety

  double tauradxbunch = tauradx / gamma;
  double tauradybunch = taurady / gamma;
  double tauradsbunch = taurads / gamma;
  double tauradmax = max(max(tauradxbunch, tauradybunch), tauradsbunch);

  //---- Calculate transverse temperature as 2*P*X',
  //     i.e., assume the transverse energy is temperature/2.
  double qion = fabs(charge);
  double etrans = ft8 * (gamma * en0 - amass) * (ex / bxbar);
  double tempev = two * etrans;

  //---- Calculate beam volume to get density (in cm**-3).
  double sigxcm = ot2 * sqrt(ex * bxbar + (dxbar * dxbar * sige * sige));
  double sigycm = ot2 * sqrt(ey * bybar + (dybar * dybar * sige * sige));
  double sigtcm = ot2 * sigt;

  double vol = eight * sqrt(pi * pi * pi) * sigxcm * sigycm * sigtcm;
  double densty = pnumber / vol;

  //---- Calculate RMAX as smaller of SIGXCM and DEBYE length.
  double debyel = fac1 * sqrt(tempev / densty) / qion;
  double rmax = min(sigxcm, debyel);

  //---- Calculate RMIN as larger of classical distance of closest approach
  //     or quantum mechanical diffraction limit from nuclear radius.
  double rmincl = fac2 * qion * qion / tempev;
  double rminqm = hbar * c * ot5 / (two * sqrt(ttm3 * etrans * amass));
  double rmintailcut =
      1.0 / sqrt(pnumber * pi * tauradmax * c * gamma * sqrt(ex / bxbar));
  double rmin = max(max(rmincl, rminqm), rmintailcut);
  double coulog = log(rmax / rmin);

  double constt = pnumber * coulog * r0 * r0 * c /
                  (eight * pi * betar * betar * betar * gamma * gamma * gamma *
                   gamma * ex * ey * sige * sigt);
  double cbunch = qion * pnumber * ec * betar * c / len;

  if (printout) {
    printf("Radius max          : %.8e\n", rmax);
    printf("Density             : %.8e\n", pow(densty, -1.0f / 3.0f));
    printf("CONST               ; %.8e\n", constt);
    printf("ENERGY              ; %.8e GeV\n", en0);
    printf("BETA                ; %.8e\n", betar);
    printf("GAMMA               ; %.8e\n", gamma);
    printf("COULOMB LOG         ; %.8e\n", coulog);
    printf("X-emittance         ; %.8e m*rad\n", ex);
    printf("Y-emittance         ; %.8e m*rad\n", ey);
    printf("Momentum spread     ; %.8e\n", sige);
    printf("Bunch length        ; %.8e m \n", sigt);
    printf("Particles per bunch ; %.8e\n", pnumber);
    printf("Bunch current       ; %.8e A\n", cbunch);
  }

  clog[0] = coulog;
  clog[1] = constt;
}