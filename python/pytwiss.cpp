//#include "../cpp/include/ibs_bits/CoulombLogFunctions.hpp"
//#include "../cpp/include/ibs_bits/NumericFunctions.hpp"
//#include "../cpp/include/ibs_bits/RadiationDamping.hpp"
//#include "../cpp/include/ibs_bits/twiss.hpp"
#include <ibs>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// PYBIND11_MAKE_OPAQUE(std::vector<double>);
namespace py = pybind11;

PYBIND11_MODULE(IBSLib, m) {
  m.doc() = "IBS library";
  /*
================================================================================
                                    TWISS
================================================================================
  */
  m.def("GetTwissHeader", &GetTwissHeader,
        "Get the twiss header as a dictionary.", py::arg("filename"));

  m.def("GetTwissTable", &GetTwissTableAsMap, "Get the twiss data table.",
        py::arg("filename"));

  m.def("updateTwiss",
        [](map<string, vector<double>> &table) {
          updateTwiss(table);
          return table;
        },
        "Extend Twiss Table with rad int, CS gamma, curly H and rho.",
        py::arg("table"));

  m.def("printTwissColumn", &printTwissMap, "Print Twiss column",
        py::arg("columnName"), py::arg("twissTableMap"));
  /*
================================================================================
                                    CONSTANTS
================================================================================
  */
  m.attr("clight") = py::float_(clight);
  m.attr("hbarGeV") = py::float_(hbar);
  m.attr("electron_mass") = py::float_(emass);
  m.attr("proton_mass") = py::float_(pmass);
  m.attr("neutron_mass") = py::float_(nmass);
  m.attr("mu_mass") = py::float_(mumass);
  m.attr("atomic_mass_unit") = py::float_(atomicmassunit);
  m.attr("pi") = py::float_(pi);
  m.attr("electric_charge") = py::float_(ec);
  m.attr("euler") = py::float_(euler);
  m.attr("electron_radius") = py::float_(erad);
  m.attr("proton_radius") = py::float_(prad);
  /*
================================================================================
                             RADIATION DAMPING
================================================================================
 */
  m.def("RadiationDampingApprox",
        [](double latticeLength, double gamma, double gammaTransition,
           double dipoleBendingRadius, double betax, double betay,
           py::array_t<double> out) {
          double *radint;
          radint = RadiationDampingApprox(latticeLength, gamma, gammaTransition,
                                          dipoleBendingRadius, betax, betay);
          // convert double pointer to numpy array
          auto buf_out = out.request();
          double *ptr_out = static_cast<double *>(buf_out.ptr);

          ptr_out[0] = radint[0];
          ptr_out[1] = radint[1];
          ptr_out[2] = radint[2];
          ptr_out[3] = radint[3];
          ptr_out[4] = radint[4];
          ptr_out[5] = radint[5];
          ptr_out[6] = radint[6];
        },
        "Radiation Damping using Ring averages.", py::arg("latticelength"),
        py::arg("gamma"), py::arg("gammaTransition"), py::arg("bendingradius"),
        py::arg("betax"), py::arg("betay"), py::arg("outputArray"));

  m.def("RadiationDampingLattice",
        [](map<string, vector<double>> &table, py::array_t<double> out) {
          double *radint;
          radint = RadiationDampingLattice(table);

          auto buf_out = out.request();
          double *ptr_out = static_cast<double *>(buf_out.ptr);
          ptr_out[0] = radint[0];
          ptr_out[1] = radint[1];
          ptr_out[2] = radint[2];
          ptr_out[3] = radint[3];
          ptr_out[4] = radint[4];
          ptr_out[5] = radint[5];
          ptr_out[6] = radint[6];
        },
        "Radiation damping weighted per element.", py::arg("twissTableMap"),
        py::arg("outputArray"));

  m.def(
      "RadiationDampingEquilibria",
      [](map<string, double> &twissheadermap, py::array_t<double> radint,
         double aatom, double qs, py::array_t<double> out) {
        double *equi;
        double r[7];
        auto buf = radint.request();
        double *ptr = static_cast<double *>(buf.ptr);
        for (int i = 0; i < 7; i++) {
          r[i] = ptr[i];
        }
        equi =
            RadiationDampingLifeTimesAndEquilibriumEmittancesWithPartitionNumbers(
                twissheadermap, r, aatom, qs);
        auto buf_out = out.request();
        double *ptr_out = static_cast<double *>(buf_out.ptr);
        ptr_out[0] = equi[0];
        ptr_out[1] = equi[1];
        ptr_out[2] = equi[2];
        ptr_out[3] = equi[3];
        ptr_out[4] = equi[4];
        ptr_out[5] = equi[5];
        ptr_out[6] = equi[6];
        ptr_out[7] = equi[7];
        ptr_out[8] = equi[8];
      },
      "Get radiation damping equilibs.", py::arg("twissHeaderMap"),
      py::arg("RadIntArray"), py::arg("AtomicMassNumber"),
      py::arg("synchrotronTune"), py::arg("outputArray"));

  m.def("get_energy_loss_per_turn", &RadiationLossesPerTurn,
        "Get the energy loss per turn.", py::arg("twissTableMap"),
        py::arg("I2"), py::arg("AtomicMassNumber"));

  /*
 ================================================================================
                              NUMERIC FUNCTIONS
 ================================================================================
  */
  m.def("beta_relativistic_from_gamma", &BetaRelativisticFromGamma,
        "Beta relativistic.", py::arg("Gamma"));

  m.def("eta", &eta, "Slip factor", py::arg("Gamma"),
        py::arg("GammaTransition"));

  m.def("particle_radius", &ParticleRadius,
        "Particle radius from charge and atomic mass.", py::arg("charge"),
        py::arg("AtomicNumber"));

  m.def("dee_to_dpp", &dee_to_dpp, "DE/E to DP/P", py::arg("dee"),
        py::arg("beta"));

  m.def("dpp_to_dee", &dpp_to_dee, "DP/P to DE/E", py::arg("dpp"),
        py::arg("beta"));

  /*
 ================================================================================
                              LONGITUDINAL / RF
 ================================================================================
  */
  m.def("sige_from_sigs", &sigefromsigs, "Energy spread from bunch length.",
        py::arg("omega0"), py::arg("sigs"), py::arg("qs"), py::arg("gamma"),
        py::arg("gammatr"));

  m.def("sigs_from_sige", &sigsfromsige, "Bunch length from energy spread.",
        py::arg("sige"), py::arg("gamma"), py::arg("gammatr"),
        py::arg("omegas"));

  m.def("rf_voltage_in_ev",
        [](double phi, double c, std::vector<double> h, std::vector<double> v) {
          return EffectiveRFVoltageInElectronVolt(phi, c, h.size(), h.data(),
                                                  v.data());
        },
        "RF energy gain per turn [eV]", py::arg("phi"), py::arg("charge"),
        py::arg("harmonicNumbers"), py::arg("voltages"));

  m.def("rf_voltage_in_ev_prime",
        [](double phi, double c, std::vector<double> h, std::vector<double> v) {
          return EffectiveRFVoltageInElectronVoltPrime(phi, c, h.size(),
                                                       h.data(), v.data());
        },
        "Derivative of RF energy gain per turn [eV]", py::arg("phi"),
        py::arg("charge"), py::arg("harmonicNumbers"), py::arg("voltages"));

  m.def("rf_voltage_in_ev_with_rad_losses",
        [](double phi, double U0, double c, std::vector<double> h,
           std::vector<double> v) {
          return VeffRFeVRadlosses(phi, c, U0, h.size(), h.data(), v.data());
        },
        "RF energy gain per turn [eV] with Radiation Losses", py::arg("phi"),
        py::arg("U0"), py::arg("charge"), py::arg("harmonicNumbers"),
        py::arg("voltages"));

  m.def("get_synchronuous_phase",
        [](double target, double init_phi, double U0, double c,
           std::vector<double> h, std::vector<double> v, double e) {
          return SynchronuousPhase(target, init_phi, U0, c, h.size(), h.data(),
                                   v.data(), e);
        },
        "Synchronuous Phase", py::arg("targetEnergyGain"), py::arg("init_phi"),
        py::arg("U0"), py::arg("charge"), py::arg("harmonicNumbers"),
        py::arg("voltages"), py::arg("epsilon"));

  m.def("rf_voltage_with_potential_well_distortion",
        [](double target, double U0, double c, std::vector<double> h,
           std::vector<double> v, double L, double N, double sigs, double pc) {
          return VeffRFeVPotentialWellDistortion(
              target, U0, c, h.size(), h.data(), v.data(), L, N, sigs, pc);
        },
        "RF energy gain [eV] with Potential Well Distortion", py::arg("phi"),
        py::arg("U0"), py::arg("charge"), py::arg("harmonicNumbers"),
        py::arg("voltages"), py::arg("L"), py::arg("N"), py::arg("sigs"),
        py::arg("pc"));

  m.def("rf_voltage_with_potential_well_distortion_prime",
        [](double target, double U0, double c, std::vector<double> h,
           std::vector<double> v, double L, double N, double sigs, double pc) {
          return VeffRFeVPotentialWellDistortionPrime(
              target, U0, c, h.size(), h.data(), v.data(), L, N, sigs, pc);
        },
        "Derivative of RF energy gain [eV] with Potential Well Distortion",
        py::arg("phi"), py::arg("U0"), py::arg("charge"),
        py::arg("harmonicNumbers"), py::arg("voltages"), py::arg("L"),
        py::arg("N"), py::arg("sigs"), py::arg("pc"));

  m.def("get_synchronuous_phase_with_potential_well_distortion",
        [](double target, double init_phi, double U0, double c,
           std::vector<double> h, std::vector<double> v, double L, double N,
           double sigs, double pc, double e) {
          return SynchronuousPhaseWithPWD(target, init_phi, U0, c, h.size(),
                                          h.data(), v.data(), L, N, sigs, pc,
                                          e);
        },
        "Synchronuous Phase with Potential Well Distortion",
        py::arg("targetEnergyGain"), py::arg("init_phi"), py::arg("U0"),
        py::arg("charge"), py::arg("harmonicNumbers"), py::arg("voltages"),
        py::arg("L"), py::arg("N"), py::arg("sigs"), py::arg("pc"),
        py::arg("epsilon"));

  m.def("get_synchrotron_tune",
        [](double omega0, double U0, double c, std::vector<double> h,
           std::vector<double> v, double phis, double eta, double pc) {
          return SynchrotronTune(omega0, U0, c, h.size(), h.data(), v.data(),
                                 phis, eta, pc);
        },
        "Synchronuous Tune", py::arg("omega0"), py::arg("U0"),
        py::arg("charge"), py::arg("harmonicNumbers"), py::arg("voltages"),
        py::arg("phis"), py::arg("eta"), py::arg("pc"));

  m.def("get_synchrotron_tune_with_potential_well",
        [](double omega0, double U0, double c, std::vector<double> h,
           std::vector<double> v, double L, double N, double sigs, double phis,
           double eta, double pc) {
          return SynchrotronTunePWD(omega0, U0, c, h.size(), h.data(), v.data(),
                                    L, N, sigs, phis, eta, pc);
        },
        "Synchronuous Tune ith Potential Well Distortion", py::arg("omega0"),
        py::arg("U0"), py::arg("charge"), py::arg("harmonicNumbers"),
        py::arg("voltages"), py::arg("L"), py::arg("N"), py::arg("sigs"),
        py::arg("phis"), py::arg("eta"), py::arg("pc"));

  m.def("sige_form_sigs_using_rf",
        [](double sigs, double U0, double c, std::vector<double> h,
           std::vector<double> v, double gamma, double gammatr, double pc,
           double circ, double phis, bool printout) {
          return SigeFromRFAndSigs(sigs, U0, c, h.size(), h.data(), v.data(),
                                   gamma, gammatr, pc, circ, phis, printout);
        },
        "Sigma E from Sigma s using RF settings.", py::arg(""), py::arg("U0"),
        py::arg("charge"), py::arg("harmonicNumbers"), py::arg("voltages"),
        py::arg("gamma"), py::arg("gammaTransition"), py::arg("pc"),
        py::arg("circ"), py::arg("phis"), py::arg("printout"));

  /*
================================================================================
                       IBS NUMERIC FUNCTIONS
================================================================================
*/
  m.def("fmohl", &fmohl, "Fmohl function");

  m.def("rds", &rds, "Nagaitsev paper rds function");

  /*
================================================================================
                       COULOMB LOG
================================================================================
*/
  m.def("twclog",
        [](double pnumber, double bx, double by, double dx, double dy,
           double ex, double ey, double r0, double gamma, double charge,
           double en0, double amass, double sige, double sigt,
           py::array_t<double> npclog) {
          double clog[2];
          twclog(pnumber, bx, by, dx, dy, ex, ey, r0, gamma, charge, en0, amass,
                 sige, sigt, clog);
          auto buf_clog = npclog.request();
          double *ptr_clog = static_cast<double *>(buf_clog.ptr);
          ptr_clog[0] = clog[0];
          ptr_clog[1] = clog[1];
        },
        "Calculate the Coulomb Log.", py::arg("pnumber"), py::arg("betax"),
        py::arg("betay"), py::arg("dispersionx"), py::arg("disperisony"),
        py::arg("emitx"), py::arg("emity"), py::arg("classicalRadius"),
        py::arg("gamma"), py::arg("charge"), py::arg("energy"),
        py::arg("restEnergy_GeV"), py::arg("sige"), py::arg("sigt"),
        py::arg("outputArray"));

  m.def("twclogtail",
        [](double pnumber, double l, double bx, double by, double dx,
           double dpx, double dy, double dpy, double ax, double ay,
           double angle, double k1l, double k1sl, double ex, double ey,
           double r0, double aatom, double gamma, double en0, double len,
           double amass, double charge, double sige, double sigt,
           py::array_t<double> npclog) {
          double clog[2];
          twclogtail(pnumber, l, bx, by, dx, dpx, dy, dpy, ax, ay, angle, k1l,
                     k1sl, ex, ey, r0, aatom, gamma, en0, len, amass, charge,
                     sige, sigt, clog);
          auto buf_clog = npclog.request();
          double *ptr_clog = static_cast<double *>(buf_clog.ptr);
          ptr_clog[0] = clog[0];
          ptr_clog[1] = clog[1];
        },
        "Calculate the Coulomb Log with Tailcut.", py::arg("pnumber"),
        py::arg("elementLength"), py::arg("betax"), py::arg("betay"),
        py::arg("dispersionx"), py::arg("dispersionx_der"),
        py::arg("dispersiony"), py::arg("dispersiony_der"), py::arg("alphax"),
        py::arg("alphay"), py::arg("angle"), py::arg("k1l"), py::arg("k1sl"),
        py::arg("emitx"), py::arg("emity"), py::arg("classicalRadius"),
        py::arg("AtomicMassNumber"), py::arg("gamma"), py::arg("energy"),
        py::arg("acceleratorLength"), py::arg("restEnergy_GeV"),
        py::arg("charge"), py::arg("sige"), py::arg("sigt"),
        py::arg("outputArray"));

  m.def("coulomblog",
        [](double pnumber, double ex, double ey,
           map<string, double> &twissheader, double sige, double sigt,
           double r0, bool printout, py::array_t<double> npclog) {
          double clog[2];
          CoulombLog(pnumber, ex, ey, twissheader, sige, sigt, r0, printout,
                     clog);
          auto buf_clog = npclog.request();
          double *ptr_clog = static_cast<double *>(buf_clog.ptr);
          ptr_clog[0] = clog[0];
          ptr_clog[1] = clog[1];
        },
        "Calculate the Coulomb Log using ring average.", py::arg("pnumber"),
        py::arg("emitx"), py::arg("emity"), py::arg("twissheader"),
        py::arg("sige"), py::arg("sigt"), py::arg("classialRadius"),
        py::arg("printout"), py::arg("outputArray"));

  m.def("coulomblogtail",
        [](double pnumber, double ex, double ey,
           map<string, double> &twissheader, double sige, double sigt,
           double tauradx, double taurady, double taurads, double r0,
           bool printout, py::array_t<double> npclog) {
          double clog[2];
          TailCutCoulombLog(pnumber, ex, ey, twissheader, sige, sigt, tauradx,
                            taurady, taurads, r0, printout, clog);
          auto buf_clog = npclog.request();
          double *ptr_clog = static_cast<double *>(buf_clog.ptr);
          ptr_clog[0] = clog[0];
          ptr_clog[1] = clog[1];
        },
        "Calculate the Coulomb Log with tailcut using ring average.",
        py::arg("pnumber"), py::arg("emitx"), py::arg("emity"),
        py::arg("twissheader"), py::arg("sige"), py::arg("sigt"),
        py::arg("tauradx"), py::arg("taurady"), py::arg("taurads"),
        py::arg("classialRadius"), py::arg("printout"), py::arg("outputArray"));
  /*
================================================================================
                    INTEGRATORS
================================================================================
*/
  m.def("integrator_simpson_decade",
        [](double a, double b, double c, double cl, double cx, double cy,
           double cprime, double cyy, double tl1, double tl2, double tx1,
           double tx2, double ty1, double ty2, py::array_t<double> tau) {
          double alpha[3];
          SimpsonDecade(a, b, c, cl, cx, cy, cprime, cyy, tl1, tl2, tx1, tx2,
                        ty1, ty2, alpha);
          auto buf_tau = tau.request();
          double *ptr_tau = static_cast<double *>(buf_tau.ptr);
          ptr_tau[0] = alpha[0];
          ptr_tau[1] = alpha[1];
          ptr_tau[2] = alpha[2];
        },
        "Simpson integrator with decade splitting.", py::arg("a"), py::arg("b"),
        py::arg("c"), py::arg("cl"), py::arg("cx"), py::arg("cy"),
        py::arg("cprime"), py::arg("cyy"), py::arg("tl1"), py::arg("tl2"),
        py::arg("tx1"), py::arg("tx2"), py::arg("ty1"), py::arg("ty2"),
        py::arg("tau"));
  m.def("integrator_twsint",
        [](double pnumber, double ex, double ey, double sigs, double sige,
           double gammas, double betax, double betay, double alx, double aly,
           double dx, double dpx, double dy, double dpy,
           py::array_t<double> tau) {
          double alpha[3];
          twsint(pnumber, ex, ey, sigs, sige, gammas, betax, betay, alx, aly,
                 dx, dpx, dy, dpy, alpha);
          auto buf_tau = tau.request();
          double *ptr_tau = static_cast<double *>(buf_tau.ptr);
          ptr_tau[0] = alpha[0];
          ptr_tau[1] = alpha[1];
          ptr_tau[2] = alpha[2];
        },
        "MADX Simpson Decade with scaling integrator", py::arg("pnumber"),
        py::arg("emitx"), py::arg("emity"), py::arg("sigs"), py::arg("sige"),
        py::arg("gamma"), py::arg("betax"), py::arg("betay"), py::arg("alphax"),
        py::arg("alphay"), py::arg("dispersionx"), py::arg("dispersionx_der"),
        py::arg("dispersiony"), py::arg("dispersiony_der"), py::arg("tau"));
  m.def("integrand", &IBSIntegralIntegrand, "IBS integral integrand ",
        py::arg("int_var"), py::arg("ax"), py::arg("bx"), py::arg("a"),
        py::arg("b"), py::arg("c"));
  m.def("integrator_simpson", &simpson,
        "standard simpson integration on variable integrand (function)"); /*
[](const std::function<double(double, double, double, double, double,
        double)> &ibsintegrand,
double ax, double bx, double a, double b, double c, double al,
double bl,
int n) { return simpson(ibsintegrand, ax, bx, a, b, c, al, bl, n);
}, "Standard Simpson integrator");*/
  m.def("integrator_simpson_bjorken_mtingwa",
        [](const std::function<double(double, double, double, double, double,
                                      double)> &ibsintegrand,
           double ax, double bx, double ay, double by, double as, double bs,
           double a, double b, double ci, py::array_t<double> tau) {
          double alpha[3];
          intSimpson(ibsintegrand, ax, bx, ay, by, as, bs, a, b, ci, alpha);
          auto buf_tau = tau.request();
          double *ptr_tau = static_cast<double *>(buf_tau.ptr);
          ptr_tau[0] = alpha[0];
          ptr_tau[1] = alpha[1];
          ptr_tau[2] = alpha[2];
        },
        "Standard Simpson integral on Bjorken-Mtingwa integrand.");
  m.def("integral_bjorken_mtingwa",
        [](double pnumber, double ex, double ey, double sigs, double sige,
           double gammas, double betx, double bety, double alx, double aly,
           double dx, double dpx, double dy, double dpy,
           py::array_t<double> tau) {
          double alpha[3];
          BjorkenMtingwaInt(pnumber, ex, ey, sigs, sige, gammas, betx, bety,
                            alx, aly, dx, dpx, dy, dpy, alpha);
          auto buf_tau = tau.request();
          double *ptr_tau = static_cast<double *>(buf_tau.ptr);
          ptr_tau[0] = alpha[0];
          ptr_tau[1] = alpha[1];
          ptr_tau[2] = alpha[2];
        },
        "Bjorken-Mtingwa integral calculated using Simpson Decade",
        py::arg("pnumber"), py::arg("emitx"), py::arg("emity"), py::arg("sigs"),
        py::arg("sige"), py::arg("gamma"), py::arg("betax"), py::arg("betay"),
        py::arg("alphax"), py::arg("alphay"), py::arg("dispersionx"),
        py::arg("dispersionx_der"), py::arg("dispersiony"),
        py::arg("dispersiony_der"), py::arg("tau"));
  m.def("integral_conte_martini",
        [](double pnumber, double ex, double ey, double sigs, double sige,
           double gammas, double betx, double bety, double alx, double aly,
           double dx, double dpx, double dy, double dpy,
           py::array_t<double> tau) {
          double alpha[3];
          ConteMartiniInt(pnumber, ex, ey, sigs, sige, gammas, betx, bety, alx,
                          aly, dx, dpx, dy, dpy, alpha);
          auto buf_tau = tau.request();
          double *ptr_tau = static_cast<double *>(buf_tau.ptr);
          ptr_tau[0] = alpha[0];
          ptr_tau[1] = alpha[1];
          ptr_tau[2] = alpha[2];
        },
        "Conte-Martini integral calculated using Simpson Decade",
        py::arg("pnumber"), py::arg("emitx"), py::arg("emity"), py::arg("sigs"),
        py::arg("sige"), py::arg("gamma"), py::arg("betax"), py::arg("betay"),
        py::arg("alphax"), py::arg("alphay"), py::arg("dispersionx"),
        py::arg("dispersionx_der"), py::arg("dispersiony"),
        py::arg("dispersiony_der"), py::arg("tau"));
  m.def("integral_zimmerman",
        [](double pnumber, double ex, double ey, double sigs, double sige,
           double gammas, double betx, double bety, double alx, double aly,
           double dx, double dpx, double dy, double dpy,
           py::array_t<double> tau) {
          double alpha[3];
          MadxInt(pnumber, ex, ey, sigs, sige, gammas, betx, bety, alx, aly, dx,
                  dpx, dy, dpy, alpha);
          auto buf_tau = tau.request();
          double *ptr_tau = static_cast<double *>(buf_tau.ptr);
          ptr_tau[0] = alpha[0];
          ptr_tau[1] = alpha[1];
          ptr_tau[2] = alpha[2];
        },
        "Zimmerman integral calculated using Simpson Decade",
        py::arg("pnumber"), py::arg("emitx"), py::arg("emity"), py::arg("sigs"),
        py::arg("sige"), py::arg("gamma"), py::arg("betax"), py::arg("betay"),
        py::arg("alphax"), py::arg("alphay"), py::arg("dispersionx"),
        py::arg("dispersionx_der"), py::arg("dispersiony"),
        py::arg("dispersiony_der"), py::arg("tau"));
  /*
================================================================================
                 IBS MODELS
================================================================================
*/
  m.def("PiwinskiSmooth",
        [](double pnumber, double ex, double ey, double sigs, double dponp,
           map<string, double> &header, double r0, py::array_t<double> out) {
          double *ibs;
          ibs = PiwinskiSmooth(pnumber, ex, ey, sigs, dponp, header, r0);

          auto buf_out = out.request();
          double *ptr_out = static_cast<double *>(buf_out.ptr);
          ptr_out[0] = ibs[0];
          ptr_out[1] = ibs[1];
          ptr_out[2] = ibs[2];
        },
        "Piwinski smooth");

  m.def("PiwinskiLattice",
        [](double pnumber, double ex, double ey, double sigs, double dponp,
           map<string, double> &header, map<string, vector<double>> &table,
           double r0, py::array_t<double> out) {
          double *ibs;
          ibs =
              PiwinskiLattice(pnumber, ex, ey, sigs, dponp, header, table, r0);

          auto buf_out = out.request();
          double *ptr_out = static_cast<double *>(buf_out.ptr);
          ptr_out[0] = ibs[0];
          ptr_out[1] = ibs[1];
          ptr_out[2] = ibs[2];
        },
        "Piwinski Lattice");

  m.def("PiwinskiLatticeModified",
        [](double pnumber, double ex, double ey, double sigs, double dponp,
           map<string, double> &header, map<string, vector<double>> &table,
           double r0, py::array_t<double> out) {
          double *ibs;
          ibs = PiwinskiLatticeModified(pnumber, ex, ey, sigs, dponp, header,
                                        table, r0);

          auto buf_out = out.request();
          double *ptr_out = static_cast<double *>(buf_out.ptr);
          ptr_out[0] = ibs[0];
          ptr_out[1] = ibs[1];
          ptr_out[2] = ibs[2];
        },
        "Piwinski Lattice Modified");

  m.def("Nagaitsev",
        [](double pnumber, double ex, double ey, double sigs, double dponp,
           map<string, double> &header, map<string, vector<double>> &table,
           double r0, py::array_t<double> out) {
          double *ibs;
          ibs = Nagaitsev(pnumber, ex, ey, sigs, dponp, header, table, r0);

          auto buf_out = out.request();
          double *ptr_out = static_cast<double *>(buf_out.ptr);
          ptr_out[0] = ibs[0];
          ptr_out[1] = ibs[1];
          ptr_out[2] = ibs[2];
        },
        "Nagaitsev");

  m.def("NagaitsevTailcut",
        [](double pnumber, double ex, double ey, double sigs, double dponp,
           map<string, double> &header, map<string, vector<double>> &table,
           double aatom, double r0, py::array_t<double> out) {
          double *ibs;
          ibs = Nagaitsevtailcut(pnumber, ex, ey, sigs, dponp, header, table,
                                 r0, aatom);

          auto buf_out = out.request();
          double *ptr_out = static_cast<double *>(buf_out.ptr);
          ptr_out[0] = ibs[0];
          ptr_out[1] = ibs[1];
          ptr_out[2] = ibs[2];
        },
        "Nagaitsev Tailcut");

  m.def("Zimmerman",
        [](double pnumber, double ex, double ey, double sigs, double sige,
           map<string, double> &header, map<string, vector<double>> &table,
           double r0, bool printout, py::array_t<double> out) {
          double *ibs;
          ibs =
              ibsmadx(pnumber, ex, ey, sigs, sige, header, table, r0, printout);

          auto buf_out = out.request();
          double *ptr_out = static_cast<double *>(buf_out.ptr);
          ptr_out[0] = ibs[0];
          ptr_out[1] = ibs[1];
          ptr_out[2] = ibs[2];
        },
        "Zimmerman");

  m.def("ZimmermanTailcut",
        [](double pnumber, double ex, double ey, double sigs, double sige,
           map<string, double> &header, map<string, vector<double>> &table,
           double r0, bool aatom, py::array_t<double> out) {
          double *ibs;
          ibs = ibsmadxtailcut(pnumber, ex, ey, sigs, sige, header, table, r0,
                               aatom);

          auto buf_out = out.request();
          double *ptr_out = static_cast<double *>(buf_out.ptr);
          ptr_out[0] = ibs[0];
          ptr_out[1] = ibs[1];
          ptr_out[2] = ibs[2];
        },
        "Zimmerman Tailcut");

  m.def("BjorkenMtingwaSimpson",
        [](double pnumber, double ex, double ey, double sigs, double dponp,
           map<string, double> &header, map<string, vector<double>> &table,
           double r0, bool aatom, py::array_t<double> out) {
          double *ibs;
          ibs =
              BjorkenMtingwa2(pnumber, ex, ey, sigs, dponp, header, table, r0);

          auto buf_out = out.request();
          double *ptr_out = static_cast<double *>(buf_out.ptr);
          ptr_out[0] = ibs[0];
          ptr_out[1] = ibs[1];
          ptr_out[2] = ibs[2];
        },
        "Bjorken-Mtingwa using Standard Simpson integration.");

  m.def("BjorkenMtingwaSimpsonDecade",
        [](double pnumber, double ex, double ey, double sigs, double sige,
           map<string, double> &header, map<string, vector<double>> &table,
           double r0, bool aatom, py::array_t<double> out) {
          double *ibs;
          ibs = BjorkenMtingwa(pnumber, ex, ey, sigs, sige, header, table, r0);

          auto buf_out = out.request();
          double *ptr_out = static_cast<double *>(buf_out.ptr);
          ptr_out[0] = ibs[0];
          ptr_out[1] = ibs[1];
          ptr_out[2] = ibs[2];
        },
        "Bjorken-Mtingwa using Simpson Decade integration.");

  m.def("BjorkenMtingwaTailcutSimpsonDecade",
        [](double pnumber, double ex, double ey, double sigs, double sige,
           map<string, double> &header, map<string, vector<double>> &table,
           double r0, bool aatom, py::array_t<double> out) {
          double *ibs;
          ibs = BjorkenMtingwatailcut(pnumber, ex, ey, sigs, sige, header,
                                      table, r0, aatom);

          auto buf_out = out.request();
          double *ptr_out = static_cast<double *>(buf_out.ptr);
          ptr_out[0] = ibs[0];
          ptr_out[1] = ibs[1];
          ptr_out[2] = ibs[2];
        },
        "Bjorken-Mtingwa using Simpson Decade integration with Tailcut");

  m.def("ConteMartiniSimpsonDecade",
        [](double pnumber, double ex, double ey, double sigs, double sige,
           map<string, double> &header, map<string, vector<double>> &table,
           double r0, bool aatom, py::array_t<double> out) {
          double *ibs;
          ibs = ConteMartini(pnumber, ex, ey, sigs, sige, header, table, r0);

          auto buf_out = out.request();
          double *ptr_out = static_cast<double *>(buf_out.ptr);
          ptr_out[0] = ibs[0];
          ptr_out[1] = ibs[1];
          ptr_out[2] = ibs[2];
        },
        "Conte-Martini using Simpson Decade integration.");

  m.def("ConteMartiniTailcutSimpsonDecade",
        [](double pnumber, double ex, double ey, double sigs, double sige,
           map<string, double> &header, map<string, vector<double>> &table,
           double r0, bool aatom, py::array_t<double> out) {
          double *ibs;
          ibs = ConteMartinitailcut(pnumber, ex, ey, sigs, sige, header, table,
                                    r0, aatom);

          auto buf_out = out.request();
          double *ptr_out = static_cast<double *>(buf_out.ptr);
          ptr_out[0] = ibs[0];
          ptr_out[1] = ibs[1];
          ptr_out[2] = ibs[2];
        },
        "Conte-Martini using Simpson Decade integration with Tailcut");

  m.def("ZimmermanSimpsonDecade",
        [](double pnumber, double ex, double ey, double sigs, double sige,
           map<string, double> &header, map<string, vector<double>> &table,
           double r0, bool printout, py::array_t<double> out) {
          double *ibs;
          ibs =
              ibsmadx(pnumber, ex, ey, sigs, sige, header, table, r0, printout);

          auto buf_out = out.request();
          double *ptr_out = static_cast<double *>(buf_out.ptr);
          ptr_out[0] = ibs[0];
          ptr_out[1] = ibs[1];
          ptr_out[2] = ibs[2];
        },
        "Zimmerman using Simpson Decade");

  m.def("runODE",
        [](map<string, double> &twiss, map<string, vector<double>> &twissdata,
           vector<double> h, vector<double> v, vector<double> &t,
           vector<double> &ex, vector<double> &ey, vector<double> &sigs,
           vector<double> sige, int model, double pnumber,
           int couplingpercentage, double threshold, string method) {
          ODE(twiss, twissdata, h.size(), h.data(), v.data(), t, ex, ey, sigs,
              sige, model, pnumber, couplingpercentage, threshold, method);
          map<string, vector<double>> res;
          res["t"] = t;
          res["ex"] = ex;
          res["ey"] = ey;
          res["sigs"] = sigs;
          return res;
        },
        "");
  m.def("runODE",
        [](map<string, double> &twiss, map<string, vector<double>> &twissdata,
           vector<double> h, vector<double> v, vector<double> &t,
           vector<double> &ex, vector<double> &ey, vector<double> &sigs,
           vector<double> sige, int model, double pnumber, int nsteps,
           double stepsize, int couplingpercentage, string method) {
          ODE(twiss, twissdata, h.size(), h.data(), v.data(), t, ex, ey, sigs,
              sige, model, pnumber, nsteps, stepsize, couplingpercentage,
              method);
          map<string, vector<double>> res;
          res["t"] = t;
          res["ex"] = ex;
          res["ey"] = ey;
          res["sigs"] = sigs;
          return res;
        },
        "");
}
