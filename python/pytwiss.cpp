//#include "../cpp/include/ibs_bits/CoulombLogFunctions.hpp"
//#include "../cpp/include/ibs_bits/NumericFunctions.hpp"
//#include "../cpp/include/ibs_bits/RadiationDamping.hpp"
//#include "../cpp/include/ibs_bits/twiss.hpp"
#include <ibs>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(IBSLib, m) {
  m.doc() = "IBS library";
  m.def("GetTwissHeader", &GetTwissHeader,
        "Get the twiss header as a dictionary.");
  m.def("GetTwissTable", &GetTwissTableAsMap, "Get the twiss data table.");
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
  m.def("sige_from_sigs", &sigefromsigs, "Energy spread from bunch length.");
  m.def("sigs_from_sige", &sigsfromsige, "Bunch length from energy spread.");
  m.def("eta", &eta, "Slip factor");
  m.def("fmohl", &fmohl, "Fmohl function");
  m.def("particle_radius", &ParticleRadius,
        "Particle radius from charge and atomic mass.");
  m.def("beta_relativistic_from_gamma", &BetaRelativisticFromGamma,
        "Beta relativistic.");
  m.def("rds", &rds, "Nagaitsev paper rds function");
  m.def("rf_voltage_in_ev",
        [](double phi, double c, std::vector<double> h, std::vector<double> v) {
          return EffectiveRFVoltageInElectronVolt(phi, c, h.size(), h.data(),
                                                  v.data());
        });
  m.def("rf_voltage_in_ev_prime",
        [](double phi, double c, std::vector<double> h, std::vector<double> v) {
          return EffectiveRFVoltageInElectronVoltPrime(phi, c, h.size(),
                                                       h.data(), v.data());
        });
  m.def("rf_voltage_in_ev_with_rad_losses",
        [](double phi, double U0, double c, std::vector<double> h,
           std::vector<double> v) {
          return VeffRFeVRadlosses(phi, c, U0, h.size(), h.data(), v.data());
        });
  m.def("get_synchronuous_phase",
        [](double target, double init_phi, double U0, double c,
           std::vector<double> h, std::vector<double> v, double e) {
          return SynchronuousPhase(target, init_phi, c, U0, h.size(), h.data(),
                                   v.data(), e);
        });
  m.def("rf_voltage_with_potential_well_distortion",
        [](double target, double U0, double c, std::vector<double> h,
           std::vector<double> v, double L, double N, double sigs, double pc) {
          return VeffRFeVPotentialWellDistortion(
              target, U0, c, h.size(), h.data(), v.data(), L, N, sigs, pc);
        });
  m.def("rf_voltage_with_potential_well_distortion_prime",
        [](double target, double U0, double c, std::vector<double> h,
           std::vector<double> v, double L, double N, double sigs, double pc) {
          return VeffRFeVPotentialWellDistortionPrime(
              target, U0, c, h.size(), h.data(), v.data(), L, N, sigs, pc);
        });
  m.def("get_synchronuous_phase_with_potential_well_distortion",
        [](double target, double init_phi, double U0, double c,
           std::vector<double> h, std::vector<double> v, double L, double N,
           double sigs, double pc, double e) {
          return SynchronuousPhaseWithPWD(target, init_phi, U0, c, h.size(),
                                          h.data(), v.data(), L, N, sigs, pc,
                                          e);
        });
  m.def("get_synchrotron_tune",
        [](double omega0, double U0, double c, std::vector<double> h,
           std::vector<double> v, double phis, double eta, double pc) {
          return SynchrotronTune(omega0, U0, c, h.size(), h.data(), v.data(),
                                 phis, eta, pc);
        });
  m.def("get_synchrotron_tune_with_potential_well",
        [](double omega0, double U0, double c, std::vector<double> h,
           std::vector<double> v, double L, double N, double sigs, double phis,
           double eta, double pc) {
          return SynchrotronTunePWD(omega0, U0, c, h.size(), h.data(), v.data(),
                                    L, N, sigs, phis, eta, pc);
        });
  m.def("sige_form_sigs_using_rf",
        [](double sigs, double U0, double c, std::vector<double> h,
           std::vector<double> v, double gamma, double gammatr, double pc,
           double circ, double phis, bool printout) {
          return SigeFromRFAndSigs(sigs, U0, c, h.size(), h.data(), v.data(),
                                   gamma, gammatr, pc, circ, phis, printout);
        });
  m.def("updateTwiss",
        [](map<string, vector<double>> &table) {
          updateTwiss(table);
          return table;
        },
        "Extend Twiss Table with rad int, CS gamma, curly H and rho.");
  m.def("printTwissColumn", &printTwissMap, "Print Twiss column");
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
        "Calculate the Coulomb Log.");
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
        "Calculate the Coulomb Log with Tailcut.");
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
        "Calculate the Coulomb Log using ring average.");
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
        "Calculate the Coulomb Log with tailcut using ring average.");
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
        "Simpson integrator with decade splitting.");
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
        "MADX Simpson Decade with scaling integrator");
  m.def("integrand", &IBSIntegralIntegrand, "IBS integral integrand ");
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
        "Bjorken-Mtingwa integral calculated using Simpson Decade");
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
        "Conte-Martini integral calculated using Simpson Decade");
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
        "Zimmerman integral calculated using Simpson Decade");
}
