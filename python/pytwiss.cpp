#include "../cpp/include/ibs_bits/CoulombLogFunctions.hpp"
#include "../cpp/include/ibs_bits/NumericFunctions.hpp"
#include "../cpp/include/ibs_bits/RadiationDamping.hpp"
#include "../cpp/include/ibs_bits/twiss.hpp"
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
}