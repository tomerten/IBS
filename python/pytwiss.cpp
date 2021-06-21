#include <../cpp/include/ibs_bits/twiss.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(IBSLib, m) {
  m.doc() = "IBS library";
  m.def("GetTwissHeader", &GetTwissHeader,
        "Get the twiss header as a dictionary.");
  m.def("GetTwissTable", &GetTwissTableAsMap, "Get the twiss data table.");
}