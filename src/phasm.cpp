#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "overlapper.h"

namespace py = pybind11;

PYBIND11_PLUGIN(overlapper) {
    py::module m("overlapper",
            "Optimised C++ all-vs-all exact overlap finder for PHASM");

    py::class_<ExactOverlapper>(m, "ExactOverlapper")
        .def(py::init())
        .def("add_sequence", &ExactOverlapper::addSequence)
        .def("overlaps", &ExactOverlapper::overlaps);

    return m.ptr();
}
