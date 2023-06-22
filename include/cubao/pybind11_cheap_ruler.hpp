// should sync
// -
// https://github.com/cubao/polyline-ruler/blob/master/src/pybind11_cheap_ruler.hpp
// -
// https://github.com/cubao/headers/tree/main/include/cubao/pybind11_cheap_ruler.hpp

#pragma once

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "cubao_inline.hpp"
#include "cheap_ruler.hpp"

namespace cubao
{
namespace py = pybind11;
using namespace pybind11::literals;
using rvp = py::return_value_policy;

CUBAO_INLINE void bind_cheap_ruler(py::module &m)
{
    py::class_<CheapRuler> cheap_ruler(m, "CheapRuler", py::module_local());

    py::enum_<CheapRuler::Unit>(cheap_ruler, "Unit")
        .value("Kilometers", CheapRuler::Unit::Kilometers)
        .value("Miles", CheapRuler::Unit::Miles)
        .value("NauticalMiles", CheapRuler::Unit::NauticalMiles)
        .value("Meters", CheapRuler::Unit::Meters)
        .value("Metres", CheapRuler::Unit::Metres)
        .value("Yards", CheapRuler::Unit::Yards)
        .value("Feet", CheapRuler::Unit::Feet)
        .value("Inches", CheapRuler::Unit::Inches)
        .export_values();

    cheap_ruler //
        .def("k", py::overload_cast<>(&CheapRuler::k, py::const_))
        .def_static(
            "_k", py::overload_cast<double, CheapRuler::Unit>(&CheapRuler::k),
            "latitude"_a, py::kw_only(), "unit"_a = CheapRuler::Unit::Meters)

        .def_readonly_static("RE", &CheapRuler::RE)
        .def_readonly_static("FE", &CheapRuler::FE)
        .def_readonly_static("E2", &CheapRuler::E2)
        .def_readonly_static("RAD", &CheapRuler::RAD)
        //
        .def(py::init<double, CheapRuler::Unit>(), "latitude"_a, py::kw_only(),
             "unit"_a = CheapRuler::Unit::Meters)
        //
        .def_static("_fromTile", &CheapRuler::fromTile, "x"_a, "y"_a)
        .def("delta", &CheapRuler::delta, "lla0"_a, "lla1"_a)
        .def("squareDistance", &CheapRuler::squareDistance, "a"_a, "b"_a)
        .def("distance", &CheapRuler::distance, "a"_a, "b"_a)
        .def("bearing", &CheapRuler::bearing, "a"_a, "b"_a)
        .def("destination", &CheapRuler::destination, "origin"_a, "dist"_a,
             "bearing"_a)
        .def("offset", &CheapRuler::offset, "origin"_a, "dx"_a, "dy"_a,
             "dz"_a = 0.0)
        .def("lineDistance", &CheapRuler::lineDistance, "points"_a)
        .def("area", &CheapRuler::area, "ring"_a)
        .def("along", &CheapRuler::along, "line"_a, "dist"_a)
        .def("pointToSegmentDistance", &CheapRuler::pointToSegmentDistance,
             "p"_a, "a"_a, "b"_a)
        .def("pointOnLine", &CheapRuler::pointOnLine, "line"_a, "p"_a)
        .def("lineSlice", &CheapRuler::lineSlice, "start"_a, "stop"_a, "line"_a)
        .def("lineSliceAlong", &CheapRuler::lineSliceAlong, "start"_a, "stop"_a,
             "line"_a)
        .def("bufferPoint", &CheapRuler::bufferPoint, "p"_a, "buffer"_a)
        .def("bufferBBox", &CheapRuler::bufferBBox, "bbox"_a, "buffer"_a)
        .def_static("_insideBBox", &CheapRuler::insideBBox, "p"_a, "bbox"_a,
                    py::kw_only(), "cheak_z"_a = false)
        .def_static("_interpolate", &CheapRuler::interpolate, "a"_a, "b"_a,
                    "t"_a)
        .def_static("_longDiff", &CheapRuler::longDiff, "a"_a, "b"_a);
}
} // namespace cubao
