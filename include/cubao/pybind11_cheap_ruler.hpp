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
    py::class_<CheapRuler> cheap_ruler(m, "CheapRuler", py::module_local(),
                                       R"docstring(
        A class for fast distance calculations and geometric operations.

        CheapRuler provides methods for various geometric calculations
        optimized for speed and simplicity, sacrificing some accuracy
        for performance.
        )docstring");

    py::enum_<CheapRuler::Unit>(cheap_ruler, "Unit",
                                R"docstring(
        Enumeration of supported distance units.
        )docstring")
        .value("Kilometers", CheapRuler::Unit::Kilometers, "Kilometers")
        .value("Miles", CheapRuler::Unit::Miles, "Miles")
        .value("NauticalMiles", CheapRuler::Unit::NauticalMiles,
               "Nautical Miles")
        .value("Meters", CheapRuler::Unit::Meters, "Meters")
        .value("Metres", CheapRuler::Unit::Metres, "Metres (alias for Meters)")
        .value("Yards", CheapRuler::Unit::Yards, "Yards")
        .value("Feet", CheapRuler::Unit::Feet, "Feet")
        .value("Inches", CheapRuler::Unit::Inches, "Inches")
        .export_values();

    cheap_ruler //
        .def("k", py::overload_cast<>(&CheapRuler::k, py::const_),
             "Get the ruler's unit conversion factor.")
        .def_static(
            "_k", py::overload_cast<double, CheapRuler::Unit>(&CheapRuler::k),
            "latitude"_a, py::kw_only(), "unit"_a = CheapRuler::Unit::Meters,
            "Get the unit conversion factor for a given latitude and unit.")

        .def_readonly_static("RE", &CheapRuler::RE,
                             "Earth's equatorial radius in meters.")
        .def_readonly_static("FE", &CheapRuler::FE, "Earth's flattening.")
        .def_readonly_static("E2", &CheapRuler::E2,
                             "Square of Earth's eccentricity.")
        .def_readonly_static("RAD", &CheapRuler::RAD,
                             "Conversion factor from degrees to radians.")
        //
        .def(py::init<double, CheapRuler::Unit>(), "latitude"_a, py::kw_only(),
             "unit"_a = CheapRuler::Unit::Meters,
             "Initialize a CheapRuler object with a given latitude and unit.")
        //
        .def_static("_fromTile", &CheapRuler::fromTile, "x"_a, "y"_a,
                    "Create a CheapRuler from tile coordinates (x, y).")
        .def("delta", &CheapRuler::delta, "lla0"_a, "lla1"_a,
             "Calculate the distance between two points in the x, y plane.")
        .def("squareDistance", &CheapRuler::squareDistance, "a"_a, "b"_a,
             "Calculate the squared distance between two points.")
        .def("distance", &CheapRuler::distance, "a"_a, "b"_a,
             "Calculate the distance between two points.")
        .def("bearing", &CheapRuler::bearing, "a"_a, "b"_a,
             "Calculate the bearing between two points.")
        .def("destination", &CheapRuler::destination, "origin"_a, "dist"_a,
             "bearing"_a,
             "Calculate the destination point given origin, distance, and "
             "bearing.")
        .def("offset", &CheapRuler::offset, "origin"_a, "dx"_a, "dy"_a,
             "dz"_a = 0.0, "Calculate a new point given origin and offsets.")
        .def("lineDistance", &CheapRuler::lineDistance, "points"_a,
             "Calculate the total distance of a line (an array of points).")
        .def("area", &CheapRuler::area, "ring"_a,
             "Calculate the area of a polygon.")
        .def("along", &CheapRuler::along, "line"_a, "dist"_a,
             "Calculate a point at a specified distance along the line.")
        .def("pointToSegmentDistance", &CheapRuler::pointToSegmentDistance,
             "p"_a, "a"_a, "b"_a,
             "Calculate the distance from a point to a line segment.")
        .def("pointOnLine", &CheapRuler::pointOnLine, "line"_a, "p"_a,
             "Calculate the closest point on a line to the given point.")
        .def("lineSlice", &CheapRuler::lineSlice, "start"_a, "stop"_a, "line"_a,
             "Get a part of the line between the start and stop points.")
        .def("lineSliceAlong", &CheapRuler::lineSliceAlong, "start"_a, "stop"_a,
             "line"_a,
             "Get a part of the line between the start and stop distances "
             "along the line.")
        .def("bufferPoint", &CheapRuler::bufferPoint, "p"_a, "buffer"_a,
             "Create a bounding box around a point.")
        .def("bufferBBox", &CheapRuler::bufferBBox, "bbox"_a, "buffer"_a,
             "Create a bounding box around another bounding box.")
        .def_static("_insideBBox", &CheapRuler::insideBBox, "p"_a, "bbox"_a,
                    py::kw_only(), "cheak_z"_a = false,
                    "Check if a point is inside a bounding box.")
        .def_static("_interpolate", &CheapRuler::interpolate, "a"_a, "b"_a,
                    "t"_a, "Interpolate linearly between two points.")
        .def_static("_longDiff", &CheapRuler::longDiff, "a"_a, "b"_a,
                    "Calculate the difference between two longitudes.");
}
} // namespace cubao
