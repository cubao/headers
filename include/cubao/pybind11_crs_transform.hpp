// should sync
// -
// https://github.com/cubao/polyline-ruler/blob/master/src/pybind11_crs_transform.hpp
// -
// https://github.com/cubao/headers/tree/main/include/cubao/pybind11_crs_transform.hpp

#pragma once

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "cubao_inline.hpp"
#include "crs_transform.hpp"

namespace cubao
{
namespace py = pybind11;
using namespace pybind11::literals;
using rvp = py::return_value_policy;

CUBAO_INLINE void bind_crs_transform(py::module &m)
{
    //
    m
        // ecef <-> lla
        .def("ecef2lla", py::overload_cast<double, double, double>(ecef2lla),
             "x"_a, "y"_a, "z"_a,
             "Convert ECEF coordinates to LLA (Longitude, Latitude, Altitude).")
        .def("ecef2lla",
             py::overload_cast<const Eigen::Ref<const RowVectors> &>(ecef2lla),
             "ecefs"_a,
             "Convert multiple ECEF coordinates to LLA (Longitude, Latitude, "
             "Altitude).")
        .def("lla2ecef", py::overload_cast<double, double, double>(lla2ecef),
             "lon"_a, "lat"_a, "alt"_a,
             "Convert LLA (Longitude, Latitude, Altitude) to ECEF coordinates.")
        .def("lla2ecef",
             py::overload_cast<const Eigen::Ref<const RowVectors> &>(lla2ecef),
             "llas"_a,
             "Convert multiple LLA (Longitude, Latitude, Altitude) to ECEF "
             "coordinates.")
        // lla <-> enu
        .def("lla2enu", &lla2enu, "llas"_a, py::kw_only(),
             CUBAO_ARGV_DEFAULT_NONE(anchor_lla), "cheap_ruler"_a = true,
             "Convert LLA (Longitude, Latitude, Altitude) to ENU (East, North, "
             "Up) coordinates.")
        .def("enu2lla", &enu2lla, "enus"_a, py::kw_only(), "anchor_lla"_a,
             "cheap_ruler"_a = true,
             "Convert ENU (East, North, Up) to LLA (Longitude, Latitude, "
             "Altitude) coordinates.")
        // enu <-> ecef
        .def("enu2ecef", &enu2ecef, "enus"_a, py::kw_only(), "anchor_lla"_a,
             "cheap_ruler"_a = false,
             "Convert ENU (East, North, Up) to ECEF coordinates.")
        .def("ecef2enu", &ecef2enu, "ecefs"_a, py::kw_only(),
             CUBAO_ARGV_DEFAULT_NONE(anchor_lla), "cheap_ruler"_a = false,
             "Convert ECEF to ENU (East, North, Up) coordinates.")
        // T_ecef_enu
        .def("R_ecef_enu", &R_ecef_enu, "lon"_a, "lat"_a,
             "Get rotation matrix from ECEF to ENU coordinate system.")
        .def("T_ecef_enu",
             py::overload_cast<double, double, double>(&T_ecef_enu), "lon"_a,
             "lat"_a, "alt"_a,
             "Get transformation matrix from ECEF to ENU coordinate system.")
        .def("T_ecef_enu",
             py::overload_cast<const Eigen::Vector3d &>(&T_ecef_enu), "lla"_a,
             "Get transformation matrix from ECEF to ENU coordinate system "
             "using LLA vector.")
        // apply transform
        .def("apply_transform", &apply_transform, "T"_a, "coords"_a,
             "Apply transformation matrix to coordinates.")
        .def("apply_transform_inplace", &apply_transform_inplace, "T"_a,
             "coords"_a, py::kw_only(), "batch_size"_a = 1000,
             "Apply transformation matrix to coordinates in-place.")
        //
        .def("cheap_ruler_k", &cheap_ruler_k, "latitude"_a,
             "Get the cheap ruler's unit conversion factor for a given "
             "latitude.")
        //
        ;
}
} // namespace cubao
