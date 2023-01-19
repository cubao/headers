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
             "x"_a, "y"_a, "z"_a)
        .def("ecef2lla",
             py::overload_cast<const Eigen::Ref<const RowVectors> &>(ecef2lla),
             "ecefs"_a)
        .def("lla2ecef", py::overload_cast<double, double, double>(lla2ecef),
             "lon"_a, "lat"_a, "alt"_a)
        .def("lla2ecef",
             py::overload_cast<const Eigen::Ref<const RowVectors> &>(lla2ecef),
             "llas"_a)
        // lla <-> enu
        .def("lla2enu", &lla2enu, "llas"_a, py::kw_only(), //
             CUBAO_ARGV_DEFAULT_NONE(anchor_lla), "cheap_ruler"_a = true)
        .def("enu2lla", &enu2lla, "enus"_a, py::kw_only(), //
             "anchor_lla"_a, "cheap_ruler"_a = true)
        // enu <-> ecef
        .def("enu2ecef", &enu2ecef, "enus"_a, py::kw_only(), //
             "anchor_lla"_a, "cheap_ruler"_a = false)
        .def("ecef2enu", &ecef2enu, "ecefs"_a, py::kw_only(), //
             CUBAO_ARGV_DEFAULT_NONE(anchor_lla), "cheap_ruler"_a = false)
        // T_ecef_enu
        .def("R_ecef_enu", &R_ecef_enu, "lon"_a, "lat"_a)
        .def("T_ecef_enu",
             py::overload_cast<double, double, double>(&T_ecef_enu), //
             "lon"_a, "lat"_a, "alt"_a)
        .def("T_ecef_enu",
             py::overload_cast<const Eigen::Vector3d &>(&T_ecef_enu), "lla"_a)
        // apply transform
        .def("apply_transform", &apply_transform, "T"_a, "coords"_a)
        .def("apply_transform_inplace", &apply_transform_inplace, //
             "T"_a, "coords"_a, py::kw_only(), "batch_size"_a = 1000)
        //
        .def("cheap_ruler_k", &cheap_ruler_k, "latitude"_a)
        //
        ;
}
} // namespace cubao
