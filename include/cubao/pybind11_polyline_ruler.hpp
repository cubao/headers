// should sync
// -
// https://github.com/cubao/polyline-ruler/blob/master/src/pybind11_polyline_ruler.hpp
// -
// https://github.com/cubao/headers/tree/main/include/cubao/pybind11_polyline_ruler.hpp

#pragma once

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "cubao_inline.hpp"
#include "polyline_ruler.hpp"

namespace cubao
{
namespace py = pybind11;
using namespace pybind11::literals;
using rvp = py::return_value_policy;

CUBAO_INLINE void bind_polyline_ruler(py::module &m)
{
    m
        //
        .def(
            "intersect_segments",
            py::overload_cast<const Eigen::Vector2d &, const Eigen::Vector2d &,
                              const Eigen::Vector2d &, const Eigen::Vector2d &>(
                &intersect_segments), //
            "a1"_a, "a2"_a, "b1"_a, "b2"_a)
        .def(
            "intersect_segments",
            py::overload_cast<const Eigen::Vector3d &, const Eigen::Vector3d &,
                              const Eigen::Vector3d &, const Eigen::Vector3d &>(
                &intersect_segments), //
            "a1"_a, "a2"_a, "b1"_a, "b2"_a)
        //
        ;

    py::class_<LineSegment>(m, "LineSegment", py::module_local())      //
        .def(py::init<const Eigen::Vector3d, const Eigen::Vector3d>(), //
             "A"_a, "B"_a)
        .def("distance", &LineSegment::distance, "P"_a)
        .def("distance2", &LineSegment::distance2, "P"_a)
        .def("intersects", &LineSegment::intersects, "other"_a)
        .def_property_readonly(
            "length",
            [](const LineSegment &self) { return std::sqrt(self.len2); })
        .def_property_readonly(
            "length2", [](const LineSegment &self) { return self.len2; })
        .def_property_readonly(
            "A",
            [](const LineSegment &self) -> const Eigen::Vector3d & {
                return self.A;
            },
            rvp::reference_internal)
        .def_property_readonly(
            "B",
            [](const LineSegment &self) -> const Eigen::Vector3d & {
                return self.B;
            },
            rvp::reference_internal)
        .def_property_readonly(
            "AB",
            [](const LineSegment &self) -> const Eigen::Vector3d & {
                return self.AB;
            },
            rvp::reference_internal)
        //
        ;

    py::class_<PolylineRuler>(m, "PolylineRuler", py::module_local()) //
        .def(py::init<const Eigen::Ref<const RowVectors> &, bool>(),  //
             "coords"_a, py::kw_only(), "is_wgs84"_a = false)
        //
        .def("polyline", &PolylineRuler::polyline, rvp::reference_internal)
        .def("N", &PolylineRuler::N)
        .def("is_wgs84", &PolylineRuler::is_wgs84)
        .def("k", &PolylineRuler::k)
        //
        .def_static(
            "_ranges",
            py::overload_cast<const Eigen::Ref<const RowVectors> &, bool>(
                &PolylineRuler::ranges),
            "polyline"_a, py::kw_only(), "is_wgs84"_a = false)
        .def("ranges", py::overload_cast<>(&PolylineRuler::ranges, py::const_),
             rvp::reference_internal)
        .def("range", py::overload_cast<int>(&PolylineRuler::range, py::const_),
             "segment_index"_a)
        .def("range",
             py::overload_cast<int, double>(&PolylineRuler::range, py::const_),
             py::kw_only(), "segment_index"_a, "t"_a)
        //
        .def("segment_index", &PolylineRuler::segment_index, "range"_a)
        .def("segment_index_t", &PolylineRuler::segment_index_t, "range"_a)
        //
        .def("length", &PolylineRuler::length)
        //
        .def_static(
            "_dirs",
            py::overload_cast<const Eigen::Ref<const RowVectors> &, bool>(
                &PolylineRuler::dirs),
            "polyline"_a, py::kw_only(), "is_wgs84"_a = false)
        .def("dirs", py::overload_cast<>(&PolylineRuler::dirs, py::const_),
             rvp::reference_internal)
        //
        .def("dir", py::overload_cast<int>(&PolylineRuler::dir, py::const_),
             py::kw_only(), "point_index"_a)
        .def("dir",
             py::overload_cast<double, bool>(&PolylineRuler::dir, py::const_),
             py::kw_only(), "range"_a, "smooth_joint"_a = true)
        .def("extended_along", &PolylineRuler::extended_along, "range"_a)
        .def("at", py::overload_cast<double>(&PolylineRuler::at, py::const_),
             py::kw_only(), "range"_a)
        .def("at", py::overload_cast<int>(&PolylineRuler::at, py::const_),
             py::kw_only(), "segment_index"_a)
        .def("at",
             py::overload_cast<int, double>(&PolylineRuler::at, py::const_),
             py::kw_only(), "segment_index"_a, "t"_a)
        .def("arrow",
             py::overload_cast<int, double>(&PolylineRuler::arrow, py::const_),
             py::kw_only(), "index"_a, "t"_a)
        .def("arrow",
             py::overload_cast<double, bool>(&PolylineRuler::arrow, py::const_),
             "range"_a, //
             py::kw_only(), "smooth_joint"_a = true)
        .def("arrows",
             py::overload_cast<const Eigen::Ref<const Eigen::VectorXd> &, bool>(
                 &PolylineRuler::arrows, py::const_),
             "ranges"_a, //
             py::kw_only(), "smooth_joint"_a = true)
        .def("arrows",
             py::overload_cast<double, bool, bool>(&PolylineRuler::arrows,
                                                   py::const_),
             "step"_a, //
             py::kw_only(), "with_last"_a = true, "smooth_joint"_a = true)
        .def("scanline", &PolylineRuler::scanline, "range"_a, //
             py::kw_only(), "min"_a, "max"_a, "smooth_joint"_a = true)
        //
        .def("local_frame", &PolylineRuler::local_frame, "range"_a,
             py::kw_only(), "smooth_joint"_a = true)
        //
        .def_static(
            "_squareDistance",
            py::overload_cast<const Eigen::Vector3d &, const Eigen::Vector3d &,
                              bool>(&PolylineRuler::squareDistance),
            "a"_a, "b"_a, py::kw_only(), "is_wgs84"_a = false)
        .def_static(
            "_distance",
            py::overload_cast<const Eigen::Vector3d &, const Eigen::Vector3d &,
                              bool>(&PolylineRuler::distance),
            "a"_a, "b"_a, py::kw_only(), "is_wgs84"_a = false)
        .def_static(
            "_lineDistance",
            py::overload_cast<const Eigen::Ref<const RowVectors> &, bool>(
                &PolylineRuler::lineDistance),
            "line"_a, py::kw_only(), "is_wgs84"_a = false)
        .def("lineDistance",
             py::overload_cast<>(&PolylineRuler::lineDistance, py::const_))
        .def_static("_along",
                    py::overload_cast<const Eigen::Ref<const RowVectors> &,
                                      double, bool>(&PolylineRuler::along),
                    "line"_a, "dist"_a, py::kw_only(), "is_wgs84"_a = false)
        .def("along",
             py::overload_cast<double>(&PolylineRuler::along, py::const_),
             "dist"_a)
        //
        .def_static(
            "_pointToSegmentDistance",
            py::overload_cast<const Eigen::Vector3d &, const Eigen::Vector3d &,
                              const Eigen::Vector3d &, bool>(
                &PolylineRuler::pointToSegmentDistance),
            "P"_a, "A"_a, "B"_a, py::kw_only(), "is_wgs84"_a = false)
        .def_static("_pointOnLine",
                    py::overload_cast<const Eigen::Ref<const RowVectors> &,
                                      const Eigen::Vector3d &, bool>(
                        &PolylineRuler::pointOnLine),
                    "line"_a, "P"_a, py::kw_only(), "is_wgs84"_a = false)
        .def("pointOnLine",
             py::overload_cast<const Eigen::Vector3d &>(
                 &PolylineRuler::pointOnLine, py::const_),
             "P"_a)
        .def_static(
            "_lineSlice",
            py::overload_cast<const Eigen::Vector3d &, const Eigen::Vector3d &,
                              const Eigen::Ref<const RowVectors> &, bool>(
                &PolylineRuler::lineSlice),
            "start"_a, "stop"_a, "line"_a, //
            py::kw_only(), "is_wgs84"_a = false)
        .def(
            "lineSlice",
            py::overload_cast<const Eigen::Vector3d &, const Eigen::Vector3d &>(
                &PolylineRuler::lineSlice, py::const_),
            //
            "start"_a, "stop"_a)
        //
        .def_static(
            "_lineSliceAlong",
            py::overload_cast<double, double,
                              const Eigen::Ref<const RowVectors> &, bool>(
                &PolylineRuler::lineSliceAlong),
            "start"_a, "stop"_a, "line"_a, //
            py::kw_only(), "is_wgs84"_a = false)
        .def("lineSliceAlong",
             py::overload_cast<double, double>(&PolylineRuler::lineSliceAlong,
                                               py::const_),
             "start"_a, "stop"_a)
        .def_static("_interpolate", &PolylineRuler::interpolate, //
                    "A"_a, "B"_a, py::kw_only(), "t"_a, "is_wgs84"_a = false)
        //
        ;

    m.def("douglas_simplify",
          py::overload_cast<const RowVectors &, double, bool,
                            bool>(&douglas_simplify), //
          "coords"_a, "epsilon"_a,                    //
          py::kw_only(),                              //
          "is_wgs84"_a = false,                       //
          "recursive"_a = true);
    m.def(
        "douglas_simplify",
        py::overload_cast<const Eigen::Ref<const RowVectorsNx2> &, double, bool,
                          bool>(&douglas_simplify), //
        "coords"_a, "epsilon"_a,                    //
        py::kw_only(),                              //
        "is_wgs84"_a = false,                       //
        "recursive"_a = true);
    m.def("douglas_simplify_mask",
          py::overload_cast<const RowVectors &, double, bool,
                            bool>(&douglas_simplify_mask), //
          "coords"_a, "epsilon"_a, py::kw_only(),          //
          "is_wgs84"_a = false,                            //
          "recursive"_a = true);
    m.def(
        "douglas_simplify_mask",
        py::overload_cast<const Eigen::Ref<const RowVectorsNx2> &, double, bool,
                          bool>(&douglas_simplify_mask), //
        "coords"_a, "epsilon"_a, py::kw_only(),          //
        "is_wgs84"_a = false,                            //
        "recursive"_a = true);
    m.def("douglas_simplify_indexes",
          py::overload_cast<const RowVectors &, double, bool,
                            bool>(&douglas_simplify_indexes), //
          "coords"_a, "epsilon"_a, py::kw_only(),             //
          "is_wgs84"_a = false,                               //
          "recursive"_a = true);
    m.def(
        "douglas_simplify_indexes",
        py::overload_cast<const Eigen::Ref<const RowVectorsNx2> &, double, bool,
                          bool>(&douglas_simplify_indexes), //
        "coords"_a, "epsilon"_a, py::kw_only(),             //
        "is_wgs84"_a = false,                               //
        "recursive"_a = true);
}
} // namespace cubao
