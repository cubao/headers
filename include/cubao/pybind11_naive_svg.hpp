// should sync
// -
// https://github.com/cubao/pybind11-naive-svg/blob/master/src/pybind11_naive_svg.cpp
// -
// https://github.com/cubao/headers/tree/main/include/cubao/pybind11_naive_svg.hpp

#pragma once

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cubao_inline.hpp"
#include "naive_svg.hpp"

namespace cubao
{
namespace py = pybind11;
using namespace pybind11::literals;
using rvp = py::return_value_policy;

using RowVectorsNx2 = Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>;

CUBAO_INLINE void bind_naive_svg(py::module &m)
{
    // https://github.com/gagan-bansal/geojson2svg
    // https://milevski.co/geojson2svg/demo/lands.svg
    // 还是转化到 ENU 下，更好。radius 的尺度是一致的, stroke 也更好调

#define SETUP_FLUENT_API_PYBIND(Klass, VarType, VarName)                       \
    .def(#VarName, [](const Klass &self) { return self.VarName(); })           \
        .def(                                                                  \
            #VarName,                                                          \
            [](Klass &self, const VarType &v) -> Klass & {                     \
                return self.VarName(v);                                        \
            },                                                                 \
            rvp::reference_internal)

    using Color = SVG::Color;
    py::class_<Color>(m, "Color",
                      py::module_local()) //
        .def(py::init<int>(), "rgb"_a = -1, "Initialize Color with RGB value")
        .def(py::init<int, int, int, float>(), "r"_a, "g"_a, "b"_a,
             "a"_a = -1.f,
             "Initialize Color with R, G, B, and optional Alpha values")
            SETUP_FLUENT_API_PYBIND(Color, int, r)
                SETUP_FLUENT_API_PYBIND(Color, int, g)
                    SETUP_FLUENT_API_PYBIND(Color, int, b)
                        SETUP_FLUENT_API_PYBIND(Color, float, a)
        .def("invalid", &Color::invalid, "Check if the color is invalid")
        .def("to_string", &Color::to_string,
             "Convert color to string representation")

        .def("clone", &Color::clone, "Create a deep copy of the Color object")
        .def(
            "__copy__",
            [](const Color &self, py::dict) -> Color {
                // always deepcopy (maybe not good?)
                return self.clone();
            },
            "Create a shallow copy of the Color object")
        .def(
            "__deepcopy__",
            [](const Color &self, py::dict) -> Color { return self.clone(); },
            "memo"_a, "Create a deep copy of the Color object")
        .def_static(
            "parse",
            [](const std::string &text) {
                int i = text.size() - 6;
                if (i < 0) {
                    return Color();
                }
                return Color(std::stoi(text.substr(i, 2), nullptr, 16),
                             std::stoi(text.substr(i + 2, 2), nullptr, 16),
                             std::stoi(text.substr(i + 4, 2), nullptr, 16));
            },
            "Parse a color from a string representation")

        //
        .def(
            "__repr__", [](const Color &self) { return self.to_string(); },
            "Return a string representation of the Color object")
        //
        ;

    using Polyline = SVG::Polyline;
    py::class_<Polyline>(m, "Polyline",
                         py::module_local()) //
        .def(py::init([](const Eigen::Ref<const RowVectorsNx2> &points) {
                 std::vector<SVG::PointType> _(points.rows());
                 Eigen::Map<RowVectorsNx2>(&_[0][0], points.rows(), 2) = points;
                 return new SVG::Polyline(_);
             }),
             "points"_a, "Initialize Polyline with a set of points")
        //
        .def(
            "to_numpy",
            [](const Polyline &self) -> RowVectorsNx2 {
                auto &points = self.points();
                return Eigen::Map<const RowVectorsNx2>(&points[0][0],
                                                       points.size(), 2);
            },
            "Convert Polyline points to NumPy array")
        .def(
            "from_numpy",
            [](Polyline &self,
               const Eigen::Ref<const RowVectorsNx2> &points) -> Polyline & {
                std::vector<SVG::PointType> _(points.rows());
                Eigen::Map<RowVectorsNx2>(&_[0][0], points.rows(), 2) = points;
                return self.points(_);
            },
            rvp::reference_internal,
            "Set Polyline points from NumPy array") //
        //
        SETUP_FLUENT_API_PYBIND(Polyline, Color, stroke)
            SETUP_FLUENT_API_PYBIND(Polyline, double, stroke_width)
                SETUP_FLUENT_API_PYBIND(Polyline, Color, fill)
                    SETUP_FLUENT_API_PYBIND(Polyline, std::string, attrs)
        //
        .def("to_string", &Polyline::to_string,
             "Convert Polyline to SVG string representation")
        .def("clone", &Polyline::clone,
             "Create a deep copy of the Polyline object")
        .def(
            "__copy__",
            [](const Polyline &self, py::dict) -> Polyline {
                return self.clone();
            },
            "Create a shallow copy of the Polyline object")
        .def(
            "__deepcopy__",
            [](const Polyline &self, py::dict) -> Polyline {
                return self.clone();
            },
            "memo"_a, "Create a deep copy of the Polyline object")
        //
        ;

    using Polygon = SVG::Polygon;
    py::class_<Polygon>(m, "Polygon", py::module_local()) //
                                                          //
        .def(py::init([](const Eigen::Ref<const RowVectorsNx2> &points) {
                 std::vector<SVG::PointType> _(points.rows());
                 Eigen::Map<RowVectorsNx2>(&_[0][0], points.rows(), 2) = points;
                 return new SVG::Polygon(_);
             }),
             "points"_a,
             "Initialize Polygon with a set of points") //
        //
        .def(
            "to_numpy",
            [](const Polygon &self) -> RowVectorsNx2 {
                auto &points = self.points();
                return Eigen::Map<const RowVectorsNx2>(&points[0][0],
                                                       points.size(), 2);
            },
            "Convert Polygon points to NumPy array")
        .def(
            "from_numpy",
            [](Polygon &self,
               const Eigen::Ref<const RowVectorsNx2> &points) -> Polygon & {
                std::vector<SVG::PointType> _(points.rows());
                Eigen::Map<RowVectorsNx2>(&_[0][0], points.rows(), 2) = points;
                return self.points(_);
            },
            rvp::reference_internal,
            "Set Polygon points from NumPy array") //
        SETUP_FLUENT_API_PYBIND(Polygon, Color, stroke)
            SETUP_FLUENT_API_PYBIND(Polygon, double, stroke_width)
                SETUP_FLUENT_API_PYBIND(Polygon, Color, fill)
                    SETUP_FLUENT_API_PYBIND(Polygon, std::string, attrs)
        //
        .def("to_string", &Polygon::to_string,
             "Convert Polygon to SVG string representation")
        .def("clone", &Polygon::clone,
             "Create a deep copy of the Polygon object")
        .def(
            "__copy__",
            [](const Polygon &self, py::dict) -> Polygon {
                return self.clone();
            },
            "Create a shallow copy of the Polygon object")
        .def(
            "__deepcopy__",
            [](const Polygon &self, py::dict) -> Polygon {
                return self.clone();
            },
            "memo"_a, "Create a deep copy of the Polygon object")
        //

        ;

    using Circle = SVG::Circle;
    py::class_<Circle>(m, "Circle", py::module_local()) //
        .def(py::init([](const Eigen::Vector2d &center, double r) {
                 return new SVG::Circle({center[0], center[1]}, r);
             }),
             "center"_a, "r"_a = 1.0,
             "Initialize Circle with center point and radius") //

        .def(
            "center",
            [](const Circle &self) -> Eigen::Vector2d {
                auto &c = self.center();
                return Eigen::Vector2d(c[0], c[1]);
            },
            "Get the center of the Circle")
        .def(
            "center",
            [](Circle &self, const Eigen::Vector2d &center) -> Circle & {
                return self.center({center[0], center[1]});
            },
            rvp::reference_internal, "Set the center of the Circle")
            SETUP_FLUENT_API_PYBIND(Circle, double, x)
                SETUP_FLUENT_API_PYBIND(Circle, double, y)
                    SETUP_FLUENT_API_PYBIND(Circle, double, r)
        //
        SETUP_FLUENT_API_PYBIND(Circle, Color, stroke)
            SETUP_FLUENT_API_PYBIND(Circle, double, stroke_width)
                SETUP_FLUENT_API_PYBIND(Circle, Color, fill)
                    SETUP_FLUENT_API_PYBIND(Circle, std::string, attrs)
        //
        .def("to_string", &Circle::to_string,
             "Convert Circle to SVG string representation")
        .def("clone", &Circle::clone, "Create a deep copy of the Circle object")
        .def(
            "__copy__",
            [](const Circle &self, py::dict) -> Circle { return self.clone(); },
            "Create a shallow copy of the Circle object")
        .def(
            "__deepcopy__",
            [](const Circle &self, py::dict) -> Circle { return self.clone(); },
            "memo"_a, "Create a deep copy of the Circle object")
        //
        ;

    using Text = SVG::Text;
    py::class_<Text>(m, "Text", py::module_local()) //
        .def(py::init([](const Eigen::Vector2d &position,
                         const std::string &text, double fontsize) {
                 return new SVG::Text({position[0], position[1]}, text,
                                      fontsize);
             }),
             "position"_a, "text"_a, "fontsize"_a = 10.0,
             "Initialize Text with position, content, and font size") //
                                                                      //
        .def(
            "position",
            [](const Text &self) -> Eigen::Vector2d {
                auto &p = self.position();
                return Eigen::Vector2d(p[0], p[1]);
            },
            "Get the position of the Text")
        .def(
            "position",
            [](Text &self, const Eigen::Vector2d &position) -> Text & {
                return self.position({position[0], position[1]});
            },
            rvp::reference_internal, "Set the position of the Text")
        //
        SETUP_FLUENT_API_PYBIND(Text, std::string, text)
            SETUP_FLUENT_API_PYBIND(Text, std::vector<std::string>, lines)
                SETUP_FLUENT_API_PYBIND(Text, double, fontsize)
        //
        SETUP_FLUENT_API_PYBIND(Text, Color, stroke)
            SETUP_FLUENT_API_PYBIND(Text, double, stroke_width)
                SETUP_FLUENT_API_PYBIND(Text, Color, fill)
                    SETUP_FLUENT_API_PYBIND(Text, std::string, attrs)
        //
        .def("to_string", &Text::to_string,
             "Convert Text to SVG string representation")
        .def("clone", &Text::clone, "Create a deep copy of the Text object")
        .def(
            "__copy__",
            [](const Text &self, py::dict) -> Text { return self.clone(); },
            "Create a shallow copy of the Text object")
        .def(
            "__deepcopy__",
            [](const Text &self, py::dict) -> Text { return self.clone(); },
            "memo"_a, "Create a deep copy of the Text object")
        //
        .def_static("html_escape", &Text::html_escape, "text"_a,
                    "Escape special characters in the text for HTML")
        //
        ;

    py::class_<SVG>(m, "SVG", py::module_local())
        .def(py::init<double, double>(), "width"_a, "height"_a,
             "Initialize SVG with width and height")
        //
        .def("clone", &SVG::clone, "Create a deep copy of the SVG object")
        .def(
            "__copy__", [](const SVG &self, py::dict) { return self.clone(); },
            "Create a shallow copy of the SVG object")
        .def(
            "__deepcopy__",
            [](const SVG &self, py::dict) { return self.clone(); }, "memo"_a,
            "Create a deep copy of the SVG object")
        //
        SETUP_FLUENT_API_PYBIND(SVG, double, width)                 //
        SETUP_FLUENT_API_PYBIND(SVG, double, height)                //
        SETUP_FLUENT_API_PYBIND(SVG, std::vector<double>, view_box) //
        SETUP_FLUENT_API_PYBIND(SVG, double, grid_step)             //
        SETUP_FLUENT_API_PYBIND(SVG, std::vector<double>, grid_x)   //
        SETUP_FLUENT_API_PYBIND(SVG, std::vector<double>, grid_y)   //
        SETUP_FLUENT_API_PYBIND(SVG, Color, grid_color)             //
        SETUP_FLUENT_API_PYBIND(SVG, Color, background)             //
        SETUP_FLUENT_API_PYBIND(SVG, std::string, attrs)            //
                                                                    //
        .def("add", py::overload_cast<const Polyline &>(&SVG::add),
             "polyline"_a, rvp::reference_internal, "Add a Polyline to the SVG")
        .def("add", py::overload_cast<const Polygon &>(&SVG::add), //
             "polygon"_a, rvp::reference_internal, "Add a Polygon to the SVG")
        .def("add", py::overload_cast<const Circle &>(&SVG::add), //
             "circle"_a, rvp ::reference_internal, "Add a Circle to the SVG")
        .def("add", py::overload_cast<const Text &>(&SVG::add), //
             "text"_a, rvp::reference_internal, "Add a Text to the SVG")
        //
        .def(
            "add_polyline",
            [](SVG &self,
               const Eigen::Ref<const RowVectorsNx2> &points) -> Polyline & {
                std::vector<SVG::PointType> _(points.rows());
                Eigen::Map<RowVectorsNx2>(&_[0][0], points.rows(), 2) = points;
                return self.add_polyline(_);
            },
            "points"_a, rvp::reference_internal,
            "Add a Polyline to the SVG using NumPy array of points")
        .def(
            "add_polygon",
            [](SVG &self,
               const Eigen::Ref<const RowVectorsNx2> &points) -> Polygon & {
                std::vector<SVG::PointType> _(points.rows());
                Eigen::Map<RowVectorsNx2>(&_[0][0], points.rows(), 2) = points;
                return self.add_polygon(_);
            },
            "points"_a, rvp::reference_internal,
            "Add a Polygon to the SVG using NumPy array of points")
        .def(
            "add_circle",
            [](SVG &self, const Eigen::Vector2d &center, double r) -> Circle & {
                return self.add_circle({center[0], center[1]}, r);
            },
            "center"_a, py::kw_only(), "r"_a = 1.0, rvp::reference_internal,
            "Add a Circle to the SVG")
        .def(
            "add_text",
            [](SVG &self, const Eigen::Vector2d &position,
               const std::string &text, double fontsize) -> Text & {
                return self.add_text({position[0], position[1]}, text,
                                     fontsize);
            },
            "position"_a, py::kw_only(), "text"_a, "fontsize"_a = 10.0,
            rvp::reference_internal, "Add a Text to the SVG")
        //
        .def("num_elements", &SVG::num_elements,
             "Get the number of elements in the SVG")
        .def("empty", &SVG::empty, "Check if the SVG is empty")
        .def("pop", &SVG::pop, "Remove the last added element from the SVG")
        //
        .def("is_polyline", &SVG::is_polyline,
             "Check if the element at the given index is a Polyline")
        .def("is_polygon", &SVG::is_polygon,
             "Check if the element at the given index is a Polygon")
        .def("is_circle", &SVG::is_circle,
             "Check if the element at the given index is a Circle")
        .def("is_text", &SVG::is_text,
             "Check if the element at the given index is a Text")
        //
        .def("as_polyline", py::overload_cast<int>(&SVG::as_polyline),
             "index"_a, rvp::reference_internal,
             "Get the element at the given index as a Polyline")
        .def("as_polygon", py::overload_cast<int>(&SVG::as_polygon), "index"_a,
             rvp::reference_internal,
             "Get the element at the given index as a Polygon")
        .def("as_circle", py::overload_cast<int>(&SVG::as_circle), "index"_a,
             rvp::reference_internal,
             "Get the element at the given index as a Circle")
        .def("as_text", py::overload_cast<int>(&SVG::as_text), "index"_a,
             rvp::reference_internal,
             "Get the element at the given index as a Text")
        //
        .def("to_string", &SVG::to_string,
             "Convert the SVG to a string representation")
        .def("dump", &SVG::dump, "path"_a, "Save the SVG to a file")
        //
        ;
}
} // namespace cubao
