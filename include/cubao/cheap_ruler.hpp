#ifndef CUBAO_CHEAP_RULER_HPP
#define CUBAO_CHEAP_RULER_HPP

// should sync
// - https://github.com/cubao/polyline-ruler/blob/master/src/cheap_ruler.hpp
// - https://github.com/cubao/headers/tree/main/include/cubao/cheap_ruler.hpp

// based on https://github.com/mapbox/cheap-ruler-cpp

// https://github.com/microsoft/vscode-cpptools/issues/9692
#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <Eigen/Core>

#include <cassert>
#include <cstdint>
#include <limits>
#include <array>
#include <tuple>
#include <utility>
#include <vector>

#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace cubao
{
using RowVectors = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

namespace cheap_ruler
{

using box = std::pair<Eigen::Vector3d, Eigen::Vector3d>;
using line_string = RowVectors;
using linear_ring = RowVectors;
using multi_line_string = RowVectors;
using point = Eigen::Vector3d;
using polygon = RowVectors;

class CheapRuler
{

  public:
    // Values that define WGS84 ellipsoid model of the Earth
    static constexpr double RE = 6378.137;            // equatorial radius
    static constexpr double FE = 1.0 / 298.257223563; // flattening

    static constexpr double E2 = FE * (2 - FE);
    static constexpr double RAD = M_PI / 180.0;

    enum Unit
    {
        Kilometers,
        Miles,
        NauticalMiles,
        Meters,
        Metres = Meters,
        Yards,
        Feet,
        Inches
    };

    // k := [kx, ky, kz]
    // delta(lon) * kx -> delta east/west in your metric/unit (e.g. meters)
    // delta(lat) * ky -> delta north/south in your metric/unit (e.g. meters)
    // delta(alt) * kz -> delta up/down in your metric/unit (e.g. meters)
    Eigen::Vector3d k() const { return Eigen::Vector3d(kx, ky, kz); }

    static Eigen::Vector3d k(double latitude, Unit unit = Meters)
    {
        return CheapRuler(latitude, unit).k();
    }

    //
    // A collection of very fast approximations to common geodesic measurements.
    // Useful for performance-sensitive code that measures things on a city
    // scale. Point coordinates are in the [x = longitude, y = latitude] form.
    //
    explicit CheapRuler(double latitude, Unit unit = Meters)
    {
        double m = 0.;

        switch (unit) {
        case Kilometers:
            m = 1.;
            break;
        case Miles:
            m = 1000. / 1609.344;
            break;
        case NauticalMiles:
            m = 1000. / 1852.;
            break;
        case Meters:
            m = 1000.;
            break;
        case Yards:
            m = 1000. / 0.9144;
            break;
        case Feet:
            m = 1000. / 0.3048;
            break;
        case Inches:
            m = 1000. / 0.0254;
            break;
        }

        // Curvature formulas from
        // https://en.wikipedia.org/wiki/Earth_radius#Meridional
        double mul = RAD * RE * m;
        double coslat = std::cos(latitude * RAD);
        double w2 = 1 / (1 - E2 * (1 - coslat * coslat));
        double w = std::sqrt(w2);

        // multipliers for converting longitude and latitude degrees into
        // distance
        kx = mul * w * coslat;        // based on normal radius of curvature
        ky = mul * w * w2 * (1 - E2); // based on meridonal radius of curvature
        kz = m / 1000.0;              // altitute always in meter
    }

    static CheapRuler fromTile(uint32_t y, uint32_t z)
    {
        assert(z < 32);
        double n = M_PI * (1. - 2. * (y + 0.5) / double(uint32_t(1) << z));
        double latitude = std::atan(std::sinh(n)) / RAD;

        return CheapRuler(latitude);
    }

    // Note that:
    //      should always pass by const reference for eigen types!!!
    //      https://eigen.tuxfamily.org/dox/group__TopicPassingByValue.html

    point delta(const point &lla0, const point &lla1) const
    {
        auto dx = longDiff(lla1[0], lla0[0]) * kx;
        auto dy = (lla1[1] - lla0[1]) * ky;
        auto dz = (lla1[2] - lla0[2]) * kz;
        return point(dx, dy, dz);
    }

    double squareDistance(const point &a, const point &b) const
    {
        auto dx = longDiff(a[0], b[0]) * kx;
        auto dy = (a[1] - b[1]) * ky;
        auto dz = (a[2] - b[2]) * kz;
        return dx * dx + dy * dy + dz * dz;
    }

    //
    // Given two points of the form [x = longitude, y = latitude], returns the
    // distance.
    //
    double distance(const point &a, const point &b) const
    {
        return std::sqrt(squareDistance(a, b));
    }

    //
    // Returns the bearing between two points in angles.
    //        0
    //   -45         45
    // -90              90
    // -135         135
    //       180
    double bearing(const point &a, const point &b) const
    {
        auto dx = longDiff(b[0], a[0]) * kx;
        auto dy = (b[1] - a[1]) * ky;

        return std::atan2(dx, dy) / RAD;
    }

    //
    // Returns a new point given distance and bearing from the starting point.
    //
    point destination(const point &origin, double dist, double bearing) const
    {
        auto a = bearing * RAD;

        return offset(origin, std::sin(a) * dist, std::cos(a) * dist);
    }

    //
    // Returns a new point given easting and northing offsets from the starting
    // point.
    //
    point offset(const point &origin, double dx, double dy, double dz = 0) const
    {
        return origin + point(dx / kx, dy / ky, dz / kz);
    }

    //
    // Given a line (an array of points), returns the total line distance.
    //
    double lineDistance(const Eigen::Ref<const line_string> &points)
    {
        double total = 0.;

        for (int i = 1; i < points.rows(); ++i) {
            total += distance(points.row(i - 1), points.row(i));
        }

        return total;
    }

    //
    // Given a polygon (an array of rings, where each ring is an array of
    // points), returns the area.
    //
    double area(const Eigen::Ref<const polygon> &ring) const
    {
        double sum = 0.;

        for (unsigned j = 0, len = ring.rows(), k = len - 1; j < len; k = j++) {
            sum += longDiff(ring(j, 0), ring(k, 0)) * (ring(j, 1) + ring(k, 1));
        }

        return (std::abs(sum) / 2.) * kx * ky;
    }

    //
    // Returns the point at a specified distance along the line.
    //
    point along(const Eigen::Ref<const line_string> &line, double dist) const
    {
        double sum = 0.;

        if (!line.rows()) {
            return {};
        }

        if (dist <= 0.) {
            return line.row(0);
        }

        for (unsigned i = 0; i < line.rows() - 1; ++i) {
            auto p0 = line.row(i);
            auto p1 = line.row(i + 1);
            auto d = distance(p0, p1);

            sum += d;

            if (sum > dist) {
                return interpolate(p0, p1, (dist - (sum - d)) / d);
            }
        }

        return line.row(line.rows() - 1);
    }

    //
    // Returns the distance from a point `p` to a line segment `a` to `b`.
    //
    double pointToSegmentDistance(const point &p, const point &a,
                                  const point &b) const
    {
        auto t = 0.0;
        auto x = a[0];
        auto y = a[1];
        auto z = a[2];
        auto dx = longDiff(b[0], x) * kx;
        auto dy = (b[1] - y) * ky;
        auto dz = (b[2] - z) * kz;

        if (dx != 0.0 || dy != 0.0 || dz != 0.0) {
            // t = DOT(ap, ab)
            t = (longDiff(p[0], x) * kx * dx + (p[1] - y) * ky * dy +
                 (p[2] - z) * kz * dz) /
                (dx * dx + dy * dy + dz * dz);
            if (t > 1.0) {
                x = b[0];
                y = b[1];
                z = b[2];
            } else if (t > 0.0) {
                x += (dx / kx) * t;
                y += (dy / ky) * t;
                z += (dz / kz) * t;
            }
        }
        return distance(p, {x, y, z});
    }

    //
    // Returns a tuple of the form <point, index, t> where point is closest
    // point on the line from the given point, index is the start index of the
    // segment with the closest point, and t is a parameter from 0 to 1 that
    // indicates where the closest point is on that segment.
    //
    std::tuple<point, int, double>
    pointOnLine(const Eigen::Ref<const line_string> &line, const point &p) const
    {
        double minDist = std::numeric_limits<double>::infinity();
        double minX = 0., minY = 0., minZ = 0, minI = 0., minT = 0.;

        if (!line.rows()) {
            return std::make_tuple(p, -1, 0.);
        }

        for (unsigned i = 0; i < line.rows() - 1; ++i) {
            auto t = 0.;
            auto x = line(i, 0);
            auto y = line(i, 1);
            auto z = line(i, 2);
            auto dx = longDiff(line(i + 1, 0), x) * kx;
            auto dy = (line(i + 1, 1) - y) * ky;
            auto dz = (line(i + 1, 2) - z) * kz;

            if (dx != 0. || dy != 0. || dz != 0.) {
                // t = DOT(ap, ab)
                t = (longDiff(p[0], x) * kx * dx + (p[1] - y) * ky * dy +
                     (p[2] - z) * kz * dz) /
                    (dx * dx + dy * dy + dz * dz);
                if (t > 1) {
                    x = line(i + 1, 0);
                    y = line(i + 1, 1);
                    z = line(i + 1, 2);
                } else if (t > 0) {
                    x += (dx / kx) * t;
                    y += (dy / ky) * t;
                    z += (dz / kz) * t;
                }
            }

            auto sqDist = squareDistance(p, {x, y, z});

            if (sqDist < minDist) {
                minDist = sqDist;
                minX = x;
                minY = y;
                minZ = z;
                minI = i;
                minT = t;
            }
        }

        return std::make_tuple(point(minX, minY, minZ), minI,
                               ::fmax(0., ::fmin(1., minT)));
    }

    //
    // Returns a part of the given line between the start and the stop points
    // (or their closest points on the line).
    //
    line_string lineSlice(const point &start, const point &stop,
                          const Eigen::Ref<const line_string> &line) const
    {
        auto getPoint = [](auto &tuple) -> const Eigen::Vector3d & {
            return std::get<0>(tuple);
        };
        auto getIndex = [](auto &tuple) -> int { return std::get<1>(tuple); };
        auto getT = [](auto &tuple) -> double { return std::get<2>(tuple); };
        auto same_point = [](const Eigen::Vector3d &p1,
                             const Eigen::Vector3d &p2) { return p1 == p2; };

        auto p1 = pointOnLine(line, start);
        auto p2 = pointOnLine(line, stop);

        if (getIndex(p1) > getIndex(p2) ||
            (getIndex(p1) == getIndex(p2) && getT(p1) > getT(p2))) {
            auto tmp = p1;
            p1 = p2;
            p2 = tmp;
        }

        auto slice = std::vector<Eigen::Vector3d>{getPoint(p1)};

        auto l = getIndex(p1) + 1;
        auto r = getIndex(p2);

        if (!same_point(line.row(l), slice[0]) && l <= r) {
            slice.push_back(line.row(l));
        }

        for (int i = l + 1; i <= r; ++i) {
            slice.push_back(line.row(i));
        }

        if (!same_point(line.row(r), getPoint(p2))) {
            slice.push_back(getPoint(p2));
        }

        return line_string::Map(slice[0].data(), (Eigen::Index)slice.size(), 3);
    }

    //
    // Returns a part of the given line between the start and the stop points
    // indicated by distance along the line.
    //
    line_string lineSliceAlong(double start, double stop,
                               const Eigen::Ref<line_string> &line) const
    {
        double sum = 0.;
        std::vector<Eigen::Vector3d> slice;

        for (int i = 1; i < line.rows(); ++i) {
            auto p0 = line.row(i - 1);
            auto p1 = line.row(i);
            auto d = distance(p0, p1);

            sum += d;

            if (sum > start && slice.size() == 0) {
                slice.push_back(interpolate(p0, p1, (start - (sum - d)) / d));
            }

            if (sum >= stop) {
                slice.push_back(interpolate(p0, p1, (stop - (sum - d)) / d));
                return line_string::Map(slice[0].data(),
                                        (Eigen::Index)slice.size(), 3);
            }

            if (sum > start) {
                slice.push_back(p1);
            }
        }

        return line_string::Map(slice[0].data(), (Eigen::Index)slice.size(), 3);
    }

    //
    // Given a point, returns a bounding box object ([w, s, e, n])
    // created from the given point buffered by a given distance.
    //
    box bufferPoint(const point &p, double buffer) const
    {
        auto v = buffer / ky;
        auto h = buffer / kx;
        auto e = buffer / kz;

        return box(point(p[0] - h, p[1] - v, p[2] - e),
                   point(p[0] + h, p[1] + v, p[2] + e));
    }

    //
    // Given a bounding box, returns the box buffered by a given distance.
    //
    box bufferBBox(const box &bbox, double buffer) const
    {
        auto v = buffer / ky;
        auto h = buffer / kx;
        auto e = buffer / kz;

        return box(
            point(bbox.first[0] - h, bbox.first[1] - v, bbox.first[2] - e),
            point(bbox.second[0] + h, bbox.second[1] + v, bbox.second[2] + e));
    }

    //
    // Returns true if the given point is inside in the given bounding box,
    // otherwise false.
    //
    static bool insideBBox(const point &p, const box &bbox,
                           bool check_z = false)
    {
        bool inside2d = p[1] >= bbox.first[1] && p[1] <= bbox.second[1] &&
                        longDiff(p[0], bbox.first[0]) >= 0 &&
                        longDiff(p[0], bbox.second[0]) <= 0;
        if (!check_z) {
            return inside2d;
        }
        return inside2d && bbox.first[2] <= p[2] && p[2] <= bbox.second[2];
    }

    static point interpolate(const point &a, const point &b, double t)
    {
        double dx = longDiff(b[0], a[0]);
        double dy = b[1] - a[1];
        double dz = b[2] - a[2];

        return point(a[0] + dx * t, a[1] + dy * t, a[2] + dz * t);
    }

    static double longDiff(double a, double b)
    {
        return std::remainder(a - b, 360);
    }

  private:
    double ky;
    double kx;
    double kz;
};

} // namespace cheap_ruler

using CheapRuler = cheap_ruler::CheapRuler;
} // namespace cubao

#endif
