#ifndef CUBAO_CRS_TRANSFORM_HPP
#define CUBAO_CRS_TRANSFORM_HPP

// should sync
// - https://github.com/cubao/polyline-ruler/blob/master/src/crs_transform.hpp
// - https://github.com/cubao/headers/tree/main/include/cubao/crs_transform.hpp

// https://github.com/microsoft/vscode-cpptools/issues/9692
#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <optional>

#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace cubao
{

namespace internal
{
// https://github.com/planet36/ecef-geodetic/blob/main/ecef-to-geodetic-funcs.hpp#L3649-L3735
inline void ecef_to_geodetic(const double x, const double y, const double z,
                             double &lon_rad, double &lat_rad, double &ht)
{
    const auto w2 = x * x + y * y;
    const auto w = std::sqrt(w2);
    const auto z2 = z * z;
    lon_rad = std::atan2(y, x);

    constexpr auto a = 6378137.0;
    constexpr auto e2 = 6.6943799901377997e-3;
    constexpr auto a1 = a * e2;
    constexpr auto a2 = a1 * a1;
    constexpr auto a3 = a1 * e2 / 2;
    constexpr auto a4 = 2.5 * a2;
    constexpr auto a5 = a1 + a3;

    const auto r2 = w2 + z2;
    const auto r = std::sqrt(r2);

    const auto s2 = z2 / r2;
    const auto c2 = w2 / r2;
    auto u = a2 / r;
    auto v = a3 - a4 / r;

    double s;
    double c;
    double ss;

    // cos(45 deg)^2 == 0.5
    if (c2 > 0.5) // Equatorial
    {
        s = (z / r) * (1 + c2 * (a1 + u + s2 * v) / r);
        lat_rad = std::asin(s);
        ss = s * s;
        c = std::sqrt(1 - ss);
    } else // Polar
    {
        c = (w / r) * (1 - s2 * (a5 - u - c2 * v) / r);
        lat_rad = std::acos(c);
        ss = 1 - c * c;
        s = std::sqrt(ss);

        if (z < 0) {
            lat_rad = -lat_rad;
            s = -s;
        }
    }

    const auto d2 = 1 - e2 * ss;
    const auto Rn = a / std::sqrt(d2);
    const auto Rm = (1 - e2) * Rn / d2;
    const auto rf = (1 - e2) * Rn;
    u = w - Rn * c;
    v = z - rf * s;
    const auto f = c * u + s * v;
    const auto m = c * v - s * u;
    const auto p = m / (Rm + f);

    lat_rad += p;

    ht = f + m * p / 2;
}

inline void geodetic_to_ecef(const double lon_rad, const double lat_rad,
                             const double ht, //
                             double &x, double &y, double &z)
{
    constexpr double a = 6378137.0;
    constexpr double b = 6356752.314245;
    constexpr double E = (a * a - b * b) / (a * a);
    double coslat = std::cos(lat_rad);
    double sinlat = std::sin(lat_rad);
    double coslong = std::cos(lon_rad);
    double sinlong = std::sin(lon_rad);
    double N = a / (std::sqrt(1 - E * sinlat * sinlat));
    double NH = N + ht;
    x = NH * coslat * coslong;
    y = NH * coslat * sinlong;
    z = (b * b * N / (a * a) + ht) * sinlat;
}
} // namespace internal

// note that:
//  lla: lon, lat, alt, in degrees

inline double to_degrees(double r) { return 180.0 / M_PI * r; }
inline double to_radians(double d) { return M_PI / 180.0 * d; }

inline Eigen::Vector3d ecef2lla(double x, double y, double z)
{
    double lon_rad, lat_rad, ht;
    internal::ecef_to_geodetic(x, y, z, lon_rad, lat_rad, ht);
    return {to_degrees(lon_rad), to_degrees(lat_rad), ht};
}
inline Eigen::Vector3d lla2ecef(double lon, double lat, double alt)
{
    Eigen::Vector3d xyz;
    internal::geodetic_to_ecef(to_radians(lon), to_radians(lat), alt, //
                               xyz[0], xyz[1], xyz[2]);
    return xyz;
}

using RowVectors = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
inline RowVectors ecef2lla(const Eigen::Ref<const RowVectors> &ecefs)
{
    const int N = ecefs.rows();
    if (!N) {
        return RowVectors(0, 3);
    }
    RowVectors ret = ecefs;
    for (int i = 0; i < N; ++i) {
        internal::ecef_to_geodetic(ret(i, 0), ret(i, 1), ret(i, 2), //
                                   ret(i, 0), ret(i, 1), ret(i, 2));
    }
    ret.col(0) *= 180.0 / M_PI;
    ret.col(1) *= 180.0 / M_PI;
    return ret;
}
inline RowVectors lla2ecef(const Eigen::Ref<const RowVectors> &llas)
{
    const int N = llas.rows();
    if (!N) {
        return RowVectors(0, 3);
    }
    RowVectors ret = llas;
    ret.col(0) *= M_PI / 180.0;
    ret.col(1) *= M_PI / 180.0;
    for (int i = 0; i < N; ++i) {
        internal::geodetic_to_ecef(ret(i, 0), ret(i, 1), ret(i, 2), //
                                   ret(i, 0), ret(i, 1), ret(i, 2));
    }
    return ret;
}

inline Eigen::Matrix3d R_ecef_enu(double lon, double lat)
{
    lon *= M_PI / 180.0;
    lat *= M_PI / 180.0;
    double s1 = std::sin(M_PI / 4.0 + lon / 2.0);
    double c1 = std::cos(M_PI / 4.0 + lon / 2.0);
    double s2 = std::sin(M_PI / 4.0 - lat / 2.0);
    double c2 = std::cos(M_PI / 4.0 - lat / 2.0);
    return Eigen::Quaterniond{c1 * c2, -c1 * s2, -s1 * s2, -s1 * c2}
        .toRotationMatrix()
        .transpose();
}

// T_ecef_enu: Transform ecef <-- enu
inline Eigen::Matrix4d T_ecef_enu(double lon, double lat, double alt)
{
    Eigen::Matrix4d T;
    T.setIdentity();
    T.topLeftCorner(3, 3) = R_ecef_enu(lon, lat);
    T.topRightCorner(3, 1) = lla2ecef(lon, lat, alt);
    return T;
}
inline Eigen::Matrix4d T_ecef_enu(const Eigen::Vector3d &lla)
{
    return T_ecef_enu(lla[0], lla[1], lla[2]);
}

inline RowVectors apply_transform(const Eigen::Matrix4d &T,
                                  const Eigen::Ref<const RowVectors> &points)
{
    return ((T.topLeftCorner<3, 3>() * points.transpose()).colwise() +
            T.topRightCorner<3, 1>())
        .transpose();
}

inline void apply_transform_inplace(const Eigen::Matrix4d &T,
                                    Eigen::Ref<RowVectors> points,
                                    int batch_size = 1000)
{
    int R = points.rows(), r = 0;
    while (r < R) {
        batch_size = std::min(batch_size, R - r);
        points.middleRows(r, batch_size) =
            apply_transform(T, points.middleRows(r, batch_size));
        r += batch_size;
    }
}

inline Eigen::Vector3d cheap_ruler_k(double latitude)
{
    // based on https://github.com/mapbox/cheap-ruler-cpp
    static constexpr double RE = 6378.137;
    static constexpr double FE = 1.0 / 298.257223563;
    static constexpr double E2 = FE * (2 - FE);
    static constexpr double RAD = M_PI / 180.0;
    static constexpr double MUL = RAD * RE * 1000.;
    double coslat = std::cos(latitude * RAD);
    double w2 = 1 / (1 - E2 * (1 - coslat * coslat));
    double w = std::sqrt(w2);
    return Eigen::Vector3d(MUL * w * coslat, MUL * w * w2 * (1 - E2), 1.0);
}

inline RowVectors lla2enu(const Eigen::Ref<const RowVectors> &llas,
                          std::optional<Eigen::Vector3d> anchor_lla = {},
                          bool cheap_ruler = true)
{
    if (!llas.rows()) {
        return RowVectors(0, 3);
    }
    if (!anchor_lla) {
        anchor_lla = llas.row(0);
    }
    if (!cheap_ruler) {
        return apply_transform(T_ecef_enu(*anchor_lla).inverse(),
                               lla2ecef(llas));
    }
    auto k = cheap_ruler_k((*anchor_lla)[1]);
    RowVectors enus = llas;
    for (int i = 0; i < 3; ++i) {
        enus.col(i).array() -= (*anchor_lla)[i];
        enus.col(i).array() *= k[i];
    }
    return enus;
}
inline RowVectors enu2lla(const Eigen::Ref<const RowVectors> &enus,
                          const Eigen::Vector3d &anchor_lla,
                          bool cheap_ruler = true)
{
    if (!enus.rows()) {
        return RowVectors(0, 3);
    }
    if (!cheap_ruler) {
        return ecef2lla(apply_transform(T_ecef_enu(anchor_lla), enus));
    }
    auto k = cheap_ruler_k(anchor_lla[1]);
    RowVectors llas = enus;
    for (int i = 0; i < 3; ++i) {
        llas.col(i).array() /= k[i];
        llas.col(i).array() += anchor_lla[i];
    }
    return llas;
}

inline RowVectors enu2ecef(const Eigen::Ref<const RowVectors> &enus,
                           const Eigen::Vector3d &anchor_lla,
                           bool cheap_ruler = false)
{
    if (!enus.rows()) {
        return RowVectors(0, 3);
    }
    if (cheap_ruler) {
        return lla2ecef(enu2lla(enus, anchor_lla, cheap_ruler));
    }
    // lossless and and even faster than cheap-ruler
    return apply_transform(T_ecef_enu(anchor_lla).inverse(), enus);
}
inline RowVectors ecef2enu(const Eigen::Ref<const RowVectors> &ecef,
                           std::optional<Eigen::Vector3d> anchor_lla = {},
                           bool cheap_ruler = false)
{
    if (!ecef.rows()) {
        return RowVectors(0, 3);
    }
    if (!anchor_lla) {
        anchor_lla = ecef2lla(ecef(0, 0), ecef(0, 1), ecef(0, 2));
    }
    if (cheap_ruler) {
        return lla2enu(ecef2lla(ecef), anchor_lla, cheap_ruler);
    }
    // lossless and and even faster than cheap-ruler
    return apply_transform(T_ecef_enu(*anchor_lla).inverse(), ecef);
}

} // namespace cubao

#endif
