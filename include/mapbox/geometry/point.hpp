#pragma once

namespace mapbox {
namespace geometry {

template <typename T>
struct point
{
    using coordinate_type = T;

    constexpr point()
        : x(), y(), z()
    {
    }
    constexpr point(T x_, T y_, T z_ = 0)
        : x(x_), y(y_), z(z_)
    {
    }

    T x;
    T y;
    T z;
};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

template <typename T>
constexpr bool operator==(point<T> const& lhs, point<T> const& rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

#pragma GCC diagnostic pop

template <typename T>
constexpr bool operator!=(point<T> const& lhs, point<T> const& rhs)
{
    return !(lhs == rhs);
}

} // namespace geometry
} // namespace mapbox
