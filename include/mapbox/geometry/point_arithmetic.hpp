#pragma once

namespace mapbox {
namespace geometry {

template <typename T>
point<T> operator+(point<T> const& lhs, point<T> const& rhs)
{
    return point<T>(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
}

template <typename T>
point<T> operator+(point<T> const& lhs, T const& rhs)
{
    return point<T>(lhs.x + rhs, lhs.y + rhs, lhs.z + rhs);
}

template <typename T>
point<T> operator-(point<T> const& lhs, point<T> const& rhs)
{
    return point<T>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

template <typename T>
point<T> operator-(point<T> const& lhs, T const& rhs)
{
    return point<T>(lhs.x - rhs, lhs.y - rhs, lhs.z - rhs);
}

template <typename T>
point<T> operator*(point<T> const& lhs, point<T> const& rhs)
{
    return point<T>(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z);
}

template <typename T>
point<T> operator*(point<T> const& lhs, T const& rhs)
{
    return point<T>(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs);
}

template <typename T>
point<T> operator/(point<T> const& lhs, point<T> const& rhs)
{
    return point<T>(lhs.x / rhs.x, lhs.y / rhs.y, lhs.z / rhs.z);
}

template <typename T>
point<T> operator/(point<T> const& lhs, T const& rhs)
{
    return point<T>(lhs.x / rhs, lhs.y / rhs, lhs.z / rhs);
}

template <typename T>
point<T>& operator+=(point<T>& lhs, point<T> const& rhs)
{
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    lhs.z += rhs.z;
    return lhs;
}

template <typename T>
point<T>& operator+=(point<T>& lhs, T const& rhs)
{
    lhs.x += rhs;
    lhs.y += rhs;
    lhs.z += rhs;
    return lhs;
}

template <typename T>
point<T>& operator-=(point<T>& lhs, point<T> const& rhs)
{
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    lhs.z -= rhs.z;
    return lhs;
}

template <typename T>
point<T>& operator-=(point<T>& lhs, T const& rhs)
{
    lhs.x -= rhs;
    lhs.y -= rhs;
    lhs.z -= rhs;
    return lhs;
}

template <typename T>
point<T>& operator*=(point<T>& lhs, point<T> const& rhs)
{
    lhs.x *= rhs.x;
    lhs.y *= rhs.y;
    lhs.z *= rhs.z;
    return lhs;
}

template <typename T>
point<T>& operator*=(point<T>& lhs, T const& rhs)
{
    lhs.x *= rhs;
    lhs.y *= rhs;
    lhs.z *= rhs;
    return lhs;
}

template <typename T>
point<T>& operator/=(point<T>& lhs, point<T> const& rhs)
{
    lhs.x /= rhs.x;
    lhs.y /= rhs.y;
    lhs.z /= rhs.z;
    return lhs;
}

template <typename T>
point<T>& operator/=(point<T>& lhs, T const& rhs)
{
    lhs.x /= rhs;
    lhs.y /= rhs;
    lhs.z /= rhs;
    return lhs;
}

} // namespace geometry
} // namespace mapbox
