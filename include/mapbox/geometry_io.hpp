#pragma once

#include <mapbox/geometry/empty.hpp>
#include <mapbox/feature.hpp>

#include <iostream>
#include <string>

namespace mapbox {
namespace geometry {

std::ostream& operator<<(std::ostream& os, const empty&)
{
    return os << "[]";
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

template <typename T>
std::ostream& operator<<(std::ostream& os, const point<T>& point)
{
    os << "[" << point.x << "," << point.y;
    if (point.z != 0) {
        os << "," << point.z;
    }
    os << "]";
    return os;
}

#pragma GCC diagnostic pop

template <typename T, template <class, class...> class C, class... Args>
std::ostream& operator<<(std::ostream& os, const C<T, Args...>& cont)
{
    os << "[";
    for (auto it = cont.cbegin();;)
    {
        os << *it;
        if (++it == cont.cend())
        {
            break;
        }
        os << ",";
    }
    return os << "]";
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const line_string<T>& geom)
{
    return os << static_cast<typename line_string<T>::container_type>(geom);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const linear_ring<T>& geom)
{
    return os << static_cast<typename linear_ring<T>::container_type>(geom);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const polygon<T>& geom)
{
    return os << static_cast<typename polygon<T>::container_type>(geom);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const multi_point<T>& geom)
{
    return os << static_cast<typename multi_point<T>::container_type>(geom);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const multi_line_string<T>& geom)
{
    return os << static_cast<typename multi_line_string<T>::container_type>(geom);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const multi_polygon<T>& geom)
{
    return os << static_cast<typename multi_polygon<T>::container_type>(geom);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const geometry<T>& geom)
{
    geometry<T>::visit(geom, [&](const auto& g) { os << g; });
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const geometry_collection<T>& geom)
{
    return os << static_cast<typename geometry_collection<T>::container_type>(geom);
}

} // namespace geometry

namespace feature {

std::ostream& operator<<(std::ostream& os, const null_value_t&)
{
    return os << "[]";
}

} // namespace feature
} // namespace mapbox
