#pragma once

#include <mapbox/geometry.hpp>

#include <mapbox/variant.hpp>

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

namespace mapbox {
namespace feature {

template <class T>
struct feature
{
    using coordinate_type = T;
    using geometry_type = mapbox::geometry::geometry<T>; // Fully qualified to avoid GCC -fpermissive error.

    geometry_type geometry;
    property_map properties;
    identifier id;
    property_map custom_properties;

    feature()
        : geometry(),
          properties(),
          id() {}
    feature(geometry_type const& geom_)
        : geometry(geom_),
          properties(),
          id() {}
    feature(geometry_type&& geom_)
        : geometry(std::move(geom_)),
          properties(),
          id() {}
    feature(geometry_type const& geom_, property_map const& prop_)
        : geometry(geom_), properties(prop_), id() {}
    feature(geometry_type&& geom_, property_map&& prop_)
        : geometry(std::move(geom_)),
          properties(std::move(prop_)),
          id() {}
    feature(geometry_type const& geom_, property_map const& prop_, identifier const& id_)
        : geometry(geom_),
          properties(prop_),
          id(id_) {}
    feature(geometry_type&& geom_, property_map&& prop_, identifier&& id_)
        : geometry(std::move(geom_)),
          properties(std::move(prop_)),
          id(std::move(id_)) {}
    feature(geometry_type const& geom_, property_map const& prop_, identifier const& id_, property_map const &custom_prop_)
        : geometry(geom_),
          properties(prop_),
          id(id_),
          custom_properties(custom_prop_) {}
    feature(geometry_type&& geom_, property_map&& prop_, identifier&& id_, property_map&& custom_prop_)
        : geometry(std::move(geom_)),
          properties(std::move(prop_)),
          id(std::move(id_)),
          custom_properties(std::move(custom_prop_)) {}
    feature(geometry_type&& geom_, property_map&& prop_, property_map&& custom_prop_)
        : geometry(std::move(geom_)),
          properties(std::move(prop_)),
          custom_properties(std::move(custom_prop_)) {}
};

template <class T>
constexpr bool operator==(feature<T> const& lhs, feature<T> const& rhs)
{
    return lhs.id == rhs.id && lhs.geometry == rhs.geometry && lhs.properties == rhs.properties && lhs.custom_properties == rhs.custom_properties;
}

template <class T>
constexpr bool operator!=(feature<T> const& lhs, feature<T> const& rhs)
{
    return !(lhs == rhs);
}

template <class T, template <typename...> class Cont = std::vector>
struct feature_collection : Cont<feature<T>>
{
    using coordinate_type = T;
    using feature_type = feature<T>;
    using container_type = Cont<feature_type>;
    using size_type = typename container_type::size_type;

    property_map custom_properties;
    feature_collection() = default;
    feature_collection& operator=(feature_collection const &other) = default;
    // template <class... Args>
    // feature_collection(Args&&... args) : container_type(std::forward<Args>(args)...)
    // {
    // }
    feature_collection(std::initializer_list<feature_type> args)
        : container_type(std::move(args)) {}
};

} // namespace feature
} // namespace mapbox
