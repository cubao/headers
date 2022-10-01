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

// set MAPBOX_GEOMETRY_ENABLE_CUSTOM_PROPERTIES to 1, then you can load
//      {
//            "type": "Feature",
//            "geometry": { "type": "Point", "coordinates": [1, 2] },
//            "properties": {},
//            "key": "value",
//            "list": [],
//            ...
//          }
//      }
// with custom_properties := {
//      "key": "value",
//      "list": [],
//      ...
// }
// just like in geobuf: https://github.com/mapbox/geobuf/blob/master/geobuf.proto#L25
#if MAPBOX_GEOMETRY_ENABLE_CUSTOM_PROPERTIES
    property_map custom_properties;
#endif

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
#if MAPBOX_GEOMETRY_ENABLE_CUSTOM_PROPERTIES
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
#endif
};

template <class T>
constexpr bool operator==(feature<T> const& lhs, feature<T> const& rhs)
{
#if MAPBOX_GEOMETRY_ENABLE_CUSTOM_PROPERTIES
    return lhs.id == rhs.id && lhs.geometry == rhs.geometry && lhs.properties == rhs.properties && lhs.custom_properties == rhs.custom_properties;
#else
    return lhs.id == rhs.id && lhs.geometry == rhs.geometry && lhs.properties == rhs.properties;
#endif
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

#if MAPBOX_GEOMETRY_ENABLE_CUSTOM_PROPERTIES
    property_map custom_properties;
    feature_collection() = default;
    feature_collection& operator=(feature_collection const &other) = default;
    feature_collection(std::initializer_list<feature_type> args)
        : container_type(std::move(args)) {}
#else
    template <class... Args>
    feature_collection(Args&&... args) : container_type(std::forward<Args>(args)...)
    {
    }
    feature_collection(std::initializer_list<feature_type> args)
        : container_type(std::move(args)) {}
#endif
};

} // namespace feature
} // namespace mapbox
