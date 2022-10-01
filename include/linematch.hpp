#pragma once

#include <Eigen/Core>
#include <tl/optional.hpp>

namespace linematch
{
// https://github.com/anvaka/isect/blob/80832e75bf8f197845e52ea52c6ca72935abb24a/build/isect.js#L869
// (buggy?)
// https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect/1968345#1968345
// https://coliru.stacked-crooked.com/a/624e6e0eabc8a103
// returns [x, y, t, s]
// Note that: won't handle seg1 == seg2
inline tl::optional<Eigen::Vector4d>
intersect_segments(const Eigen::Vector2d &a1, const Eigen::Vector2d &a2,
                   const Eigen::Vector2d &b1, const Eigen::Vector2d &b2)
{
    double p0_x = a1[0], p0_y = a1[1];
    double p2_x = b1[0], p2_y = b1[1];
    double s1_x = a2[0] - a1[0];
    double s1_y = a2[1] - a1[1];
    double s2_x = b2[0] - b1[0];
    double s2_y = b2[1] - b1[1];
    double div = s1_x * s2_y - s2_x * s1_y;
    if (div == 0.0) {
        return {};
    }
    double inv = 1.0 / div;
    double s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) * inv;
    if (s < 0.0 || s > 1.0) {
        return {};
    }
    double t = (s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) * inv;
    if (t < 0.0 || t > 1.0) {
        return {};
    }
    return Eigen::Vector4d(p0_x + (t * s1_x), p0_y + (t * s1_y), t, s);
}
} // namespace linematch