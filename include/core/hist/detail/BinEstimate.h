// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsAxes.h>
#include <hist/detail/data/WidthControllers.h>
#include <settings/HistogramSettings.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <ranges>

namespace ausaxs::hist::detail {
    namespace bin_estimate {
        constexpr int min_bin_count = 10; // minimum number of bins for all returned histograms
        constexpr int headroom = 2;       // extra bins on top of the geometric bound

        // a point that stores its position as a member, as atoms and waters do
        template<typename T>
        concept PointLike = requires(const T& t) {t.coordinates().x();};

        // a point that is itself a position. atoms forward x()/y()/z() to their coordinates, so they satisfy
        // the requirement too and must be excluded here to keep the two ranges below unambiguous
        template<typename T>
        concept VectorLike = !PointLike<T> && requires(const T& t) {t.x(); t.y(); t.z();};

        template<typename T>
        concept PointRange = std::ranges::input_range<T> && PointLike<std::ranges::range_value_t<T>>;

        template<typename T>
        concept VectorRange = std::ranges::input_range<T> && VectorLike<std::ranges::range_value_t<T>>;

        // a set of coordinates that exposes its positions component-wise
        template<typename T>
        concept CoordinateSet = requires(const T& t) {t.size(); t[0].value.pos;};

        // invoke f(x, y, z) for every point in the set
        template<typename F, CoordinateSet Coords>
        void for_each_point(F& f, const Coords& coords) {
            int size = static_cast<int>(coords.size());
            for (int i = 0; i < size; ++i) {
                const auto& p = coords[i].value.pos;
                f(static_cast<double>(p.x()), static_cast<double>(p.y()), static_cast<double>(p.z()));
            }
        }

        template<typename F, PointRange Range>
        void for_each_point(F& f, const Range& points) {
            for (const auto& point : points) {
                const auto& coordinates = point.coordinates();
                f(static_cast<double>(coordinates.x()), static_cast<double>(coordinates.y()), static_cast<double>(coordinates.z()));
            }
        }

        template<typename F, VectorRange Range>
        void for_each_point(F& f, const Range& points) {
            for (const auto& point : points) {
                f(static_cast<double>(point.x()), static_cast<double>(point.y()), static_cast<double>(point.z()));
            }
        }

        // a container of sets - possibly nested, as in the per-body symmetry data
        template<typename F, std::ranges::input_range Range> requires (!CoordinateSet<Range> && !PointRange<Range> && !VectorRange<Range>)
        void for_each_point(F& f, const Range& sets) {
            for (const auto& set : sets) {for_each_point(f, set);}
        }

        /**
         * @brief A strict upper bound on the maximum distance between any two of the given points.
         *
         * Two independent O(N) bounds are evaluated and the tighter one is returned: the diagonal of the bounding box, which is tight for 
         * elongated structures, and twice the largest distance from the box centre to any point, which is tight for globular ones. 
         */
        template<typename... Sets>
        double max_distance(const Sets&... sets) {
            std::array<double, 3> lo = {
                std::numeric_limits<double>::max(),
                std::numeric_limits<double>::max(),
                std::numeric_limits<double>::max()
            };
            std::array<double, 3> hi = {
                std::numeric_limits<double>::lowest(),
                std::numeric_limits<double>::lowest(),
                std::numeric_limits<double>::lowest()
            };
            auto expand = [&lo, &hi] (double x, double y, double z) {
                const std::array<double, 3> p = {x, y, z};
                for (int k = 0; k < 3; ++k) {
                    lo[k] = std::min(lo[k], p[k]);
                    hi[k] = std::max(hi[k], p[k]);
                }
            };
            (for_each_point(expand, sets), ...);
            if (hi[0] < lo[0]) {return 0;} // no points were seen

            const std::array<double, 3> centre = {(lo[0]+hi[0])/2, (lo[1]+hi[1])/2, (lo[2]+hi[2])/2};
            double r2_max = 0;
            auto radius = [&r2_max, &centre] (double x, double y, double z) {
                double dx = x-centre[0], dy = y-centre[1], dz = z-centre[2];
                r2_max = std::max(r2_max, dx*dx + dy*dy + dz*dz);
            };
            (for_each_point(radius, sets), ...);

            double dx = hi[0]-lo[0], dy = hi[1]-lo[1], dz = hi[2]-lo[2];
            return std::min(std::sqrt(dx*dx + dy*dy + dz*dz), 2*std::sqrt(r2_max));
        }

    }

    /**
     * @brief The number of distance bins required to histogram every pairwise distance within the given sets, as a strict upper bound.
     * @param sets Any number of coordinate sets, or (possibly nested) containers of them.
     */
    template<bool variable_bin_width, typename... Sets>
    int required_bin_count(const Sets&... sets) {
        auto inv_bin_width = static_cast<double>(WidthController<variable_bin_width>::get_inv_width());
        double bins = std::ceil(bin_estimate::max_distance(sets...)*inv_bin_width) + bin_estimate::headroom;
        assert(std::isfinite(bins) && 0 < bins && "Determined bin count is not finite.");
        return std::max<int>(static_cast<int>(bins), bin_estimate::min_bin_count);
    }

}
