// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsAxes.h>
#include <hist/detail/data/WidthControllers.h>
#include <settings/HistogramSettings.h>
#include <utility/Console.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>

namespace ausaxs::hist::detail {
    namespace bin_estimate {
        /**
         * @brief The smallest histogram the managers may allocate.
         *
         * The trim loops in the managers start at index size()-1 and stop at index 10, so
         * anything smaller than 11 bins either skips the loop entirely or - for an empty
         * histogram - underflows the unsigned index.
         */
        constexpr unsigned int min_bin_count = 11;

        /**
         * @brief Extra bins added on top of the geometric bound.
         *
         * The calculators bin as std::round(inv_width*d), which can carry the index up to
         * half a bin past the bound, and they evaluate d in single precision. Two bins of
         * headroom covers both with room to spare.
         */
        constexpr unsigned int headroom = 2;

        /**
         * @brief Warn - once - that the configured bin count was too small for the structure being calculated.
         */
        inline void warn_bin_count_raised(unsigned int required) {
            static bool warned = false;
            if (warned) {return;}
            warned = true;
            console::print_warning(
                "hist::detail::estimate_bin_count: the configured bin count (" + std::to_string(settings::axes::bin_count) + ") "
                "is too small for this structure, which requires at least " + std::to_string(required) + " bins. "
                "The histogram allocation has been raised accordingly; the configured value would have caused an out-of-bounds write."
            );
        }

        struct extent_state {
            double lo[3] = {
                std::numeric_limits<double>::max(),
                std::numeric_limits<double>::max(),
                std::numeric_limits<double>::max()
            };
            double hi[3] = {
                std::numeric_limits<double>::lowest(),
                std::numeric_limits<double>::lowest(),
                std::numeric_limits<double>::lowest()
            };
            double sum[3] = {0, 0, 0};
            std::size_t n = 0;
        };

        // first pass: axis-aligned bounding box and centroid
        template<typename Coords>
        void accumulate_bounds(extent_state& s, const Coords& coords) {
            std::size_t size = coords.size();
            for (std::size_t i = 0; i < size; ++i) {
                const auto& p = coords[i].value.pos;
                const double xyz[3] = {
                    static_cast<double>(p.x()), static_cast<double>(p.y()), static_cast<double>(p.z())
                };
                for (int k = 0; k < 3; ++k) {
                    s.lo[k] = std::min(s.lo[k], xyz[k]);
                    s.hi[k] = std::max(s.hi[k], xyz[k]);
                    s.sum[k] += xyz[k];
                }
            }
            s.n += size;
        }

        // second pass: largest squared distance from the centroid
        template<typename Coords>
        void accumulate_radius(double& r2_max, const double (&centre)[3], const Coords& coords) {
            std::size_t size = coords.size();
            for (std::size_t i = 0; i < size; ++i) {
                const auto& p = coords[i].value.pos;
                double dx = static_cast<double>(p.x()) - centre[0];
                double dy = static_cast<double>(p.y()) - centre[1];
                double dz = static_cast<double>(p.z()) - centre[2];
                r2_max = std::max(r2_max, dx*dx + dy*dy + dz*dz);
            }
        }
    }

    /**
     * @brief Compute a strict upper bound on the maximum distance between any two points in the given coordinate sets.
     *
     * Two independent O(N) bounds are evaluated and the tighter one is returned:
     *  - the diagonal of the axis-aligned bounding box enclosing every point, and
     *  - twice the largest distance from the centroid to any point.
     * The first is tight for elongated structures, the second for globular ones.
     *
     * @param coords Any number of coordinate sets. Any type with size() and an operator[] returning
     *               something with a .value.pos member works, i.e. every CompactCoordinates* specialisation.
     */
    template<typename... Coords>
    double estimate_max_distance(const Coords&... coords) {
        bin_estimate::extent_state s;
        (bin_estimate::accumulate_bounds(s, coords), ...);
        if (s.n == 0) {return 0;}

        double dx = s.hi[0] - s.lo[0], dy = s.hi[1] - s.lo[1], dz = s.hi[2] - s.lo[2];
        double d_box = std::sqrt(dx*dx + dy*dy + dz*dz);

        double centre[3] = {
            s.sum[0]/static_cast<double>(s.n),
            s.sum[1]/static_cast<double>(s.n),
            s.sum[2]/static_cast<double>(s.n)
        };
        double r2_max = 0;
        (bin_estimate::accumulate_radius(r2_max, centre, coords), ...);
        double d_sphere = 2*std::sqrt(r2_max);

        return std::min(d_box, d_sphere);
    }

    /**
     * @brief Compute the number of distance bins required to histogram every pairwise distance
     *        within the given coordinate sets, as a strict upper bound.
     *
     * This is intended as a drop-in replacement for settings::axes::bin_count at the allocation
     * sites of the histogram managers. Because it is an upper bound, every bin above the true
     * maximum distance is zero - exactly as it is when allocating the full default range - so the
     * trim-to-last-non-zero step the managers already perform yields a bit-identical histogram.
     *
     * @param inv_bin_width The reciprocal bin width the calculators will use.
     */
    template<typename... Coords>
    unsigned int estimate_bin_count(double inv_bin_width, const Coords&... coords) {
        double d_max = estimate_max_distance(coords...);
        double bins = std::ceil(d_max*inv_bin_width) + bin_estimate::headroom;

        // guard against a degenerate bin width or a non-finite coordinate: fall back to the
        // configured range rather than allocating something unusable.
        if (!std::isfinite(bins) || bins <= 0) {
            return std::max<unsigned int>(settings::axes::bin_count, bin_estimate::min_bin_count);
        }

        // Note that the estimate deliberately overrides a *too small* configured bin_count.
        // Allocating fewer bins than the structure needs is not a slower calculation, it is an
        // out-of-bounds write - see the warning in settings::axes::bin_count. Since the estimate
        // is an upper bound, honouring it costs nothing and removes that failure mode entirely.
        if (bins > static_cast<double>(settings::axes::bin_count)) {
            bin_estimate::warn_bin_count_raised(static_cast<unsigned int>(bins));
        }
        return std::max<unsigned int>(static_cast<unsigned int>(bins), bin_estimate::min_bin_count);
    }

    /**
     * @brief Convenience overload picking up the bin width the calculators will actually use.
     *
     * @tparam variable_bin_width The same flag the manager and its coordinate sets are templated on.
     */
    template<bool variable_bin_width, typename... Coords>
    unsigned int required_bin_count(const Coords&... coords) {
        return estimate_bin_count(static_cast<double>(WidthController<variable_bin_width>::get_inv_width()), coords...);
    }
}
