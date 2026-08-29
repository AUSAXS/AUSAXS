// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsAxes.h>
#include <hist/detail/data/WidthControllers.h>
#include <settings/Flags.h>
#include <settings/HistogramSettings.h>
#include <utility/Console.h>
#include <utility/Exceptions.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ranges>
#include <string>
#include <string_view>

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
         * @brief A set of coordinates that can be scanned directly, i.e. every CompactCoordinates* specialisation.
         *
         * Anything indexable yielding a .value.pos qualifies. Containers of such sets do not, and are
         * handled by the range overloads below instead - which is what distinguishes the two.
         */
        template<typename T>
        concept CoordinateSet = requires(const T& t) {
            t.size();
            t[0].value.pos;
        };

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
        template<CoordinateSet Coords>
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

        // a container of coordinate sets - possibly nested, as in the per-body symmetry data
        template<std::ranges::input_range Range> requires (!CoordinateSet<Range>)
        void accumulate_bounds(extent_state& s, const Range& sets) {
            for (const auto& set : sets) {accumulate_bounds(s, set);}
        }

        /**
         * @brief The centroid of everything accumulated so far. Only meaningful if s.n != 0.
         */
        inline std::array<double, 3> centroid(const extent_state& s) {
            return {
                s.sum[0]/static_cast<double>(s.n),
                s.sum[1]/static_cast<double>(s.n),
                s.sum[2]/static_cast<double>(s.n)
            };
        }

        // second pass: largest squared distance from the centroid
        template<CoordinateSet Coords>
        void accumulate_radius(double& r2_max, const std::array<double, 3>& centre, const Coords& coords) {
            std::size_t size = coords.size();
            for (std::size_t i = 0; i < size; ++i) {
                const auto& p = coords[i].value.pos;
                double dx = static_cast<double>(p.x()) - centre[0];
                double dy = static_cast<double>(p.y()) - centre[1];
                double dz = static_cast<double>(p.z()) - centre[2];
                r2_max = std::max(r2_max, dx*dx + dy*dy + dz*dz);
            }
        }

        template<std::ranges::input_range Range> requires (!CoordinateSet<Range>)
        void accumulate_radius(double& r2_max, const std::array<double, 3>& centre, const Range& sets) {
            for (const auto& set : sets) {accumulate_radius(r2_max, centre, set);}
        }

        /**
         * @brief Combine the two passes into a strict upper bound on the maximum pairwise distance.
         *
         * @param s        The bounding box and centroid, from accumulate_bounds.
         * @param r2_max   The largest squared distance from that centroid, from accumulate_radius.
         */
        inline double max_distance(const extent_state& s, double r2_max) {
            if (s.n == 0) {return 0;}
            double dx = s.hi[0] - s.lo[0], dy = s.hi[1] - s.lo[1], dz = s.hi[2] - s.lo[2];
            double d_box = std::sqrt(dx*dx + dy*dy + dz*dz);
            double d_sphere = 2*std::sqrt(r2_max);
            return std::min(d_box, d_sphere);
        }

        /**
         * @brief Convert an upper bound on the maximum distance into a bin count.
         */
        inline unsigned int bin_count_from_distance(double d_max, double inv_bin_width) {
            double bins = std::ceil(d_max*inv_bin_width) + headroom;

            // guard against a degenerate bin width or a non-finite coordinate: fall back to the
            // full default range rather than allocating something unusable.
            if (!std::isfinite(bins) || bins <= 0) {
                return std::max<unsigned int>(constants::axes::d_axis.bins, min_bin_count);
            }
            return std::max<unsigned int>(static_cast<unsigned int>(bins), min_bin_count);
        }

        /**
         * @brief As above, picking up the bin width the calculators will actually use.
         */
        template<bool variable_bin_width>
        unsigned int bin_count_from_distance(double d_max) {
            return bin_count_from_distance(d_max, static_cast<double>(WidthController<variable_bin_width>::get_inv_width()));
        }

        /**
         * @brief The bin count for managers that cannot deduce one, i.e. settings::flags::max_bin_count.
         *
         * The partial managers cache histograms across body movements, so a bound deduced from the
         * conformation present at initialization can be exceeded by a later one. They therefore need a
         * fixed upper bound chosen by whoever drives the optimization, rather than deducing their own.
         *
         * @param context The calling manager, for the error message.
         * @throws except::invalid_operation if the flag has not been set.
         */
        inline unsigned int configured_bin_count(std::string_view context) {
            if (settings::flags::max_bin_count == 0) {
                throw except::invalid_operation(
                    std::string(context) + ": settings::flags::max_bin_count has not been set. "
                    "Managers that cache histograms across body movements cannot deduce their own bin count, so the "
                    "caller driving the movements must supply an upper bound - see settings/Flags.h. Subprograms "
                    "performing such movements normally set it to constants::axes::d_axis.bins during initialization."
                );
            }

            static bool warned = false;
            if (!warned && settings::axes::bin_width*settings::flags::max_bin_count < constants::axes::d_axis.max) {
                warned = true;
                console::print_warning(
                    std::string(context) + ": the current bin width (" + std::to_string(settings::axes::bin_width) + "Å) and "
                    "bin count (" + std::to_string(settings::flags::max_bin_count) + ") only cover distances up to " +
                    std::to_string(int(settings::axes::bin_width*settings::flags::max_bin_count)) + "Å, which is less than the "
                    "recommended " + std::to_string(int(constants::axes::d_axis.max)) + "Å. Larger structures will be truncated."
                );
            }
            return std::max(settings::flags::max_bin_count, min_bin_count);
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
     * @param coords Any number of coordinate sets, or (possibly nested) containers of them.
     */
    template<typename... Coords>
    double estimate_max_distance(const Coords&... coords) {
        bin_estimate::extent_state s;
        (bin_estimate::accumulate_bounds(s, coords), ...);
        if (s.n == 0) {return 0;}

        auto centre = bin_estimate::centroid(s);
        double r2_max = 0;
        (bin_estimate::accumulate_radius(r2_max, centre, coords), ...);
        return bin_estimate::max_distance(s, r2_max);
    }

    /**
     * @brief Compute the number of distance bins required to histogram every pairwise distance
     *        within the given coordinate sets, as a strict upper bound.
     *
     * Because it is an upper bound, every bin above the true maximum distance is zero, so the
     * trim-to-last-non-zero step the managers already perform yields an identical histogram to
     * one computed over a larger allocation.
     *
     * @param inv_bin_width The reciprocal bin width the calculators will use.
     */
    template<typename... Coords>
    unsigned int estimate_bin_count(double inv_bin_width, const Coords&... coords) {
        return bin_estimate::bin_count_from_distance(estimate_max_distance(coords...), inv_bin_width);
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
