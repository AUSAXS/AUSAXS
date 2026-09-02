// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsAxes.h>
#include <hist/detail/data/WidthControllers.h>
#include <settings/Flags.h>
#include <settings/HistogramSettings.h>
#include <utility/Console.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ranges>
#include <string>

namespace ausaxs::hist::detail {
    namespace bin_estimate {
        constexpr unsigned int min_bin_count = 10; // minimum number of bins for all returned histograms
        constexpr unsigned int headroom = 2;       // extra bins on top of the geometric bound

        // a set of coordinates that exposes its positions component-wise
        template<typename T>
        concept CoordinateSet = requires(const T& t) {t.size(); t.x(0u); t.y(0u); t.z(0u);};

        // invoke f(x, y, z) for every point in the set
        template<typename F, CoordinateSet Coords>
        void for_each_point(F& f, const Coords& coords) {
            unsigned int size = static_cast<unsigned int>(coords.size());
            for (unsigned int i = 0; i < size; ++i) {
                f(static_cast<double>(coords.x(i)), static_cast<double>(coords.y(i)), static_cast<double>(coords.z(i)));
            }
        }

        // a container of sets - possibly nested, as in the per-body symmetry data
        template<typename F, std::ranges::input_range Range> requires (!CoordinateSet<Range>)
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
            auto expand = [&lo, &hi] (double x, double y, double z) {
                const double p[3] = {x, y, z};
                for (int k = 0; k < 3; ++k) {
                    lo[k] = std::min(lo[k], p[k]);
                    hi[k] = std::max(hi[k], p[k]);
                }
            };
            (for_each_point(expand, sets), ...);
            if (hi[0] < lo[0]) {return 0;} // no points were seen

            const double centre[3] = {(lo[0]+hi[0])/2, (lo[1]+hi[1])/2, (lo[2]+hi[2])/2};
            double r2_max = 0;
            auto radius = [&r2_max, &centre] (double x, double y, double z) {
                double dx = x-centre[0], dy = y-centre[1], dz = z-centre[2];
                r2_max = std::max(r2_max, dx*dx + dy*dy + dz*dz);
            };
            (for_each_point(radius, sets), ...);

            double dx = hi[0]-lo[0], dy = hi[1]-lo[1], dz = hi[2]-lo[2];
            return std::min(std::sqrt(dx*dx + dy*dy + dz*dz), 2*std::sqrt(r2_max));
        }

        /**
         * @brief The bin count for managers that cannot deduce one.
         * @return settings::flag::max_bin_count
         */
        inline unsigned int configured_bin_count() {
            assert(settings::flags::max_bin_count != 0 && "max_bin_count has not been set.");

            static bool warned = false;
            if (!warned && settings::axes::bin_width*settings::flags::max_bin_count < constants::axes::d_axis.max) {
                warned = true;
                console::print_warning(
                    "The current bin width (" + std::to_string(settings::axes::bin_width) + "Å) and bin count (" + 
                    std::to_string(settings::flags::max_bin_count) + ") only cover distances up to " + 
                    std::to_string(int(settings::axes::bin_width*settings::flags::max_bin_count)) + "Å, which is less than the recommended " + 
                    std::to_string(int(constants::axes::d_axis.max)) + "Å. Larger structures will cause segfaults."
                );
            }
            return std::max(settings::flags::max_bin_count, min_bin_count);
        }
    }

    /**
     * @brief The number of distance bins required to histogram every pairwise distance within the given sets, as a strict upper bound.
     * @param sets Any number of coordinate sets, or (possibly nested) containers of them.
     */
    template<bool variable_bin_width, typename... Sets>
    unsigned int required_bin_count(const Sets&... sets) {
        double inv_bin_width = static_cast<double>(WidthController<variable_bin_width>::get_inv_width());
        double bins = std::ceil(bin_estimate::max_distance(sets...)*inv_bin_width) + bin_estimate::headroom;
        assert(std::isfinite(bins) && 0 < bins && "Determined bin count is not finite.");
        return std::max<unsigned int>(static_cast<unsigned int>(bins), bin_estimate::min_bin_count);
    }
}
