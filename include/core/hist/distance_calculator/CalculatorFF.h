// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactorType.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/Calculator.h>
#include <utility/observer_ptr.h>

#include <cassert>
#include <vector>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Queues the form factor-resolved pairwise distance histograms on a Calculator.
     */
    template<bool weighted_bins, bool variable_bin_width>
    class CalculatorFF {
        using CompactCoordinates_t = hist::detail::CompactCoordinates<variable_bin_width>;
        using PartitionedCoordinates_t = std::vector<CompactCoordinates_t>;
        public:
            /**
             * @brief Queue onto @a calculator, which must outlive this helper.
             */
            explicit CalculatorFF(Calculator<weighted_bins, variable_bin_width>& calculator)
                : calculator(&calculator), n_ff(form_factor::get_active_count()) 
            {}

            /**
             * @brief Queue a self-correlation calculation.
             *        The data reference must be valid until run() is called. 
             * 
             * @param a The data to calculate the self-correlation for, partitioned by form factor.
             * @param first The first handle to accumulate into. The n_ff rows from first are used. 
             */
            void enqueue_self_by_ff(const PartitionedCoordinates_t& a, int first) {
                assert(static_cast<int>(a.size()) == n_ff && "CalculatorFF::enqueue_self_by_ff: expected one set per active form factor.");
                for (int ff1 = 0; ff1 < n_ff; ++ff1) {
                    calculator->enqueue_calculate_self(a[ff1], first + ff1*n_ff + ff1);
                    for (int ff2 = ff1+1; ff2 < n_ff; ++ff2) {
                        calculator->enqueue_calculate_cross(a[ff1], a[ff2], first + ff1*n_ff + ff2, 2);
                    }
                }
            }

            /**
             * @brief Queue a cross-correlation calculation.
             *        The data references must be valid until run() is called. 
             * 
             * @param a The first set of data to calculate the cross-correlation for, partitioned by form factor.
             * @param b The second set of data to calculate the cross-correlation for, partitioned by form factor.
             * @param first The first handle to accumulate into. The n_ff rows from first are used. 
             */
            void enqueue_cross_by_ff(const PartitionedCoordinates_t& a, const CompactCoordinates_t& b, int first) {
                assert(static_cast<int>(a.size()) == n_ff && "CalculatorFF::enqueue_cross_by_ff: expected one set per active form factor.");
                for (int ff1 = 0; ff1 < n_ff; ++ff1) {
                    calculator->enqueue_calculate_cross(a[ff1], b, first + ff1, 1);
                }
            }

            /**
             * @brief Queue the self-correlation of @a b, ignoring form factors, into row @a h.
             */
            void enqueue_self_flat(const CompactCoordinates_t& b, int h) {
                calculator->enqueue_calculate_self(b, h);
            }

        private:
            observer_ptr<Calculator<weighted_bins, variable_bin_width>> calculator;
            int n_ff;
    };
}
