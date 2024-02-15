#pragma once

#include <mini/Scan.h>

#include <limits>

namespace mini {
    /**
     * @brief A scan that may terminate early if a condition is met. The scan will start at the maximum x value and work its way to the minimum.
     * 
     * The scan will terminate early (but never before half of the iterations have been performed) if one of the following conditions have been met:
     *      1. The function value exceeds the limit. The limit can either be a multiplier of the minima or a fixed value.
     *      2. The function value is on an uphill slope. This is determined by calculating the average of the last 7 function values. 
     *          If the current value is greater than the average, the counter is incremented. If the counter reaches 3, the scan is terminated.
     */
    class LimitedScan : public Scan {
        public:
            using Scan::Scan;

            /**
             * @brief Destructor.
             */
            ~LimitedScan() override;

            /**
             * @brief Set the maximum function value before terminating.
             */
            void set_limit(double limit, bool minima_multiplier = false) noexcept;

            /**
             * @brief Generate a landscape of the function.
             *        The scan will start at the maximum x value and work its way to the minimum.
             *        This will terminate early if the function value exceeds the limit.
             */
            mini::Landscape landscape(unsigned int evals) override;

        private: 
            double limit = std::numeric_limits<double>::max();
            bool limit_is_minima_multiplier = false;
    };
}
