#pragma once

#include <mini/Scan.h>

#include <limits>
#include <list>

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
            ~LimitedScan() override = default;

            /**
             * @brief Set the maximum function value before terminating.
             */
            void set_limit(double limit, bool minima_multiplier = false) noexcept {
                this->limit = limit;
                this->limit_is_minima_multiplier = minima_multiplier;
            }

            /**
             * @brief Generate a landscape of the function.
             *        The scan will start at the maximum x value and work its way to the minimum.
             *        This will terminate early if the function value exceeds the limit.
             */
            mini::Landscape landscape(unsigned int evals) override {
                // check if the minimizer has already been called
                if (!evaluations.evals.empty()) {
                    return evaluations; // if so, we can just reuse its result
                }

                double current_min = std::numeric_limits<double>::max();
                std::function<bool(double)> stop_condition;
                if (limit_is_minima_multiplier) {
                    stop_condition = [&](double val) {
                        return current_min*limit < val;
                    };
                } else {
                    stop_condition = [&](double val) {
                        return limit < val;
                    };
                }

                if (parameters.size() == 1) {
                    const Limit& bounds = parameters[0].bounds.value();
                    unsigned int c = 0;
                    std::list<double> last_evals;
                    unsigned int count = 0;
                    for (double val = bounds.max; bounds.min < val; val -= bounds.span()/evals) {
                        double fval = function({val});
                        current_min = std::min(current_min, fval);

                        // add the evaluation to the list
                        if (last_evals.size() < 7) {
                            last_evals.push_front(fval);
                        } else {
                            last_evals.pop_back();
                            last_evals.push_front(fval);
                        }
                        
                        // calculate average of list
                        double avg = std::accumulate(last_evals.begin(), last_evals.end(), 0.0) / last_evals.size();

                        // if we are more than half-way through the scan, we check for the stop condition
                        if (evals*0.7 < ++c) {
                            // if the fval is greater than the limit and we are not converging on a solution, we stop
                            if (stop_condition(avg) && stop_condition(fval)) {
                                ++count;
                                // three consecutive fvals greater than the avg also means we stop
                                if (3 == count) {
                                    break;
                                }
                            } else {
                                count = 0;
                            }
                        }
                    }
                    return get_evaluated_points();
                } 
                
                else { // parameters.size() <= 2 
                    throw except::unexpected("LimitedScan::landscape: Using more than two parameters is currently not implemented.");
                } 
            }

        private: 
            double limit = std::numeric_limits<double>::max();
            bool limit_is_minima_multiplier = false;
    };
}
