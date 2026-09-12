// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <mini/LimitedScan.h>

#include <mini/detail/Parameter.h>  // IWYU pragma: keep
#include <utility/Limit.h>

#include <limits>
#include <list>
#include <numeric>

using namespace ausaxs;
using namespace ausaxs::mini;

LimitedScan::~LimitedScan() = default;

void LimitedScan::set_limit(double limit, bool minima_multiplier) noexcept {
    this->limit = limit;
    this->limit_is_minima_multiplier = minima_multiplier;
}

mini::Landscape LimitedScan::landscape(int evals) {
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
        const auto& bounds_opt = parameters[0].bounds;
        if (!bounds_opt.has_value()) {throw except::bad_order("LimitedScan::landscape: The scanned parameter must be bounded.");}
        const Limit& bounds = *bounds_opt;
        int c = 0;
        std::list<double> last_evals;
        int count = 0;
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
            double avg = std::accumulate(last_evals.begin(), last_evals.end(), 0.0) / static_cast<double>(last_evals.size());

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

    // parameters.size() <= 2 
    throw ausaxs::except::runtime_error("LimitedScan::landscape: Using more than two parameters is currently not implemented.");    
}