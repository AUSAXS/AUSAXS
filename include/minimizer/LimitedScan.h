#pragma once

#include <minimizer/Scan.h>
#include <limits>

namespace mini {
    class LimitedScan : public Scan {
        public:
            using Scan::Scan;

            /**
             * @brief Set the maximum function value before terminating.
             */
            void set_limit(double limit) noexcept {
                this->limit = limit;
            }

            /**
             * @brief Generate a landscape of the function.
             *        The scan will start at the maximum x value and work its way to the minimum.
             *        This will terminate early if the function value exceeds the limit.
             */
            Dataset2D landscape(unsigned int evals) override {
                // check if the minimizer has already been called
                if (!evaluations.empty()) {
                    // if so, we can just reuse its result
                    Dataset2D data;
                    std::for_each(evaluations.begin(), evaluations.end(), [&data] (const Evaluation& eval) {data.push_back(eval.vals[0], eval.fval);});
                    return data;
                }

                if (parameters.size() == 1) {
                    const Limit& bounds = parameters[0].bounds.value();
                    unsigned int c = 0;
                    for (double val = bounds.max; bounds.min < val; val -= bounds.span()/evals) {
                        double fval = function(&val);

                        // if we are more than half-way through the scan and the function value is greater than the limit, stop
                        if (evals/2 < c++ && limit < fval) {break;}
                    }
                    return get_evaluated_points();
                } 
                
                else if (parameters.size() == 2) {
                    throw except::unexpected("Error in LimitedScan::landscape: Not implemented.");
                } 
                
                else {
                    throw except::unexpected("Error in LimitedScan::landscape: Not implemented.");
                }                
            }

        private: 
            double limit = std::numeric_limits<double>::max();
    };
}
