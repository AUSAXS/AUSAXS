#pragma once

#include <minimizer/Scan.h>
#include <limits>
#include <list>

namespace mini {
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
                    std::list<double> last_evals;
                    unsigned int count = 0;
                    for (double val = bounds.max; bounds.min < val; val -= bounds.span()/evals) {
                        double fval = function(&val);

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
                        if (evals*0.8 < c++) {
                            std::cout << limit << " " << fval << " " << avg << std::endl;
                            // if the fval is greater than the limit and we are not converging on a solution, we stop
                            if (limit < fval && avg < fval) {
                                count++;
                                // three consecutive fvals greater than the avg means we stop
                                if (3 < count) {
                                    break;
                                }
                            } else {
                                count = 0;
                            }
                        }
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
