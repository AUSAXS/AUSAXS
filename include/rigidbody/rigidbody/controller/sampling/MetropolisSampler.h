#pragma once

#include <rigidbody/detail/Configuration.h>
#include <utility/Random.h>

#include <cmath>
#include <random>
#include <functional>

namespace ausaxs::rigidbody::sampling {
    class MetropolisSampler {
        MetropolisSampler() {}

        private:
            std::function<double(double)> chi2;
            std::function<void(detail::Configuration& best, double new_chi2)> accept_update;
            std::function<void(const detail::Configuration&)> reject_update;

            double pdf(double x) {
                return std::exp(-chi2(x));
            }

            double random_uniform() {
                static std::uniform_real_distribution<double> dist(0, 1);
                return dist(random::generator());
            }

            double P_ratio(double x_old, double x_new) {
                return pdf(x_new)/pdf(x_old)*std::exp(-x_new + x_old);
            }

            bool accept(double x_old, double x_new) {
                double p = P_ratio(x_old, x_new);
                return random_uniform() < p;
            }
    };
}