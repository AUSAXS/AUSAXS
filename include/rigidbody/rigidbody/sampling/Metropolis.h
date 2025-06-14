#pragma once

#include <cmath>
#include <random>
#include <functional>

namespace ausaxs::math::sampling {
    class Metropolis {
        Metropolis() {
            gen = std::mt19937(rd());
        }

        private:
            std::function<double(double)> pdf;
            inline static std::random_device rd;
            inline static std::mt19937 gen;

            double random_uniform() {
                static std::uniform_real_distribution<double> dist(0, 1);
                return dist(gen);
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