// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cassert>
#include <cmath>
#include <vector>

namespace ausaxs {
    struct MovingAverage {
        public: 
            /**
             * @brief Calculate the moving average of the vector.
             *        For each point in the vector, the average of the @a window_size points in the current window is calculated.
             *        Edge points are calculated with smaller windows. 
             * 
             * @param data The data to average. Container must support indexing and size().
             * @param window_size The window size. 
             */
            template<typename T>
            static std::vector<double> average(const T& data, int window_size) {
                validate_input(data.size(), window_size);
                return weighted_average(data, std::vector<double>(window_size, 1));
            }

            /**
             * @brief Calculate the moving average of the vector.
             *        For each point in the vector, the weighted average of the @a window_size points in the current window is calculated.
             *        The weight is defined as (1/2)^i, where i is the index distance from the current point.
             *        Edge points are calculated with smaller windows. 
             * 
             * @param data The data to average. Container must support indexing [] and size().
             * @param window_size The window size. 
             */
            template<typename T>
            static std::vector<double> average_half(const T& data, int window_size) {
                validate_input(data.size(), window_size);
                std::vector<double> weights(window_size, 1);

                // define weights
                int steps = (window_size-1) / 2;
                int mid = steps;
                for (int i = 1; i < steps+1; i++) {
                    double val = 1/std::pow(1.5, i);
                    weights[mid-i] = val;
                    weights[mid+i] = val;
                }

                return weighted_average(data, std::move(weights));
            }

        private: 
            static void validate_input(int N, int window_size);

            template<typename T>
            static std::vector<double> weighted_average(const T& data, std::vector<double> weights) {
                std::size_t window_size = weights.size();
                std::vector<double> averages(data.size());
                int mid = static_cast<int>((weights.size()-1)/2);

                auto average = [&data, &weights, &mid] (int index, int steps) {
                    double sum = 0;
                    double w_sum = 0;
                    for (int j = -steps; j < steps+1; j++) {
                        sum += data[index+j] * weights[mid+j];
                        w_sum += weights[mid+j];
                    }
                    assert(w_sum != 0 && "MovingAverage::average: Division by zero.");
                    return sum/w_sum;
                };

                int steps = static_cast<int>((window_size-1)/2);
                for (int i = steps; i < static_cast<int>(data.size()) - steps; i++) {
                    averages[i] = average(i, steps);
                }

                for (int i = 0; i < steps; i++) {
                    averages[i] = average(i, i);
                    averages[data.size()-1-i] = average(data.size()-1-i, i);
                }

                return averages;
            }
    };
}