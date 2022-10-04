#pragma once

#include <utility/Exceptions.h>

#include <vector>
#include <math.h>

struct MovingAverage {
    public: 
        /**
         * @brief Calculate the moving average of the vector.
         *        For each point in the vector, the average of the @a window_size points in the current window is calculated.
         *        The first and last @a (window_size-1)/2 points are skipped, so the averaged vector is @a window_size-1 smaller than the input.
         * 
         * @param data The data to average. Container must support indexing and size().
         * @param window_size The window size. 
         */
        template<typename T>
        static std::vector<double> average(const T& data, unsigned int window_size) {
            validate_input(window_size);
            return weighted_average(data, std::vector<double>(window_size, 1));
        }

        /**
         * @brief Calculate the moving average of the vector.
         *        For each point in the vector, the weighted average of the @a window_size points in the current window is calculated.
         *        The weight is defined as (1/2)^i, where i is the index distance from the current point.
         *        The first and last @a (window_size-1)/2 points are skipped, so the averaged vector is @a window_size-1 smaller than the input.
         * 
         * @param data The data to average. Container must support indexing [] and size().
         * @param window_size The window size. 
         */
        template<typename T>
        static std::vector<double> average_half(const T& data, unsigned int window_size) {
            validate_input(window_size);
            std::vector<double> weights(window_size, 1);

            // define weights
            unsigned int steps = (window_size-1) / 2;
            unsigned int mid = steps;
            for (unsigned int i = 1; i < steps+1; i++) {
                double val = 1/std::pow(2, i);
                weights[mid-i] = val;
                weights[mid+i] = val;
            }

            return weighted_average(data, weights);
        }

    private: 
        static void validate_input(unsigned int window_size) {
            if (window_size % 2 == 0) {
                throw except::invalid_operation("Error in MovingAverager::validate_input: Window_size must be odd");
            }
        }

        template<typename T>
        static std::vector<double> weighted_average(const T& data, std::vector<double> weights) {
            unsigned int window_size = weights.size();
            std::vector<double> averages(data.size() - window_size + 1);

            unsigned int offset = (window_size-1)/2;
            for (unsigned int i = offset; i < data.size() - offset; i++) {
                double sum = 0;
                for (unsigned int j = i - offset; j < i + offset+1; j++) {
                    sum += data[j] * weights[j - (i - offset)];
                }
                averages[i - offset] = sum/window_size;
            }

            return averages;
        }
};