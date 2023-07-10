#pragma once

#include <math/slices/Slice.h>
#include <utility/Concepts.h>

#include <vector>
#include <algorithm>
#include <numeric>
#include <type_traits>
#include <initializer_list>
#include <math.h>
#include <iostream>
#include <concepts>

namespace stats {
    /**
     * @brief Calculate the weighted mean of a container class.
     */
    template<container_type T, container_type Q>
    double weighted_mean(const T& x, const Q& xerr) noexcept {
        double sum_wx = 0;
        double sum_w = 0;
        for (unsigned int i = 0; i < x.size(); i++) {
            double w = 1.0/(xerr[i]*xerr[i]);
            sum_wx += w*x[i];
            sum_w += w;
        }
        return sum_wx/sum_w;
    }

    /**
     * @brief Calculate the error in the weighted mean of a container class.
     */
    template<container_type T>
    double weighted_mean_error(const T& xerr) noexcept {
        double sum = 0;
        for (unsigned int i = 0; i < xerr.size(); i++) {
            double w = 1.0/(xerr[i]*xerr[i]);
            sum += w;
        }
        return std::sqrt(1.0/(sum));
    }

    /**
     * @brief Calculate the mean of a container class.
     */
    template<container_type T>
    double mean(const T& v) noexcept {
        double sum = 0;
        for (unsigned int i = 0; i < v.size(); i++) {
            sum += v[i];
        }
        return sum/v.size();
    }

    /**
     * @brief Calculate the variance of a container class.
     */
    template<container_type T>
    double var(const T& v, unsigned int ddof = 1) noexcept {
        double mu = mean(v);
        double sum = 0;
        for (unsigned int i = 0; i < v.size(); i++) {
            sum += std::pow(v[i] - mu, 2);
        }
        return sum/(v.size() - ddof);
    }

    /**
     * @brief Calculate the standard deviation of a container class.
     */
    template<container_type T>
    double std(const T& v, unsigned int ddof = 1) noexcept {
        return std::sqrt(var(v, ddof));
    }

    /**
     * @brief Calculate the mode of a vector.
     */
    template<numeric T>
    T mode(const std::vector<T>& v) {
        if (v.empty()) {
            throw std::invalid_argument("stats::mode: Vector is empty.");
        }
        std::vector<T> v_copy = v;
        std::sort(v_copy.begin(), v_copy.end());
        T mode = v_copy[0];
        unsigned int max_count = 0;
        unsigned int count = 1;
        for (unsigned int i = 1; i < v_copy.size(); i++) {
            if (v_copy[i] == v_copy[i - 1]) {
                count++;
            } else {
                if (max_count < count) {
                    max_count = count;
                    mode = v_copy[i-1];
                }
            }
        }

        // handle edge-case where last sorted element is the mode
        if (max_count < count) {
            max_count = count;
            mode = v_copy.back();
        }
        return mode;
    }

    template<numeric T>
    struct Measurement {
        Measurement() {}
        Measurement(const std::vector<T>& vals) : vals(vals) {}

        double mean() const noexcept {return mean(vals);}
        double std() const noexcept {return std(vals);}
        double var() const noexcept {return var(vals);}

        std::vector<T> vals;
    };
}