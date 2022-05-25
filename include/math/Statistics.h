#pragma once

#include <utility/Utility.h>

#include <vector>
#include <algorithm>
#include <numeric>
#include <type_traits>
#include <initializer_list>

#include <iostream>
namespace stats {
    /**
     * @brief Calculate the mean of a vector.
     */
    // template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
    template<typename T>
    double mean(std::vector<T> v) noexcept {
        return std::accumulate(v.begin(), v.end(), 0.0, [] (double sum, T x) {return sum + x;})/v.size();
    }

    /**
     * @brief Calculate the mean of a list.
     */
    // template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
    template<typename T>
    double mean(std::initializer_list<T> v) noexcept {
        return mean(std::vector<T>(v));
    }

    /**
     * @brief Calculate the variance of a vector.
     */
    // template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
    template<typename T>
    double var(std::vector<T> v, double ddof = 1) noexcept {
        double mu = mean(v);
        double sum = std::accumulate(v.begin(), v.end(), 0.0, [mu] (double sum, T x) {return sum + std::pow(x - mu, 2);});
        return sum/(v.size() - ddof);
    }

    /**
     * @brief Calculate the variance of a list.
     */
    // template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
    template<typename T>
    double var(std::initializer_list<T> v, double ddof = 1) noexcept {
        return var(std::vector<T>(v), ddof);
    }


    /**
     * @brief Calculate the standard deviation of a vector.
     */
    // template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
    template<typename T>
    double std(std::vector<T> v, double ddof = 1) noexcept {
        return std::sqrt(var(v, ddof));
    }

    /**
     * @brief Calculate the standard deviation of a list.
     */
    // template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
    template<typename T>
    double std(std::initializer_list<T> v, double ddof = 1) noexcept {
        return std(std::vector<T>(v), ddof);
    }
}