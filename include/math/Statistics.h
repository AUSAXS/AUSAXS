#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include <type_traits>
#include <initializer_list>
#include <math.h>
#include <iostream>
#include <math/ConstSlice.h>
#include <math/MutableSlice.h>

namespace stats {
    namespace detail {
        /**
         * @brief Calculate the weighted mean of a container class.
         *        Requires operator[] and size() to be defined.
         */
        template<typename T, typename Q>
        double weighted_mean(T& x, Q& xerr) noexcept {
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
         *        Requires operator[] and size() to be defined.
         */
        template<typename T>
        double weighted_mean_error(T& xerr) noexcept {
            double sum = 0;
            for (unsigned int i = 0; i < xerr.size(); i++) {
                double w = 1.0/(xerr[i]*xerr[i]);
                sum += w;
            }
            return std::sqrt(1.0/(sum));
        }

        /**
         * @brief Calculate the mean of a container class.
         *        Requires operator[] and size() to be defined.
         */
        template<typename T>
        double mean(T& v) noexcept {
            double sum = 0;
            for (unsigned int i = 0; i < v.size(); i++) {
                sum += v[i];
            }
            return sum/v.size();
        }

        /**
         * @brief Calculate the variance of a container class.
         *        Requires operator[] and size() to be defined.
         */
        template<typename T>
        double var(T& v, unsigned int ddof) noexcept {
            double mu = mean(v);
            double sum = 0;
            for (unsigned int i = 0; i < v.size(); i++) {
                sum += std::pow(v[i] - mu, 2);
            }
            return sum/(v.size() - ddof);
        }
    }

    /**
     * @brief Calculate the mode of a vector.
     */
    template<typename T>
    T mode(const std::vector<T>& v) {
        if (v.empty()) {
            throw std::invalid_argument("Error in stats::mode: Vector is empty.");
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

    /**
     * @brief Calculate the mean of a vector.
     */
    template<typename T>
    double weighted_mean_error(std::vector<T> xerr) noexcept {
        return detail::weighted_mean_error(xerr);
    }

    /**
     * @brief Calculate the mean of a vector.
     */
    template<typename T>
    double weighted_mean_error(ConstSlice<T> xerr) noexcept {
        return detail::weighted_mean_error(xerr);
    }

    /**
     * @brief Calculate the mean of a vector.
     */
    template<typename T>
    double weighted_mean_error(MutableSlice<T> xerr) noexcept {
        return detail::weighted_mean_error(xerr);
    }

    /**
     * @brief Calculate the weighted mean of a vector.
     */
    template<typename T, typename Q>
    double weighted_mean(std::vector<T> x, std::vector<Q> xerr) noexcept {
        return detail::weighted_mean(x, xerr);
    }

    /**
     * @brief Calculate the weighted mean of a Slice.
     */
    template<typename T, typename Q>
    double weighted_mean(ConstSlice<T> x, ConstSlice<Q> xerr) noexcept {
        return detail::weighted_mean(x, xerr);
    }

    /**
     * @brief Calculate the weighted mean of a Slice.
     */
    template<typename T, typename Q>
    double weighted_mean(MutableSlice<T> x, MutableSlice<Q> xerr) noexcept {
        return detail::weighted_mean(x, xerr);
    }

    /**
     * @brief Calculate the mean of a vector.
     */
    template<typename T>
    double mean(std::vector<T> v) noexcept {
        return detail::mean(v);
    }

    /**
     * @brief Calculate the mean of a Slice.
     */
    template<typename T>
    double mean(ConstSlice<T> s) noexcept {
        return detail::mean(s);
    }

    /**
     * @brief Calculate the mean of a Slice.
     */
    template<typename T>
    double mean(MutableSlice<T> s) noexcept {
        return detail::mean(s);
    }

    /**
     * @brief Calculate the mean of a list.
     */
    template<typename T>
    double mean(std::initializer_list<T> v) noexcept {
        return mean(std::vector<T>(v));
    }

    /**
     * @brief Calculate the variance of a Slice.
     */
    template<typename T>
    double var(ConstSlice<T> s, unsigned int ddof = 1) noexcept {
        return detail::var(s, ddof);
    }

    /**
     * @brief Calculate the variance of a Slice.
     */
    template<typename T>
    double var(MutableSlice<T> s, unsigned int ddof = 1) noexcept {
        return detail::var(s, ddof);
    }

    /**
     * @brief Calculate the variance of a vector.
     */
    template<typename T>
    double var(std::vector<T> v, unsigned int ddof = 1) noexcept {
        return detail::var(v, ddof);
    }

    /**
     * @brief Calculate the variance of a list.
     */
    template<typename T>
    double var(std::initializer_list<T> v, unsigned int ddof = 1) noexcept {
        return var(std::vector<T>(v), ddof);
    }

    /**
     * @brief Calculate the standard deviation of a Slice.
     */
    template<typename T>
    double std(ConstSlice<T> s, unsigned int ddof = 1) noexcept {
        return std::sqrt(var(s, ddof));
    }

    /**
     * @brief Calculate the standard deviation of a Slice.
     */
    template<typename T>
    double std(MutableSlice<T> s, unsigned int ddof = 1) noexcept {
        return std::sqrt(var(s, ddof));
    }

    /**
     * @brief Calculate the standard deviation of a vector.
     */
    template<typename T>
    double std(std::vector<T> v, unsigned int ddof = 1) noexcept {
        return std::sqrt(var(v, ddof));
    }

    /**
     * @brief Combine a list of errors to a single error of the mean. 
     */
    template<typename T>
    double std(std::initializer_list<T> v, unsigned int ddof = 1) noexcept {
        return std(std::vector<T>(v), ddof);
    }

    template<typename T>
    struct Measurement {
        Measurement() {}
        Measurement(std::vector<T> vals) : vals(vals) {}

        double mean() const noexcept {return mean(vals);}
        double std() const noexcept {return std(vals);}
        double var() const noexcept {return var(vals);}

        std::vector<T> vals;
    };
}