#pragma once

#include <utility/Utility.h>
#include <utility/PointSet.h>
#include <utility/Exceptions.h>

#include <vector>
#include <algorithm>
#include <numeric>
#include <type_traits>
#include <initializer_list>
#include <math.h>
#include <iostream>

namespace stats {
    /**
     * @brief Calculate the mode of a vector.
     */
    template<typename T, typename = typename std::enable_if<std::is_arithmetic_v<T>>::type>
    T mode(const std::vector<T>& v) {
        if (v.empty()) {
            throw except::invalid_operation("Error in stats::mode: Vector is empty.");
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
    template<typename T, typename = typename std::enable_if<std::is_arithmetic_v<T>, T>::type>
    double mean(std::vector<T> v) noexcept {
        return std::accumulate(v.begin(), v.end(), 0.0, [] (double sum, T x) {return sum + x;})/v.size();
    }

    /**
     * @brief Calculate the mean of a list.
     */
    template<typename T, typename = typename std::enable_if<std::is_arithmetic_v<T>, T>::type>
    double mean(std::initializer_list<T> v) noexcept {
        return mean(std::vector<T>(v));
    }

    /**
     * @brief Calculate the variance of a vector.
     */
    template<typename T, typename = typename std::enable_if<std::is_arithmetic_v<T>, T>::type>
    double var(std::vector<T> v, double ddof = 1) noexcept {
        double mu = mean(v);
        double sum = std::accumulate(v.begin(), v.end(), 0.0, [mu] (double sum, T x) {return sum + std::pow(x - mu, 2);});
        return sum/(v.size() - ddof);
    }

    /**
     * @brief Calculate the variance of a list.
     */
    template<typename T, typename = typename std::enable_if<std::is_arithmetic_v<T>, T>::type>
    double var(std::initializer_list<T> v, double ddof = 1) noexcept {
        return var(std::vector<T>(v), ddof);
    }


    /**
     * @brief Calculate the standard deviation of a vector.
     */
    template<typename T, typename = typename std::enable_if<std::is_arithmetic_v<T>, T>::type>
    double std(std::vector<T> v, double ddof = 1) noexcept {
        return std::sqrt(var(v, ddof));
    }

    /**
     * @brief Combine a list of errors to a single error of the mean. 
     */
    template<typename T, typename = typename std::enable_if<std::is_arithmetic_v<T>, T>::type>
    double std(std::initializer_list<T> v, double ddof = 1) noexcept {
        return std(std::vector<T>(v), ddof);
    }

    template<typename T, typename = typename std::enable_if<std::is_arithmetic_v<T>, T>::type>
    struct Measurement {
        Measurement() {}
        Measurement(std::vector<T> vals) : vals(vals) {}

        double mean() const noexcept {return mean(vals);}
        double std() const noexcept {return std(vals);}
        double var() const noexcept {return var(vals);}

        std::vector<T> vals;
    };

    struct MeasurementSeries {
        MeasurementSeries() {}
        MeasurementSeries(unsigned int size) : data(size) {}
        MeasurementSeries(std::vector<double> means, std::vector<double> stds) : data(means.size()) {
            if (means.size() != stds.size()) {throw except::size_error("Error in MeasurementSeries::MeasurementSeries: The two vectors \"means\" and \"stds\" must be of same length.");}
            std::transform(means.begin(), means.end(), stds.begin(), data.begin(), [] (double mean, double std) {return Point1D(mean, std);});
        }

        double weighted_mean() const noexcept {
            double num = std::accumulate(data.begin(), data.end(), 0.0, [] (double sum, const Point1D m) {return sum + m.x*m.xerr;});
            double denom = std::accumulate(data.begin(), data.end(), 0.0, [] (double sum, const Point1D m) {return sum + m.xerr;});
            return num/denom;
        }

        double error_of_mean() const noexcept {
            double dev = std(get_means());
            return dev/std::sqrt(data.size());
        }

        std::vector<double> get_means() const noexcept {
            std::vector<double> v(data.size());
            std::transform(data.begin(), data.end(), v.begin(), [] (const Point1D m) {return m.x;});
            return v;
        }

        std::vector<double> get_stds() const noexcept {
            std::vector<double> v(data.size());
            std::transform(data.begin(), data.end(), v.begin(), [] (const Point1D m) {return m.xerr;});
            return v;
        }

        void push_back(const Point1D& measurement) noexcept {data.push_back(measurement);}
        void push_back(double mean, double std) noexcept {data.push_back({mean, std});}

        
		const std::vector<Point1D>::const_iterator begin() const {return data.begin();}
        const std::vector<Point1D>::const_iterator end() const {return data.end();}
        std::vector<Point1D>::iterator begin() {return data.begin();}
        std::vector<Point1D>::iterator end() {return data.end();}

        std::vector<Point1D> data;
    };
}