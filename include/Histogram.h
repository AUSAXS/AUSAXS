#pragma once

#include <vector>
#include <string>
#include <utility>
#include <memory>

#include "data/Axis.h"

using std::vector, std::string, std::shared_ptr, std::unique_ptr;

/**
 * @brief \class Histogram. 
 * 
 * A representation of a histogram. 
 */
class Histogram {
  public:
    /**
     * @brief Default constructor.
     */
    Histogram() {}

    /**
     * @brief Constructor. 
     * 
     * Construct a new histogram based on a list of bin values. 
     * Note that the axis will not be initialized. 
     * 
     * @param p The bin values. 
     */
    Histogram(const vector<double>& p) : p(p) {}

    /**
     * @brief Constructor.
     * 
     * Construct a new histogram based on a list of bin values and the axis it spans. 
     * 
     * @param p The bin values. 
     * @param axis The axis they span. 
     */
    Histogram(const vector<double>& p, const Axis& axis) : p(p), axis(axis) {}

    /**
     * @brief Add another Histogram to this one.
     */
    Histogram& operator+=(const Histogram& rhs) {
        std::transform(p.begin(), p.end(), rhs.p.begin(), p.begin(), std::plus<double>());
        return *this;
    }

    /**
     * @brief Subtract another Histogram from this one.
     */
    Histogram& operator-=(const Histogram& rhs) {
        std::transform(p.begin(), p.end(), rhs.p.begin(), p.begin(), std::minus<double>());
        return *this;
    }

    /**
     * @brief Reduce the view axis to show only the non-zero area. 
     *        Minimum size is 10 units.
     */
    void shorten_axis() {
        int max_bin = 10; // minimum size is 10
        for (int i = axis.bins-1; i >= 10; i--) {
            if (p[i] != 0) {
                max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
                break;
            }
        }
        p.resize(max_bin);
        double width = axis.width();
        axis = Axis{max_bin, 0, max_bin*width};
    }

    vector<double> p; // The bin values. 
    Axis axis;        // The axis spanned by this histogram. 
};