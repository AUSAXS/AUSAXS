#pragma once

#include <vector>
#include <string>
#include <utility>
#include <memory>

#include "data/Axis.h"

using std::vector, std::string, std::shared_ptr, std::unique_ptr;

class Histogram {
  public:
    Histogram() {}
    Histogram(const vector<double>& p) : p(p) {}
    Histogram(const vector<double>& p, const Axis& axis) : p(p), axis(axis) {}

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

    vector<double> p;
    Axis axis;
};