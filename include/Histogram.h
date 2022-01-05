#pragma once

#include <vector>
#include <string>
#include <utility>
#include <memory>

#include "data/Axes.h"

using std::vector, std::string, std::shared_ptr, std::unique_ptr;

class Histogram {
    public:
        Histogram() {}
        Histogram(const vector<double>& p) : p(p) {}
        Histogram(const vector<double>& p, const Axes axes) : p(p), axes(axes) {}

        /**
         * @brief Reduce the view axis to show only the non-zero area. 
         *        Minimum size is 10 units.
         */
        void shorten_axes() {
            int max_bin = 10; // minimum size is 10
            for (int i = axes.bins-1; i >= 10; i--) {
                if (p[i] != 0) {
                    max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
                    break;
                }
            }
            p.resize(max_bin);
            double width = (double) (axes.xmax-axes.xmin)/axes.bins;
            axes = {max_bin, 0, int(max_bin*width)};
        }

        vector<double> p;
        Axes axes;
};