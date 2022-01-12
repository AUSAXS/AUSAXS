#pragma once

#include <vector>
#include <array>
#include <initializer_list>

struct Axes {
    public: 
        Axes() {}
        Axes(const std::vector<double>& axes) : bins(axes[0]), xmin(axes[1]), xmax(axes[2]) {}
        Axes(const double& bins, const double& xmin, const double& xmax) : bins(bins), xmin(xmin), xmax(xmax)  {}

        Axes& operator=(std::initializer_list<int> list) {
            std::vector<int> d = list;
            bins = d[0]; 
            xmin = d[1];
            xmax = d[2];
            return *this;
        }

        Axes& operator=(std::initializer_list<double> list) {
            std::vector<double> d = list;
            bins = d[0]; 
            xmin = d[1];
            xmax = d[2];
            return *this;
        }

        size_t bins;
        double xmin, xmax;
};