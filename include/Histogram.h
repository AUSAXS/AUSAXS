#pragma once

#include <vector>
#include <string>
#include <utility>
#include <memory>

using std::vector, std::string, std::shared_ptr, std::unique_ptr;

class Histogram {
    public:
        Histogram() {}
        Histogram(const vector<double>& p) : p(p) {}
        Histogram(const vector<double>& p, const vector<int> axes) : p(p), axes(axes) {}

        vector<double> p;
        vector<int> axes;
};