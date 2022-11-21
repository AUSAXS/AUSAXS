#pragma once

#include <vector>
#include <string>

namespace mini {
    struct Evaluation {
        Evaluation() = default;

        Evaluation(std::vector<double> vals, double fval) : vals(vals), fval(fval) {}

        std::string to_string() const;

        std::vector<double> vals;
        double fval = 0;
    };
}