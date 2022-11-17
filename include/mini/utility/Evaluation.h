#pragma once

#include <vector>

namespace mini {
    struct Evaluation {
        Evaluation(std::vector<double> vals, double fval) : vals(vals), fval(fval) {}

        std::vector<double> vals;
        double fval;
    };
}