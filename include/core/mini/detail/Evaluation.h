#pragma once

#include <vector>
#include <string>

namespace mini {
    struct Evaluation {
        Evaluation() noexcept = default;

        Evaluation(std::vector<double> vals, double fval) : vals(std::move(vals)), fval(fval) {}

        std::string to_string() const;

        bool operator==(const Evaluation& other) const = default;

        std::vector<double> vals;
        double fval = 0;
    };
}