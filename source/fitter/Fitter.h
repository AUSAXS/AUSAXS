#pragma once

#include <string>
#include <vector>
#include <map>

class Fitter {
public:
    struct Fit {
        Fit(std::map<string, double>& params) : params(params) {}
        std::map<string, double> params;
    };

    virtual ~Fitter() {}
    virtual Fit fit() const = 0;
};