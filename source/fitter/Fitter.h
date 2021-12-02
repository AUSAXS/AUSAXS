#pragma once

#include <string>
#include <vector>
#include <map>

class Fitter {
public:
    struct Fit {
        Fit() {}
        Fit(std::map<string, double>& params, const double& chi2, const int& dof) : params(params), chi2(chi2), dof(dof) {}
        std::map<string, double> params;
        double chi2;
        double dof;
    };

    virtual ~Fitter() {}
    virtual Fit fit() = 0;
};