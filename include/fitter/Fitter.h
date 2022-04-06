#pragma once

#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <utility>
#include <memory>

using std::cout, std::endl, std::string, std::shared_ptr;

class Fitter {
  public:
    struct Fit {
        Fit() {}

        Fit(std::map<string, double>& params, std::map<string, double>& errs, const double chi2, const int dof, const int calls, const bool converged) : 
            params(params), errs(errs), chi2(chi2), dof(dof), calls(calls), converged(converged) {}

        void print() const {
            cout << "\n############################################################" << endl;
            cout << (converged ? "Successfully " : "Failed to ") << "fit the data with " << calls << " function evaluations." << endl;
            cout << "\u03C7^2 = " << std::setprecision(6) << chi2 << ", dof = " << dof << ", \u03C7^2/dof = " << chi2/dof << endl;
            for (const auto& e : params) {
                cout << "\t" << std::setprecision(6) << e.first << " = "  << std::setw(10) << e.second << " ± " << std::setw(10) << errs.at(e.first) << endl;
            }
            cout << "############################################################\n" << endl;
        }

        std::map<string, double> params;
        std::map<string, double> errs;
        double chi2;
        double dof;
        double calls;
        bool converged;
    };

    virtual ~Fitter() {}
    virtual std::shared_ptr<Fit> fit() = 0;
};