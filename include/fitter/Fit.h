#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <fitter/Fitter.h>

class Fit {
    public:
        Fit() {}

        Fit(Fitter& fitter, const ROOT::Math::Minimizer* const minimizer, double chi2) : chi2(chi2) {
            unsigned int vars = minimizer->NDim();
            const double* result = minimizer->X();
            const double* errs = minimizer->Errors();
            for (unsigned int i = 0; i < vars; i++) {
                std::cout << "VARIABLE " << i << ": " << minimizer->VariableName(i) << std::endl;
                params.insert({minimizer->VariableName(i), result[i]});
                errors.insert({minimizer->VariableName(i), errs[i]});
            }

            normal_plot = fitter.plot();
            residual_plot = fitter.plot_residuals();

            converged = minimizer->Status() == 0;
            calls = minimizer->NCalls();
            dof = fitter.dof() - vars;
        }

        Fit(std::map<std::string, double>& params, std::map<std::string, double>& errs, const double chi2, const int dof, const int calls, const bool converged) : 
            params(params), errors(errs), chi2(chi2), dof(dof), calls(calls), converged(converged) {}

        void print() const {
            std::cout << "\n############################################################" << std::endl;
            std::cout << (converged ? "Successfully " : "Failed to ") << "fit the data with " << calls << " function evaluations." << std::endl;
            std::cout << "\u03C7^2 = " << std::setprecision(6) << chi2 << ", dof = " << dof << ", \u03C7^2/dof = " << chi2/dof << std::endl;
            for (const auto& e : params) {
                std::cout << "\t" << std::setprecision(6) << e.first << " = "  << std::setw(10) << e.second << " ± " << std::setw(10) << errors.at(e.first) << std::endl;
            }
            std::cout << "############################################################\n" << std::endl;
        }

        std::vector<std::shared_ptr<TGraph>> normal_plot;
        std::shared_ptr<TGraph> residual_plot;
        std::map<std::string, double> params;
        std::map<std::string, double> errors;
        double chi2;
        unsigned int dof;
        unsigned int calls;
        bool converged;
};