#pragma once

#include <string>
#include <vector>

#include <fitter/Fitter.h>
#include <utility/Dataset.h>
#include <utility/Multiset.h>

#include <Math/Minimizer.h>

class Fit {
    public:
        Fit() {}

        /**
         * @brief Smart constructor.
         * 
         * Construct a new Fit object based on a Fitter and Minimizer. 
         */
        Fit(Fitter& fitter, const ROOT::Math::Minimizer* const minimizer, double chi2);

        /**
         * @brief Constructor.
         * 
         * Construct a new Fit object. 
         */
        Fit(std::map<std::string, double>& params, std::map<std::string, double>& errs, const double chi2, const int dof, const int calls, const bool converged);

        /**
         * @brief Get a string representation of this object. 
         */
        std::string to_string() const;

        Dataset evaluated_points;
        Multiset figures;
        Dataset residuals;
        std::map<std::string, double> params;
        std::map<std::string, double> errors;
        double chi2;
        unsigned int dof;
        unsigned int calls;
        bool converged;
};