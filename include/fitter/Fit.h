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
         * Construct a new Fit object based on a Fitter and Minimizer. This is meant to be used when the minimizer optimizes some parameter 
         * changing the histograms being fitted by the fitter.
         */
        Fit(Fitter& fitter, const ROOT::Math::Minimizer* const minimizer, double chi2);

        /**
         * @brief Constructor.
         * 
         * Create a new Fit object based on a minimizer. If no chi2 value is provided, it is assumed to be equal to the function minimum. 
         * 
         * @param minimizer 
         * @param chi2 
         */
        Fit(const ROOT::Math::Minimizer* const minimizer, double chi2 = -1);

        /**
         * @brief Constructor.
         * 
         * Construct a new Fit object. 
         */
        Fit(std::map<std::string, double>& params, std::map<std::string, double>& errs, const double chi2, const int dof, const int calls, const bool converged);

        /**
         * @brief Add the parameters from another fit to this one. This is useful to add the parameters from an inner fit to an outer one. 
         */
        void add_fit(Fitter& fit);

        /**
         * @brief Add the parameters from another fit to this one. This is useful to add the parameters from an inner fit to an outer one. 
         */
        void add_fit(std::shared_ptr<Fit> fit);

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