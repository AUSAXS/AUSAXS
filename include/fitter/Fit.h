#pragma once

#include <string>
#include <vector>

#include <utility/Dataset.h>
#include <utility/Multiset.h>
#include <minimizer/Utility.h>

class Fitter;

class Fit {
    public:
        struct Plots {
            Dataset intensity;              // The full intensity line
            Dataset intensity_interpolated; // The intensity line interpolated at the data points. 
            Dataset data;                   // The data itself
        };

        Fit() {}

        /**
         * @brief Constructor.
         * 
         * Create a new Fit object based on a fitter and a minimizer result.
         */
        Fit(Fitter& fitter, const mini::Result& res, double chi2);

        /**
         * @brief Constructor.
         * 
         * Create a new Fit object based on a minimizer result.
         */
        Fit(const mini::Result& res, double chi2, double dof);
        
        /**
         * @brief Add the parameters from another fit to this one. This is useful to add the parameters from an inner fit to an outer one. 
         */
        void add_fit(Fitter& fit);

        /**
         * @brief Add the parameters from another fit to this one. This is useful to add the parameters from an inner fit to an outer one. 
         */
        void add_fit(std::shared_ptr<Fit> fit);

        /**
         * @brief Get a parameter by index.
         */
        mini::FittedParameter get_parameter(unsigned int index) const;

        /**
         * @brief Get a parameter by name.
         */
        mini::FittedParameter get_parameter(std::string name) const;

        /**
         * @brief Get a string representation of this object. 
         */
        std::string to_string() const;

        Dataset evaluated_points;
        Plots figures;
        Dataset residuals;
        std::vector<mini::FittedParameter> params;
        double chi2;
        unsigned int dof;
        unsigned int calls;
        bool converged;
};