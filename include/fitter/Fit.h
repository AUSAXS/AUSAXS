#pragma once

#include <string>
#include <vector>

#include <utility/Dataset.h>
#include <utility/Multiset.h>
#include <minimizer/Utility.h>

class Fitter;

class Fit : public mini::Result {
    public:
        struct Plots {
            SimpleDataset intensity;              // The full intensity line
            SimpleDataset intensity_interpolated; // The intensity line interpolated at the data points. 
            SimpleDataset data;                   // The data itself
        };

        Fit() noexcept {}

        /**
         * @brief Constructor.
         * 
         * Create a new Fit object based on a fitter and a minimizer result.
         */
        Fit(Fitter& fitter, const mini::Result& res, double chi2) noexcept;

        /**
         * @brief Constructor.
         * 
         * Create a new Fit object based on a minimizer result.
         */
        Fit(const mini::Result& res, double chi2, double dof) noexcept;
        
        /**
         * @brief Add the parameters from another fit to this one. This is useful to add the parameters from an inner fit to an outer one. 
         */
        void add_fit(Fitter& fit) noexcept;

        /**
         * @brief Add the parameters from another fit to this one. This is useful to add the parameters from an inner fit to an outer one. 
         */
        void add_fit(std::shared_ptr<Fit> fit) noexcept;

        /**
         * @brief Get a string representation of this object. 
         */
        std::string to_string() const noexcept;

        SimpleDataset evaluated_points;
        Plots figures;
        SimpleDataset residuals;
        unsigned int dof;
};