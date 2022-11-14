#pragma once

#include <string>
#include <vector>

#include <utility/Dataset.h>
#include <utility/Multiset.h>
#include <mini/Utility.h>

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
         * @brief Add the parameters from another fit to this one. Each parameter will count as an additional degree of freedom. 
         */ 
        void add_fit(Fitter& fit) noexcept;

        /**
         * @brief Add the parameters from another fit to this one. Each parameter will count as an additional degree of freedom. 
         */
        void add_fit(std::shared_ptr<Fit> fit) noexcept;

        /**
         * @brief Add plots to this fit.
         */
        void add_plots(Fitter& fitter);

        /**
         * @brief Get a string representation of this object. 
         */
        virtual std::string to_string() const noexcept;

        SimpleDataset evaluated_points;
        Plots figures;
        SimpleDataset residuals;
        unsigned int dof;
};

struct EMFit : public Fit {
    using Fit::Fit;

    std::string to_string() const noexcept override;

    double level;
};

template<typename C>
concept FitType = std::derived_from<C, Fit>;