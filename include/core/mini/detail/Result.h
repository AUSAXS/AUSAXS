// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsFitParameters.h>
#include <mini/detail/FittedParameter.h>

#include <string>
#include <vector>

namespace ausaxs::mini {
    struct Result {
        Result() noexcept = default;
        Result(const FittedParameter& param, double fval, unsigned int fevals) noexcept;
        Result(const std::vector<FittedParameter>& params, double fval, unsigned int fevals) noexcept;
        virtual ~Result() = default;

        /**
         * @brief Get a parameter from this result.
         */
        const FittedParameter& get_parameter(const std::string& name) const;
        FittedParameter& get_parameter(const std::string& name); //< @copydoc get_parameter
        const FittedParameter& get_parameter(ausaxs::constants::fit::Parameters param) const; //< @copydoc get_parameter
        FittedParameter& get_parameter(ausaxs::constants::fit::Parameters param); //< @copydoc get_parameter
        const FittedParameter& get_parameter(unsigned int index) const; //< @copydoc get_parameter
        FittedParameter& get_parameter(unsigned int index); //< @copydoc get_parameter

        /**
         * @brief Get all parameters from this result.
         */
        const std::vector<FittedParameter>& get_parameters() const;
        std::vector<FittedParameter> get_parameters(); //< @copydoc get_parameters

        /**
         * @brief Get all parameter values from this result.
         */
        std::vector<double> get_parameter_values() const;

        /**
         * @brief Get a parameter based on its index from this result.
         */
        FittedParameter& operator[](unsigned int index);

        /**
         * @brief Get a parameter based on its index from this result.
         */
        const FittedParameter& operator[](unsigned int index) const;

        /**
         * @brief Add a parameter to this result.
         */
        void add_parameter(const FittedParameter& param) noexcept;

        /**
         * @brief Get the number of parameters in this result.
         */
        unsigned int size() const noexcept;

        /**
         * @brief Get the number of parameters in this result.
         */
        unsigned int dim() const noexcept;

        std::vector<FittedParameter> parameters; // The fitted parameters
        double fval;                             // The minimum function value
        unsigned int fevals;                     // The number of function evaluations
        int status = 0;                          // The minimization status. Anything but 0 indicates an error with the fit. 
    };
}