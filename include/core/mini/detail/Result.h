#pragma once

#include <mini/detail/FittedParameter.h>

#include <string>
#include <vector>

namespace mini {
    struct Result {
        Result() noexcept = default;
        Result(const FittedParameter& param, double fval, unsigned int fevals) noexcept;
        Result(const std::vector<FittedParameter>& params, double fval, unsigned int fevals) noexcept;
        virtual ~Result() = default;

        /**
         * @brief Get a parameter based on its name from this result.
         */
        const FittedParameter& get_parameter(const std::string& name) const;

        /**
         * @brief Get a parameter based on its name from this result.
         */
        FittedParameter& get_parameter(const std::string& name);

        /**
         * @brief Get a parameter based on its index from this result.
         */
        const FittedParameter& get_parameter(unsigned int index) const;

        /**
         * @brief Get a parameter based on its index from this result.
         */
        FittedParameter& get_parameter(unsigned int index);

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