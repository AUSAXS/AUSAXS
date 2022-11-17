#pragma once

#include <string>
#include <vector>

#include <mini/utility/FittedParameter.h>

namespace mini {
    struct Result {
        /**
         * @brief Default constructor.
         */
        Result() noexcept {}

        /**
         * @brief Construct a Result. 
         * 
         * @param params The fitted parameters.
         * @param fval The function value.
         * @param fevals The number of function evaluations.
         */
        Result(const FittedParameter& param, double fval, unsigned int fevals) noexcept;

        /**
         * @brief Construct a Result. 
         * 
         * @param params The fitted parameters.
         * @param fval The function value.
         * @param fevals The number of function evaluations.
         */
        Result(const std::vector<FittedParameter>& params, double fval, unsigned int fevals) noexcept;

        /**
         * @brief Get a parameter based on its name from this result.
         */
        FittedParameter const& get_parameter(std::string name) const;

        /**
         * @brief Get a parameter based on its name from this result.
         */
        FittedParameter& get_parameter(std::string name);

        /**
         * @brief Get a parameter based on its index from this result.
         */
        FittedParameter const& get_parameter(unsigned int index) const;

        /**
         * @brief Get a parameter based on its index from this result.
         */
        FittedParameter& get_parameter(unsigned int index);

        /**
         * @brief Add a parameter to this result.
         */
        void add_parameter(const FittedParameter& param) noexcept;

        /**
         * @brief Get the number of parameters in this result.
         */
        size_t size() const noexcept;

        /**
         * @brief Get the number of parameters in this result.
         */
        size_t dim() const noexcept;

        std::vector<FittedParameter> parameters; // The fitted parameters
        double fval;                             // The minimum function value
        unsigned int fevals;                     // The number of function evaluations
        int status = 0;                          // The minimization status. Anything but 0 indicates an error with the fit. 
    };
}