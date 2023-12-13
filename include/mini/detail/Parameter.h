#pragma once

#include <mini/MiniFwd.h>
#include <utility/Limit.h>

#include <string>
#include <vector>
#include <optional>

namespace mini {
    /**
     * @brief A representation of a parameter.
     */
    struct Parameter {
        Parameter() noexcept;

        /**
         * @brief Create a Parameter with a guess value and bounds.
         * 
         * @param name The name of the parameter.
         * @param guess The guess value.
         * @param bounds The bounds. 
         */
        Parameter(const std::string& name, double guess = 0) noexcept;

        /**
         * @brief Create a Parameter with a guess value and bounds.
         * 
         * @param name The name of the parameter.
         * @param guess The guess value.
         * @param bounds The bounds. 
         */
        Parameter(const std::string& name, double guess, const Limit& bounds) noexcept;

        /**
         * @brief Create a Parameter without a guess value.
         * 
         * @param name The name of the parameter.
         * @param bounds The bounds.
         */
        Parameter(const std::string& name, const Limit& bounds) noexcept;

        /**
         * @brief Create a Parameter from a FittedParameter.
         */
        Parameter(const mini::FittedParameter& p) noexcept;

        ~Parameter();

        /**
         * @brief Check if this parameter is bounded.
         */
        [[nodiscard]]
        bool has_bounds() const noexcept;

        /**
         * @brief Check if this parameter has a guess value.
         */
        [[nodiscard]]
        bool has_guess() const noexcept;

        /**
         * @brief Check if this parameter is named.
         */
        [[nodiscard]]
        bool has_name() const noexcept;

        /**
         * @brief Check if this parameter has been initialized.
         */
        [[nodiscard]]
        bool empty() const noexcept;

        /**
         * @brief Set this equal to a fit parameter.
         */
        Parameter& operator=(const mini::FittedParameter& other) noexcept;

        /**
         * @brief Get a string representation of this parameter.
         */
        std::string to_string() const noexcept;

        /**
         * @brief Stream output operator. 
         * 
         * Allows this object to easily be output to a given stream. 
         */
        friend std::ostream& operator<<(std::ostream& os, const Parameter& param) noexcept {os << param.to_string(); return os;}

        std::string name;            // The name of this parameter.
        std::optional<double> guess; // The guess value.
        std::optional<Limit> bounds; // The bounds of this parameter. 
    };
}