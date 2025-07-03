// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <mini/MiniFwd.h>
#include <utility/Limit.h>
#include <dataset/PointSet.h>

#include <string>

namespace ausaxs::mini {
    struct FittedParameter {
        FittedParameter() noexcept = default;

        /**
         * @brief Create a FittedParameter with asymmetric errors.
         * 
         * @param name Name of this parameter. 
         * @param value Optimal value of this parameter.
         * @param error Asymmetrical errors of this parameter.
         */
        FittedParameter(const std::string& name, double value, const Limit& error) noexcept;

        /**
         * @brief Create a FittedParameter with symmetric errors.
         * 
         * @param name Name of this parameter. 
         * @param value Optimal value of this parameter.
         * @param error Symmetrical errors of this parameter.
         */
        FittedParameter(const std::string& name, double value, double error) noexcept;

        /**
         * @brief Create a FittedParameter from a Parameter with asymmetric errors.
         * 
         * @param param The fitted parameter.
         * @param value Optimal value of this parameter.
         * @param error Asymmetrical errors of this parameter.
         */
        FittedParameter(const Parameter& param, double value, const Limit& error) noexcept;

        /**
         * @brief Create a FittedParameter from a Parameter with symmetric errors.
         * 
         * @param param The fitted parameter.
         * @param value Optimal value of this parameter.
         * @param error Symmetrical errors of this parameter.
         */
        FittedParameter(const Parameter& param, double value, double error) noexcept;

        operator Point1D() const noexcept {return Point1D(value, mean_error());}

        /**
         * @brief Convenience method for implicit conversion to double. 
         */
        operator double() const noexcept {return value;}

        /**
         * @brief Get a string representation of this parameter.
         */
        std::string to_string() const noexcept;

        /**
         * @brief Get the mean error. 
         *        If the errors are asymmetric, this returns their mean. If not, the error is returned. 
         */
        double mean_error() const noexcept;

        /**
         * @brief Set a symmetric error. 
         */
        void set_error(double error) noexcept;

        /**
         * @brief Set an asymmetric error. 
         */
        void set_error(double min, double max) noexcept;

        /**
		 * @brief Stream output operator. 
		 * 
		 * Allows this object to easily be output to a given stream. 
		 */
		friend std::ostream& operator<<(std::ostream& os, const FittedParameter& param) noexcept {os << param.to_string(); return os;}

        std::string name; // The name of this parameter.
        double value;     // The optimal value of this parameter.
        Limit error;      // The error on this parameter. 
    };
}